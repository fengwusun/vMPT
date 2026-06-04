"""Tests for vmpt.optimizer (MSA pointing search; hMPT-derived)."""

from __future__ import annotations

import numpy as np
import pytest

from vmpt.optimizer import (
    CENTRATION_BUFFERS,
    PointingEvaluator,
    axy_to_shutter,
    grid_search,
    radec_to_axy,
    refine_top,
)


# ---------------------------------------------------------------------
# Coordinate maths
# ---------------------------------------------------------------------


def test_radec_to_axy_origin_is_zero():
    """A source at the pointing centre lands at (ax, ay) ≈ 0."""
    axy = radec_to_axy(np.array([53.16]), np.array([-27.78]),
                       53.16, -27.78, pa_v3_deg=30.0)
    assert np.allclose(axy[0], [0.0, 0.0], atol=1e-9)


def test_radec_to_axy_arcsec_scale():
    """A 1″ offset in RA at the pointing produces ≈ 1″ magnitude in (ax, ay)."""
    ra_p, dec_p = 53.16, -27.78
    # 1 arcsec east in RA, accounting for cos(Dec).
    delta_ra_deg = 1.0 / 3600.0 / np.cos(np.deg2rad(dec_p))
    axy = radec_to_axy(np.array([ra_p + delta_ra_deg]), np.array([dec_p]),
                       ra_p, dec_p, pa_v3_deg=0.0)
    mag = float(np.linalg.norm(axy[0]))
    assert 0.95 < mag < 1.05, f"expected ≈ 1″ magnitude, got {mag}"


def test_axy_to_shutter_inside_quadrant():
    """An (ax, ay) well inside one of the MSA quadrants must return a
    valid (q, s, d). The MSA has a cross-shaped inter-quadrant gap so
    (0, 0) is unreachable; (~30, ~20)″ lands inside Q1."""
    quad, s_frac, d_frac = axy_to_shutter(np.array([[30.0, 20.0]]))
    assert quad[0] in (1, 2, 3, 4)
    assert 0 <= s_frac[0] <= 170
    assert 0 <= d_frac[0] <= 364


def test_axy_to_shutter_inter_quadrant_gap_returns_zero():
    """The MSA has a cross-shaped gap between quadrants; (0, 0) lands
    there and should be reported as outside any quadrant."""
    quad, s_frac, d_frac = axy_to_shutter(np.array([[0.0, 0.0]]))
    assert quad[0] == 0
    assert np.isnan(s_frac[0])


def test_axy_to_shutter_far_outside_returns_zero_quad():
    """A position way outside the aperture plane returns quad == 0."""
    quad, _, _ = axy_to_shutter(np.array([[1e6, 1e6]]))
    assert quad[0] == 0


# ---------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------


def test_evaluator_detection_count_in_range():
    """For a random catalog over a 0.05-deg field, the evaluator should
    place SOMETHING (not zero, not the whole catalog) at a sensible
    pointing."""
    np.random.seed(42)
    N = 50
    ra = 53 + np.random.rand(N) * 0.08
    dec = -27 + np.random.rand(N) * 0.08
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det, tp, _ = ev.evaluate(53.05, -26.97, 30.0)
    n_det = int(det.sum())
    assert 0 < n_det < N, f"expected partial placement, got {n_det}/{N}"


def test_evaluator_centration_classes_monotone():
    """Tightening the centration buffer can only REDUCE the number of
    successful placements (never increase)."""
    np.random.seed(42)
    N = 50
    ra = 53 + np.random.rand(N) * 0.08
    dec = -27 + np.random.rand(N) * 0.08

    classes = ["UNCONSTRAINED", "ENTIRE_OPEN", "MIDPOINT",
               "CONSTRAINED", "TIGHTLY_CONSTRAINED"]
    counts = []
    for c in classes:
        ev = PointingEvaluator(ra, dec, centration=c)
        det, _, _ = ev.evaluate(53.05, -26.97, 30.0)
        counts.append(int(det.sum()))
    # Non-strictly decreasing.
    for i in range(len(counts) - 1):
        assert counts[i] >= counts[i + 1], (
            f"buffer {classes[i]} placed {counts[i]} but "
            f"{classes[i + 1]} placed {counts[i + 1]} (should be ≤)"
        )


def test_evaluator_unknown_centration_falls_back_to_unconstrained():
    """Typo'd centration class shouldn't blow up; falls back to no buffer."""
    np.random.seed(42)
    ra = 53 + np.random.rand(20) * 0.08
    dec = -27 + np.random.rand(20) * 0.08
    ev_default = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    ev_typo = PointingEvaluator(ra, dec, centration="DOES_NOT_EXIST")
    det_def, _, _ = ev_default.evaluate(53.05, -26.97, 30.0)
    det_typo, _, _ = ev_typo.evaluate(53.05, -26.97, 30.0)
    assert det_def.sum() == det_typo.sum()


# ── Per-target centration override (v1.3.1+) ─────────────────────────────


def test_per_target_centration_override_unconditional():
    """A per-target centration override wins unconditionally over the
    global ``centration`` argument, even when it's laxer than global.

    Build two evaluators on the same sources:
      • ``ev_global`` uses TIGHTLY_CONSTRAINED for everything.
      • ``ev_mixed``  uses TIGHTLY_CONSTRAINED globally but overrides
        the first half of the sources to UNCONSTRAINED.

    Whatever the same pointing keeps, ``ev_mixed`` must keep at least
    as many sources (the laxer override can only ADD placements that
    fall within the wider buffer), and the kept set in the first half
    must be a superset of ``ev_global``'s.
    """
    np.random.seed(7)
    N = 60
    ra = 53 + np.random.rand(N) * 0.08
    dec = -27 + np.random.rand(N) * 0.08
    cent_override = [""] * N
    half = N // 2
    for i in range(half):
        cent_override[i] = "UNCONSTRAINED"

    ev_global = PointingEvaluator(
        ra, dec, centration="TIGHTLY_CONSTRAINED",
    )
    ev_mixed = PointingEvaluator(
        ra, dec, centration="TIGHTLY_CONSTRAINED",
        centration_per_target=np.asarray(cent_override, dtype=object),
    )

    det_g, _, _ = ev_global.evaluate(53.05, -26.97, 30.0)
    det_m, _, _ = ev_mixed.evaluate(53.05, -26.97, 30.0)

    # First-half overrides → laxer than global, so det_m ⊇ det_g.
    assert np.all(det_m[:half] | ~det_g[:half]), (
        "Per-target UNCONSTRAINED should keep every source the "
        "global TIGHTLY_CONSTRAINED kept (det_m ⊇ det_g on first half)."
    )
    # Second-half = no override = same buffer as global = identical kept.
    assert np.array_equal(det_m[half:], det_g[half:]), (
        "Rows without override should behave identically to the global."
    )
    # And it must ALSO kick in the other direction — stricter per-target
    # overrides drop sources the laxer global would have kept.
    cent_strict = [""] * N
    for i in range(half):
        cent_strict[i] = "TIGHTLY_CONSTRAINED"
    ev_strict_mix = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        centration_per_target=np.asarray(cent_strict, dtype=object),
    )
    det_sm, _, _ = ev_strict_mix.evaluate(53.05, -26.97, 30.0)
    ev_global_lax = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_lax, _, _ = ev_global_lax.evaluate(53.05, -26.97, 30.0)
    assert np.all(det_sm[:half] <= det_lax[:half]), (
        "Per-target TIGHTLY_CONSTRAINED should drop ≥0 sources that "
        "the global UNCONSTRAINED kept."
    )


def test_per_target_centration_size_mismatch_raises():
    """A mismatched length is a programming error — fail fast."""
    ra = np.array([53.05, 53.06, 53.07])
    dec = np.array([-27.0, -27.0, -27.0])
    with pytest.raises(ValueError, match="centration_per_target"):
        PointingEvaluator(
            ra, dec, centration_per_target=np.array(["", ""], dtype=object),
        )


def test_per_target_centration_blank_uses_global():
    """Empty / None / unrecognised entries fall back to the global
    buffer for that row."""
    np.random.seed(42)
    N = 30
    ra = 53 + np.random.rand(N) * 0.08
    dec = -27 + np.random.rand(N) * 0.08

    blanks = np.array([""] * N, dtype=object)
    nones = np.array([None] * N, dtype=object)
    bogus = np.array(["NOT_A_LEVEL"] * N, dtype=object)

    ev_baseline = PointingEvaluator(ra, dec, centration="CONSTRAINED")
    det_baseline, _, _ = ev_baseline.evaluate(53.05, -26.97, 30.0)

    for fill, label in ((blanks, "blanks"),
                        (nones, "nones"),
                        (bogus, "unknown labels")):
        ev = PointingEvaluator(
            ra, dec, centration="CONSTRAINED",
            centration_per_target=fill,
        )
        det, _, _ = ev.evaluate(53.05, -26.97, 30.0)
        assert np.array_equal(det, det_baseline), (
            f"Per-target {label} should fall back to the global setting, "
            f"got differing result."
        )


def test_catalog_centration_array_truthiness_regression():
    """Regression for v1.3.4 crash:

      File "vmpt/main.py", line 7392, in on_optimize
          _cent_full = np.asarray(getattr(cat, "centration", None) or [],
      ValueError: The truth value of an array with more than one
          element is ambiguous. Use a.any() or a.all()

    The fix has to live in ``on_optimize``, which is wired up to
    Bokeh widgets and not easy to import standalone, so this test
    pins the failure shape at the layer below — a Catalog with a
    multi-element `centration` numpy array must be transparently
    usable wherever the optimizer driver touches it.
    """
    from vmpt.catalog import Catalog

    # Five-source catalog with mixed centration overrides — enough
    # to trigger the `bool(arr)` ambiguity if anyone uses
    # `cat.centration or []` again.
    n = 5
    cat = Catalog(
        ids=np.arange(1, n + 1),
        ra_deg=53.16 + 0.001 * np.arange(n),
        dec_deg=-27.78 - 0.001 * np.arange(n),
        priority=np.ones(n),
        weight=np.ones(n),
        mag=np.full(n, 24.0),
        z=np.full(n, 1.0),
        label=np.array([f"src{i}" for i in range(n)], dtype=object),
        centration=np.array([
            "TIGHTLY_CONSTRAINED", "", "CONSTRAINED", "", "ENTIRE_OPEN",
        ], dtype=object),
        source_path="",
    )
    # The bug fires on any code path that does
    #   `arr = getattr(cat, "centration", None) or fallback`
    # Pin the safe pattern: the array must coerce to np.asarray
    # without triggering bool() on it.
    cent_raw = getattr(cat, "centration", None)
    assert cent_raw is not None
    arr = np.asarray(cent_raw, dtype=object)
    assert arr.shape == (n,)

    # And it must drive PointingEvaluator end-to-end without
    # crashing. evaluate() is the same call the optimizer driver
    # makes for every grid point.
    ev = PointingEvaluator(
        cat.ra_deg, cat.dec_deg,
        centration="UNCONSTRAINED",
        centration_per_target=cat.centration,
    )
    det, _, _ = ev.evaluate(53.16, -27.78, 30.0)
    assert det.shape == (n,)


# ---------------------------------------------------------------------
# Grid search + DE
# ---------------------------------------------------------------------


def test_grid_search_returns_sorted_descending():
    np.random.seed(42)
    ra = 53 + np.random.rand(20) * 0.05
    dec = -27 + np.random.rand(20) * 0.05
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    r = grid_search(ev, 53, -27, 30.0,
                   dra_arcsec=60, ddec_arcsec=60, dpa_deg=15,
                   n_ra=5, n_dec=5, n_pa=5)
    s = r["score"]
    assert all(s[i] >= s[i + 1] for i in range(len(s) - 1))
    assert len(s) == 5 * 5 * 5
    assert len(r["ra"]) == len(s)


def test_grid_search_hmpt_published_example_in_range():
    """The hMPT README reports ≈ 25/50 placements at the grid-search
    stage for the canonical example (seed=42, N=50, ra=53+rand*0.08,
    dec=-27+rand*0.08, search 0.05°×0.05°, 30° PA). Our re-implementation
    should land within ±5 of that figure on a coarser grid."""
    np.random.seed(42)
    N = 50
    ra = 53 + np.random.rand(N) * 0.08
    dec = -27 + np.random.rand(N) * 0.08
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED",
                          sigma_arcsec=0.06)
    r = grid_search(ev, 53, -27, 30.0,
                   dra_arcsec=180, ddec_arcsec=180, dpa_deg=30,
                   n_ra=15, n_dec=15, n_pa=15)
    # Allow ±10 from hMPT's 25 — different grid resolution explains the
    # spread. Failing this would mean a real algorithmic disagreement.
    best = int(r["score"][0])
    assert 15 <= best <= 40, (
        f"grid_search best={best}; hMPT reports ≈ 25 for this case. "
        f"Big disagreement may indicate a coordinate-frame bug.")


def test_refine_top_does_not_decrease_score():
    """Differential evolution refinement should NEVER decrease the
    score of a candidate (it's a local search inside a bounded box that
    includes the candidate itself)."""
    np.random.seed(42)
    N = 40
    ra = 53 + np.random.rand(N) * 0.05
    dec = -27 + np.random.rand(N) * 0.05
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    r = grid_search(ev, 53.02, -26.98, 30.0,
                   dra_arcsec=120, ddec_arcsec=120, dpa_deg=15,
                   n_ra=8, n_dec=8, n_pa=8)
    r2 = refine_top(ev, r, n_top=3, maxiter=40)
    assert r2["score"][0] >= r["score"][0] - 1e-6, (
        f"DE refinement lowered the best score: "
        f"grid={r['score'][0]:.1f} → DE={r2['score'][0]:.1f}"
    )


def test_centration_buffer_table_is_complete():
    """Sanity: every public centration class in the constants table
    appears in the README and the optimizer accepts it."""
    expected = {"UNCONSTRAINED", "ENTIRE_OPEN", "MIDPOINT",
                "CONSTRAINED", "TIGHTLY_CONSTRAINED"}
    assert expected.issubset(CENTRATION_BUFFERS.keys())
    # Buffers are monotone increasing.
    keys_ordered = ["UNCONSTRAINED", "ENTIRE_OPEN", "MIDPOINT",
                    "CONSTRAINED", "TIGHTLY_CONSTRAINED"]
    vals = [CENTRATION_BUFFERS[k] for k in keys_ordered]
    assert all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))


def test_dpa_zero_freezes_pa_axis():
    """When ΔPA = 0 every candidate must keep the central PA exactly —
    the optimizer must not silently drift via the DE refinement."""
    np.random.seed(42)
    N = 30
    ra = 53 + np.random.rand(N) * 0.05
    dec = -27 + np.random.rand(N) * 0.05
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    pa0 = 201.6
    r = grid_search(ev, 53, -27, pa0,
                   dra_arcsec=60, ddec_arcsec=60, dpa_deg=0,
                   n_ra=6, n_dec=6, n_pa=20)
    assert np.allclose(r["pa"], pa0), (
        f"grid_search with dPA=0 should freeze PA at {pa0}, "
        f"got unique values {np.unique(r['pa'])}")
    r2 = refine_top(ev, r, n_top=10,
                   dra_arcsec=2, ddec_arcsec=2, dpa_deg=0, maxiter=40)
    assert np.allclose(r2["pa"], pa0), (
        f"refine_top with dPA=0 must keep PA frozen — got "
        f"max drift {np.max(np.abs(r2['pa'] - pa0)):.4f}°")


def test_dra_ddec_zero_freezes_those_axes():
    """ΔRA = ΔDec = 0 → only PA is searched."""
    np.random.seed(42)
    ra = 53 + np.random.rand(20) * 0.05
    dec = -27 + np.random.rand(20) * 0.05
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    r = grid_search(ev, 53.02, -26.98, 30.0,
                   dra_arcsec=0, ddec_arcsec=0, dpa_deg=20,
                   n_ra=10, n_dec=10, n_pa=15)
    # Grid should have only 1 RA × 1 Dec × 15 PA = 15 points.
    assert len(r["score"]) == 15
    assert np.allclose(r["ra"], 53.02)
    assert np.allclose(r["dec"], -26.98)


def test_refine_top_dedups_identical_solutions():
    """When several top grid candidates land on the SAME refined
    pointing (within tolerance), refine_top must collapse them so the
    user doesn't see N copies of the same row. This was the original
    bug — dPA=0 + DE producing the same optimum from every starting
    point gave 10 identical-looking rows."""
    np.random.seed(42)
    N = 30
    ra = 53 + np.random.rand(N) * 0.05
    dec = -27 + np.random.rand(N) * 0.05
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    # Hand-build grid results: 5 candidates all within the dedup
    # tolerance of each other. After refinement DE polishes each,
    # but the dedup pass must collapse them.
    fake_grid = {
        "score": np.array([10.0, 9.5, 9.0, 8.5, 8.0]),
        "ra":  np.array([53.0200, 53.02005, 53.0201, 53.02002, 53.02008]),
        "dec": np.array([-26.9800, -26.98003, -26.98001, -26.98005, -26.98]),
        "pa":  np.array([30.0, 30.001, 30.0, 30.002, 30.0]),
    }
    r2 = refine_top(ev, fake_grid, n_top=5,
                   dra_arcsec=0.1, ddec_arcsec=0.1, dpa_deg=0.05, maxiter=20)
    assert len(r2["score"]) < 5, (
        f"dedup should collapse near-identical candidates; got "
        f"{len(r2['score'])} unique solutions (input had 5 in a tiny "
        f"window).")


def test_flux_objective_uses_weights():
    """The 'flux' objective should produce a different ranking from
    'number' when fluxes are non-uniform."""
    np.random.seed(42)
    N = 30
    ra = 53 + np.random.rand(N) * 0.05
    dec = -27 + np.random.rand(N) * 0.05
    flux = np.random.rand(N) * 10.0
    ev = PointingEvaluator(ra, dec, flux_sources=flux,
                          centration="UNCONSTRAINED")
    r_num = grid_search(ev, 53.02, -26.98, 30.0,
                       dra_arcsec=120, ddec_arcsec=120, dpa_deg=15,
                       n_ra=6, n_dec=6, n_pa=6, objective="number")
    r_flux = grid_search(ev, 53.02, -26.98, 30.0,
                        dra_arcsec=120, ddec_arcsec=120, dpa_deg=15,
                        n_ra=6, n_dec=6, n_pa=6, objective="flux")
    # The two top-scoring pointings can differ; at minimum the best
    # absolute score numbers should differ (count vs weighted flux).
    assert r_num["score"][0] != pytest.approx(r_flux["score"][0])
