"""Tests for vmpt.wavelengths: per-grating dispersion model."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from vmpt.coords import MSA_V2_REF, MSA_V3_REF
from vmpt.wavelengths import GRATING_RANGES, V2_DISP_EXTENT, cutoffs


CASES = [
    ("PRISM", "CLEAR", 0.60, 5.30),
    ("G140M", "F070LP", 0.70, 1.27),
    ("G140M", "F100LP", 0.97, 1.89),
    ("G235M", "F170LP", 1.66, 3.17),
    ("G395M", "F290LP", 2.87, 5.27),
    ("G140H", "F070LP", 0.70, 1.27),
    ("G140H", "F100LP", 0.97, 1.89),
    ("G235H", "F170LP", 1.66, 3.17),
    ("G395H", "F290LP", 2.87, 5.27),
]


@pytest.mark.parametrize("disp,filt,lam_blue,lam_red", CASES)
def test_fiducial_matches_published(disp, filt, lam_blue, lam_red):
    out = cutoffs(MSA_V2_REF, MSA_V3_REF, disp, filt)
    assert out["lam_blue"] is not None
    assert out["lam_red"] is not None
    assert abs(out["lam_blue"] - lam_blue) < 0.05
    assert abs(out["lam_red"] - lam_red) < 0.05
    # Gap should be between blue and red.
    assert out["lam_gap_lo"] is not None and out["lam_gap_hi"] is not None
    assert out["lam_blue"] < out["lam_gap_lo"] < out["lam_gap_hi"] < out["lam_red"]


def test_shift_with_v2_offset_G395M():
    """After clamping, a positive V2 shift moves the blue end UP while
    the red end stays clamped to lam_max (the spectrum's red end has
    run off the right edge of the detector). Symmetric for negative."""
    lam_min, lam_max = GRATING_RANGES["G395M"]["F290LP"]
    expected_shift = 30.0 * (lam_max - lam_min) / V2_DISP_EXTENT
    assert 0.2 < abs(expected_shift) < 0.7
    # Positive V2 shift → blue moves up, red clamped at lam_max.
    shifted_pos = cutoffs(MSA_V2_REF + 30.0, MSA_V3_REF, "G395M", "F290LP")
    assert shifted_pos["lam_red"] == pytest.approx(lam_max)
    assert shifted_pos["lam_blue"] == pytest.approx(lam_min + expected_shift)
    # Negative V2 shift → blue clamped (at filter cutoff = lam_min), red moves down.
    shifted_neg = cutoffs(MSA_V2_REF - 30.0, MSA_V3_REF, "G395M", "F290LP")
    assert shifted_neg["lam_blue"] == pytest.approx(lam_min)
    assert shifted_neg["lam_red"] == pytest.approx(lam_max - expected_shift)


def test_unsupported_combo_raises():
    with pytest.raises(ValueError):
        cutoffs(MSA_V2_REF, MSA_V3_REF, "G395M", "F100LP")


def test_prism_gap_matches_msaviz_fiducial():
    """The NRS1/NRS2 detector-gap wavelengths for PRISM/CLEAR are taken
    from spacetelescope/msaviz, which integrates the pipeline PRISM
    dispersion polynomial per shutter. At the central Q1 shutter
    msaviz reports gap edges at ≈ 1.87 / 3.93 μm. PRISM dispersion is
    too non-linear for the linear shift model to capture, so we hold
    these fixed across the MSA."""
    out = cutoffs(MSA_V2_REF, MSA_V3_REF, "PRISM", "CLEAR")
    assert out["lam_gap_lo"] == pytest.approx(1.87, abs=0.01)
    assert out["lam_gap_hi"] == pytest.approx(3.93, abs=0.01)
    # Gap width should be ≈ 2 μm — much wider than the previous
    # incorrect "10 % of span" (0.47 μm) approximation.
    assert (out["lam_gap_hi"] - out["lam_gap_lo"]) > 1.5


def test_prism_gap_does_not_shift_with_v2():
    """PRISM dispersion is non-linear; the gap location does not
    follow the linear V2 shift the gratings get. A ±30″ V2 offset
    must leave PRISM gap edges unchanged (within float noise)."""
    center = cutoffs(MSA_V2_REF, MSA_V3_REF, "PRISM", "CLEAR")
    plus = cutoffs(MSA_V2_REF + 30.0, MSA_V3_REF, "PRISM", "CLEAR")
    minus = cutoffs(MSA_V2_REF - 30.0, MSA_V3_REF, "PRISM", "CLEAR")
    for k in ("lam_gap_lo", "lam_gap_hi", "lam_blue", "lam_red"):
        assert plus[k] == pytest.approx(center[k])
        assert minus[k] == pytest.approx(center[k])


def test_grating_gap_still_shifts_linearly():
    """For the gratings we keep the linear shift model — their
    dispersion IS roughly linear in V2."""
    center = cutoffs(MSA_V2_REF, MSA_V3_REF, "G395M", "F290LP")
    plus = cutoffs(MSA_V2_REF + 30.0, MSA_V3_REF, "G395M", "F290LP")
    # The gap centre must move with V2 (within clamping).
    assert plus["lam_gap_lo"] is not None and center["lam_gap_lo"] is not None
    assert plus["lam_gap_lo"] != pytest.approx(center["lam_gap_lo"])


_TABLE_PATH = (
    Path(__file__).resolve().parent.parent / "vmpt" / "data"
    / "dispersion_cutoffs.npz"
)


def _has_table() -> bool:
    return _TABLE_PATH.exists()


def test_prism_per_shutter_lookup_uses_table_when_available():
    """A central PRISM Q1 shutter should land near msaviz's
    (0.6, 1.87, 3.93, 5.3); these are not the fiducial constants the
    fallback would return."""
    if not _has_table():
        pytest.skip("dispersion_cutoffs.npz not built — run scripts/precompute_dispersion_cutoffs.py")
    central = cutoffs(MSA_V2_REF, MSA_V3_REF, "PRISM", "CLEAR",
                      q=1, s=86, d=311)
    assert central["lam_gap_lo"] == pytest.approx(1.87, abs=0.08)
    assert central["lam_gap_hi"] == pytest.approx(3.93, abs=0.08)
    assert central["lam_blue"] == pytest.approx(0.60, abs=0.05)
    assert central["lam_red"] == pytest.approx(5.30, abs=0.05)


def test_prism_per_shutter_lookup_gap_varies_across_msa():
    """For PRISM the gap location SHOULD vary materially across the
    MSA — that's the whole reason a single fiducial value was wrong.
    Only Q1+Q2 shutters disperse across the gap; Q3/Q4 PRISM spectra
    fall on a single detector."""
    if not _has_table():
        pytest.skip("dispersion_cutoffs.npz not built")
    gap_los = []
    for (q, s, d) in [
        (1, 86, 311), (1, 13, 255), (1, 82, 288),
        (2, 1, 323), (2, 58, 336),
    ]:
        out = cutoffs(MSA_V2_REF, MSA_V3_REF, "PRISM", "CLEAR",
                     q=q, s=s, d=d)
        assert out["lam_gap_lo"] is not None, (q, s, d)
        gap_los.append(out["lam_gap_lo"])
    assert max(gap_los) - min(gap_los) > 0.5, (
        f"per-shutter PRISM gap_lo should vary by >0.5 μm, got {gap_los}"
    )


def test_prism_q3_q4_shutters_have_no_gap():
    """Q3/Q4 PRISM spectra fall entirely on a single detector."""
    if not _has_table():
        pytest.skip("dispersion_cutoffs.npz not built")
    for (q, s, d) in [(3, 96, 97), (4, 86, 200)]:
        out = cutoffs(MSA_V2_REF, MSA_V3_REF, "PRISM", "CLEAR",
                     q=q, s=s, d=d)
        assert out["lam_blue"] is not None, (q, s, d)
        assert out["lam_red"] is not None, (q, s, d)
        assert out["lam_gap_lo"] is None
        assert out["lam_gap_hi"] is None


# -- Grating tables ---------------------------------------------------

GRATING_COMBOS = [
    ("G140M", "F070LP"), ("G140M", "F100LP"),
    ("G235M", "F170LP"), ("G395M", "F290LP"),
    ("G140H", "F070LP"), ("G140H", "F100LP"),
    ("G235H", "F170LP"), ("G395H", "F290LP"),
]


@pytest.mark.parametrize("disp,filt", GRATING_COMBOS)
def test_grating_table_returns_in_range_endpoints(disp, filt):
    """Each grating combo's per-shutter lookup should return values
    inside the published sci_range when the shutter has any spectrum."""
    if not _has_table():
        pytest.skip("dispersion_cutoffs.npz not built")
    lam_min, lam_max = GRATING_RANGES[disp][filt]
    # Sweep a central shutter in each quadrant.
    found_any = False
    for q in (1, 2, 3, 4):
        out = cutoffs(MSA_V2_REF, MSA_V3_REF, disp, filt,
                     q=q, s=86, d=200)
        if out["lam_blue"] is None and out["lam_red"] is None:
            continue
        found_any = True
        if out["lam_blue"] is not None:
            assert lam_min - 0.02 <= out["lam_blue"] <= lam_max + 0.02, (
                f"{disp}/{filt} q={q} blue={out['lam_blue']}")
        if out["lam_red"] is not None:
            assert lam_min - 0.02 <= out["lam_red"] <= lam_max + 0.02
    assert found_any, f"{disp}/{filt} returned None for all quadrants tested"


# ── Tilt-slope map (cross-dispersion tilt of spectral traces) ──


def _has_tilt(disp: str, filt: str) -> bool:
    """The precompute may have run partially; tolerate that."""
    if not _has_table():
        return False
    import numpy as np
    keys = np.load(_TABLE_PATH).files
    return f"{disp}_{filt}_tilt_slope" in keys


def test_tilt_slope_for_missing_table_returns_zero():
    """When the table doesn't ship tilt arrays for a combo, the
    helper returns slope = 0 so the runtime overlap check falls back
    to the original flat-row behaviour."""
    from vmpt.wavelengths import tilt_slope_for_shutter
    # An obviously fake combo — the lookup must not raise.
    k = tilt_slope_for_shutter("FAKE", "FAKE", 1, 86, 183)
    assert k == 0.0


def test_tilt_slope_within_expected_magnitude():
    """For PRISM, the tilt slope should be small (sub-row at the
    spectrum edge). v2_overlap_distance(PRISM, CLEAR) is 18″, so a
    bound of ±0.05 rows/arcsec gives ±0.9 rows max tilt — well above
    what the precompute reports (|slope| < 0.02)."""
    if not _has_tilt("PRISM", "CLEAR"):
        pytest.skip("tilt arrays not built — run scripts/precompute_trace_tilt.py")
    from vmpt.wavelengths import tilt_slope_for_shutter
    for q in (1, 2, 3, 4):
        for (s, d) in [(86, 183), (40, 80), (140, 280), (170, 50)]:
            k = tilt_slope_for_shutter("PRISM", "CLEAR", q, s, d)
            assert abs(k) < 0.05, (
                f"PRISM/Q{q} s={s} d={d}: slope {k} too large; "
                "the precompute may need re-running"
            )


def test_tilt_slope_varies_across_field():
    """The tilt slope is field-dependent — its sign or magnitude must
    vary across well-separated grid corners of a quadrant."""
    if not _has_tilt("PRISM", "CLEAR"):
        pytest.skip("tilt arrays not built")
    from vmpt.wavelengths import tilt_slope_for_shutter
    slopes = []
    for q in (1, 2, 3, 4):
        for (s, d) in [(16, 34), (16, 331), (155, 34), (155, 331)]:
            slopes.append(tilt_slope_for_shutter("PRISM", "CLEAR", q, s, d))
    # Range across the 16 sampled positions must exceed 0.005 rows/arcsec
    # (PRISM has a smooth gradient of ~0.02 across the full field).
    spread = max(slopes) - min(slopes)
    assert spread > 0.005, (
        f"PRISM tilt slope range = {spread} rows/arcsec — expected "
        f"> 0.005; slopes={slopes}"
    )


def test_tilt_slope_bilinear_at_grid_corner_matches_table():
    """Sampling the slope helper at an exact grid corner should
    return the table value directly (no interpolation artefacts)."""
    if not _has_tilt("PRISM", "CLEAR"):
        pytest.skip("tilt arrays not built")
    import numpy as np
    from vmpt.wavelengths import tilt_slope_for_shutter, tilt_slope_map
    m = tilt_slope_map("PRISM", "CLEAR")
    assert m is not None
    grid_rows, grid_cols, slope = m
    q = 3  # Q3 is fully on NRS1 for PRISM
    ri, ci = 5, 5
    sval = int(grid_rows[ri])
    dval = int(grid_cols[ci])
    k = tilt_slope_for_shutter("PRISM", "CLEAR", q, sval, dval)
    expected = float(slope[q - 1, ri, ci])
    if not np.isfinite(expected):
        pytest.skip("table NaN at chosen corner")
    assert k == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize("disp,filt", GRATING_COMBOS)
def test_grating_gap_varies_across_msa(disp, filt):
    """For every grating combo, the per-shutter gap_lo should vary
    by at least ~50 nm across well-separated shutters — confirming
    we're reading per-shutter data, not a single fiducial."""
    if not _has_table():
        pytest.skip("dispersion_cutoffs.npz not built")
    gap_los = []
    for q in (1, 2, 3, 4):
        for (s, d) in [(86, 50), (86, 200), (86, 320), (40, 200), (140, 200)]:
            out = cutoffs(MSA_V2_REF, MSA_V3_REF, disp, filt,
                         q=q, s=s, d=d)
            if out["lam_gap_lo"] is not None:
                gap_los.append(out["lam_gap_lo"])
    # H gratings tend to have many gap-spanning shutters; M gratings
    # fewer. Require at least 3 gap-spanning samples and >50 nm spread.
    if len(gap_los) < 3:
        pytest.skip(f"{disp}/{filt}: too few gap-spanning shutters in sample")
    assert max(gap_los) - min(gap_los) > 0.05, (
        f"{disp}/{filt} gap_lo should vary > 0.05 μm, got {gap_los}"
    )


# ----------------------------------------------------------------------
# v2_overlap_distance — per-combo table & realism checks


def test_v2_overlap_distance_per_combo_table_loaded():
    """The (disperser, filter) lookup table should populate at import
    and return finite, positive *full* V2 extents for every supported
    combo. Catches typos / dict-key mismatches. Full extent is the
    detector x-span of the spectrum × ~0.077 ″/V2-px; the spec-overlap
    rule uses |ΔV2| < full_extent on a SAME-detector check."""
    from vmpt.wavelengths import v2_overlap_distance
    for d, f in [
        ("PRISM", "CLEAR"),
        ("G140M", "F070LP"), ("G140M", "F100LP"),
        ("G235M", "F170LP"), ("G395M", "F290LP"),
        ("G140H", "F070LP"), ("G140H", "F100LP"),
        ("G235H", "F170LP"), ("G395H", "F290LP"),
    ]:
        v = v2_overlap_distance(d, f)
        assert 20.0 <= v <= 400.0, f"{d}/{f}: full extent {v} out of plausible range"
    # PRISM is the most compact spectrum; H-gratings the longest.
    assert v2_overlap_distance("PRISM", "CLEAR") < v2_overlap_distance("G395M", "F290LP")
    assert v2_overlap_distance("G395M", "F290LP") < v2_overlap_distance("G395H", "F290LP")


def test_v2_overlap_distance_unknown_filter_falls_back_per_disperser():
    """An unknown filter under a known disperser should still return
    a usable value — the disperser-only fallback."""
    from vmpt.wavelengths import v2_overlap_distance
    assert v2_overlap_distance("G395M", "FAKE_FILTER") == v2_overlap_distance("G395M", "F290LP")
    # Unknown disperser → safe upper bound 300″
    assert v2_overlap_distance("FAKE", "FAKE") == 300.0


def test_v2_overlap_distance_g395m_full_extent():
    """G395M's full V2 extent on detector is ~103″. Same-detector pairs
    within this V2 separation collide; cross-detector pairs (e.g. Q4↔Q2
    G395M at same d) are filtered by the primary_detector lookup, not
    by this distance value."""
    from vmpt.wavelengths import v2_overlap_distance
    full = v2_overlap_distance("G395M", "F290LP")
    # Same-quadrant max ΔV2 ≈ 60″ — must always be inside the extent.
    assert 60.0 < full, f"G395M full extent {full}\" too small to catch same-q overlaps"
    # User case Q4 d=349 vs Q2 d=349 ΔV2 ≈ 87″ is *inside* the V2 distance,
    # but the cross-detector check (primary_detector) is what suppresses
    # the false flag for that pair, NOT this distance.
    delta_v2_q4_to_q2 = 86.9
    assert delta_v2_q4_to_q2 < full, (
        f"Q4↔Q2 ΔV2={delta_v2_q4_to_q2}\" should be inside the V2 distance "
        f"({full}\"); the suppression of Q2 flagging comes from the "
        "primary-detector check, not the distance check"
    )
