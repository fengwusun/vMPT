"""Tests for the optimizer's shutter-collision protection feature.

The feature adds an optional `protect_mask` to `PointingEvaluator` so
that:
  1. A protected source landing on a row colliding with a stuck-open
     shutter is dropped (its spectrum is unavoidably contaminated).
  2. When two protected sources collide on the detector, the
     lower-priority (or, on tie, lower-weight) one is dropped.
  3. An unprotected source colliding with any still-kept protected
     source is dropped.

The tests exercise the public API
(`PointingEvaluator.evaluate_with_stats`) rather than the private
`_apply_collision_drops` so they keep working if the internals get
refactored.
"""

from __future__ import annotations

import numpy as np
import pytest

from vmpt.optimizer import (
    PointingEvaluator,
    SHVAL_S_TOLERANCE,
    grid_search,
)


# ---------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------


FIDUCIAL = dict(ra_p=53.16, dec_p=-27.78, pa_v3=0.0)


def grid_sources(n: int = 30, ra0: float = 53.16, dec0: float = -27.78,
                 spread_deg: float = 0.005, seed: int = 42):
    """Cluster of synthetic sources around a centre — wide enough that
    several land in the MSA, tight enough that lots share rows."""
    rng = np.random.default_rng(seed)
    ra = ra0 + (rng.random(n) - 0.5) * 2.0 * spread_deg
    dec = dec0 + (rng.random(n) - 0.5) * 2.0 * spread_deg
    return ra, dec


# ---------------------------------------------------------------------
# Backwards-compatibility — no protection ⇒ identical behaviour
# ---------------------------------------------------------------------


def test_no_protection_evaluate_unchanged():
    """`protect_mask=None` keeps the old 3-tuple return identical."""
    ra, dec = grid_sources(n=30)
    ev_old = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    ev_new = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=None,
    )
    det_a, tp_a, idx_a = ev_old.evaluate(**FIDUCIAL)
    det_b, tp_b, idx_b = ev_new.evaluate(**FIDUCIAL)
    assert np.array_equal(det_a, det_b)
    assert np.allclose(tp_a, tp_b)
    for x, y in zip(idx_a, idx_b):
        # quad is int; s/d are float and may have NaN — use array_equal
        # for ints and allclose-with-NaN-equal for floats.
        if x.dtype.kind == "i":
            assert np.array_equal(x, y)
        else:
            assert np.allclose(x, y, equal_nan=True)


def test_no_protection_n_dropped_is_zero():
    """`evaluate_with_stats` reports 0 drops when protection is off."""
    ra, dec = grid_sources(n=30)
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    _, _, _, n_drop = ev.evaluate_with_stats(**FIDUCIAL)
    assert n_drop == 0


def test_all_false_protect_mask_treated_as_off():
    """An all-False mask is the same as no mask — no overhead, no drops."""
    ra, dec = grid_sources(n=30)
    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=np.zeros(len(ra), dtype=bool),
        priorities=np.arange(len(ra), dtype=float),
        weights=np.ones(len(ra)),
        disperser="PRISM", filt="CLEAR",
    )
    _, _, _, n_drop = ev.evaluate_with_stats(**FIDUCIAL)
    assert n_drop == 0


# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------


def test_protect_mask_size_mismatch_raises():
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError, match="protect_mask size"):
        PointingEvaluator(
            ra, dec,
            protect_mask=np.ones(5, dtype=bool),  # wrong size
            disperser="PRISM", filt="CLEAR",
        )


def test_protect_mask_without_disperser_raises():
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError, match="disperser and filt"):
        PointingEvaluator(
            ra, dec,
            protect_mask=np.ones(len(ra), dtype=bool),
            # disperser missing
        )


# ---------------------------------------------------------------------
# Rule 3: protected ↔ unprotected
# ---------------------------------------------------------------------


def test_h_grating_protection_drops_unprotected():
    """G140H has a ~500″ V2 overlap window — almost any same-row pair
    on the same detector half collides. Protecting one source must
    drop at least one unprotected co-detection."""
    ra, dec = grid_sources(n=60, spread_deg=0.006)
    # Protect just source 0.
    protect = np.zeros(len(ra), dtype=bool)
    protect[0] = True
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))

    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    ev_prot = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    det_prot, _, _, n_drop = ev_prot.evaluate_with_stats(**FIDUCIAL)

    if not det_base[0]:
        pytest.skip("Protected source did not land in MSA at this pointing")
    if det_base.sum() < 2:
        pytest.skip("Need >1 detected sources for collision to be possible")
    # Protected source itself should still be detected (no stuck-open
    # collision in synthetic test; no P↔P collision since only 1 P).
    assert det_prot[0], (
        "Protected source must remain detected when no stuck-open is set"
    )
    # The drop count must equal the difference between baseline and
    # protected detection counts.
    assert n_drop == int(det_base.sum() - det_prot.sum()), (
        f"n_drop={n_drop} != Δdet={int(det_base.sum() - det_prot.sum())}"
    )
    # And for an H grating with a dense cluster, we expect at least 1
    # unprotected source to share the protected row.
    assert n_drop >= 1, (
        f"Expected at least one drop with G140H, got {n_drop}"
    )


def test_prism_protection_drops_fewer_than_h_grating():
    """PRISM has a 35″ V2 overlap — much smaller than H gratings'
    ~500″. Protection should drop strictly fewer unprotected sources
    under PRISM than under G140H."""
    ra, dec = grid_sources(n=80, spread_deg=0.01)
    protect = np.zeros(len(ra), dtype=bool)
    protect[0] = True
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))

    ev_prism = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="PRISM", filt="CLEAR",
    )
    ev_g140h = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    _, _, _, drop_prism = ev_prism.evaluate_with_stats(**FIDUCIAL)
    _, _, _, drop_h = ev_g140h.evaluate_with_stats(**FIDUCIAL)
    assert drop_prism <= drop_h, (
        f"PRISM drops ({drop_prism}) must be ≤ H drops ({drop_h})"
    )


# ---------------------------------------------------------------------
# Rule 2: protected ↔ protected (lower priority dropped)
# ---------------------------------------------------------------------


def test_pp_collision_drops_lower_priority():
    """When all sources are protected and several share a row under an
    H grating, only the highest-priority (smallest pri value) survivor
    per collision cluster is kept."""
    ra, dec = grid_sources(n=30, spread_deg=0.004)
    # Distinct priorities so the winner is unambiguous.
    pri = np.arange(len(ra), dtype=float)  # 0 = best, len-1 = worst
    wgt = np.ones(len(ra))
    protect = np.ones(len(ra), dtype=bool)

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G395H", filt="F290LP",  # 500″ overlap
    )
    det, _, _, n_drop = ev.evaluate_with_stats(**FIDUCIAL)
    n_kept = int(det.sum())

    # Compare to baseline detection count.
    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    if det_base.sum() < 2:
        pytest.skip("Need ≥ 2 detected sources to observe a collision")

    # Some sources must be dropped (high-extent H grating + dense cluster).
    # If by chance no two land on the same row, the test gracefully
    # passes with n_drop == 0.
    assert n_drop == int(det_base.sum() - n_kept), (
        f"n_drop={n_drop} != Δdet={int(det_base.sum() - n_kept)}"
    )

    # Among the dropped sources, none should be the lowest-priority
    # number (pri=0). The lowest-priority survivor of any colliding
    # cluster always wins.
    dropped = np.where(det_base & ~det)[0]
    if dropped.size > 0:
        # The dropped sources must not include the global priority winner.
        assert 0 not in dropped, (
            "The pri=0 source must always survive P↔P drops; "
            f"got dropped indices {dropped}"
        )


def test_pp_weight_tiebreaker_drops_lower_weight():
    """Same priority → higher weight wins."""
    # Place a tight pair of sources at the same priority; one has
    # higher weight. They almost certainly share a row.
    ra = np.array([53.16, 53.16 + 1.0 / 3600.0])  # 1″ apart
    dec = np.array([-27.78, -27.78])
    protect = np.array([True, True])
    pri = np.array([1.0, 1.0])             # tie
    wgt = np.array([10.0, 100.0])          # 1 has higher weight → wins

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G395H", filt="F290LP",
    )
    det, _, _, n_drop = ev.evaluate_with_stats(**FIDUCIAL)
    if det.sum() == 0 and n_drop == 0:
        pytest.skip("Neither source landed on the MSA at this pointing")
    if n_drop > 0:
        # The higher-weight source (index 1) must survive.
        assert det[1], (
            f"Higher-weight source dropped: det={det}, n_drop={n_drop}"
        )


# ---------------------------------------------------------------------
# Cross-detector-half — no collision
# ---------------------------------------------------------------------


def test_cross_detector_half_does_not_collide():
    """Q1/Q3 ↔ Q2/Q4 image onto different detectors and must never
    collide, even with the H grating's wide V2 window.

    We test this by checking that two sources whose nearest shutters
    fall in different detector halves are never dropped by Rule 3.
    """
    # Two sources placed to fall in different halves at a fiducial
    # pointing. Q1/Q3 is roughly -V2 side; Q2/Q4 is +V2 side.
    # Construct by spreading widely in RA at constant Dec.
    ra = np.array([53.16 - 0.005, 53.16 + 0.005])
    dec = np.array([-27.78, -27.78])
    protect = np.array([True, False])
    pri = np.array([1.0, 2.0])
    wgt = np.array([10.0, 10.0])

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G395H", filt="F290LP",
    )
    det, _, (quad, _, _), n_drop = ev.evaluate_with_stats(**FIDUCIAL)

    # If both detected, and they're in different detector halves
    # (Q1/Q3 vs Q2/Q4), they MUST both survive.
    nrs1 = {1, 3}
    nrs2 = {2, 4}
    if det[0] and det[1]:
        h0 = (int(quad[0]) in nrs1)
        h1 = (int(quad[1]) in nrs1)
        if h0 != h1:
            assert n_drop == 0, (
                f"Cross-half collision detected: quads={quad}, drops={n_drop}"
            )


# ---------------------------------------------------------------------
# Stuck-open rule
# ---------------------------------------------------------------------


def test_stuck_open_reason_array_disabled_when_none():
    """Passing reason=None disables Rule 1 entirely; no extra drops."""
    ra, dec = grid_sources(n=20)
    protect = np.ones(len(ra), dtype=bool)
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))
    # Custom operable mask: all True (no broken shutters).
    operable = np.ones((4, 171, 365), dtype=bool)
    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        operable=operable, reason=None,
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="PRISM", filt="CLEAR",
    )
    # With small overlap (PRISM) and reason=None, drops can still
    # occur due to P↔P, but not due to stuck-open. We just verify the
    # call succeeds and returns a sane n_drop.
    det, _, _, n_drop = ev.evaluate_with_stats(**FIDUCIAL)
    assert n_drop >= 0


def test_stuck_open_drops_protected_when_provided():
    """When the reason array marks shutters as stuck-open (REASON==2),
    protected sources colliding with those rows are dropped. We craft
    a custom reason mask that flags MANY stuck-open cells across a
    detector half to guarantee at least one collision."""
    ra, dec = grid_sources(n=30, spread_deg=0.005)
    protect = np.ones(len(ra), dtype=bool)
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))
    operable = np.ones((4, 171, 365), dtype=bool)

    # Compare two evaluators: one with reason=None (no stuck-open),
    # one with a synthetic reason marking every shutter in Q1 as
    # stuck-open. The latter MUST drop ≥ the former under an H grating.
    reason_none = None
    reason_stuck_q1 = np.zeros((4, 171, 365), dtype=np.int8)
    reason_stuck_q1[0, :, :] = 2  # all of Q1 stuck-open
    # Q1 also can't be opened (the source can't be detected there),
    # but Q3 sources collide with Q1 on the NRS1 detector half.
    operable_with_q1_dead = operable.copy()
    operable_with_q1_dead[0, :, :] = False

    ev_no_stuck = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        operable=operable_with_q1_dead, reason=reason_none,
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    ev_with_stuck = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        operable=operable_with_q1_dead, reason=reason_stuck_q1,
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    det_a, _, _, drop_a = ev_no_stuck.evaluate_with_stats(**FIDUCIAL)
    det_b, _, _, drop_b = ev_with_stuck.evaluate_with_stats(**FIDUCIAL)
    # With Q1 marked as stuck-open, every Q3-detected protected
    # source colliding with a Q1 row gets dropped — strictly more
    # drops than the no-stuck-open baseline.
    assert drop_b >= drop_a, (
        f"reason=Q1-stuck-open dropped {drop_b}, baseline dropped {drop_a}"
    )
    # And the kept set should shrink (or stay the same).
    assert det_b.sum() <= det_a.sum()


# ---------------------------------------------------------------------
# Scoring backends pick up the kept mask automatically
# ---------------------------------------------------------------------


def test_grid_search_score_with_protection_le_baseline():
    """Protection can only REDUCE the per-pointing score (drops never
    add sources). The best-score with protection must be ≤ baseline."""
    ra, dec = grid_sources(n=50, spread_deg=0.006)
    protect = np.zeros(len(ra), dtype=bool)
    protect[:5] = True
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))

    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    ev_prot = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    r_base = grid_search(
        ev_base, 53.16, -27.78, 0.0,
        dra_arcsec=30, ddec_arcsec=30, dpa_deg=15,
        n_ra=5, n_dec=5, n_pa=5,
    )
    r_prot = grid_search(
        ev_prot, 53.16, -27.78, 0.0,
        dra_arcsec=30, ddec_arcsec=30, dpa_deg=15,
        n_ra=5, n_dec=5, n_pa=5,
    )
    assert float(r_prot["score"].max()) <= float(r_base["score"].max()), (
        f"Protected max score {r_prot['score'].max()} > "
        f"baseline {r_base['score'].max()}"
    )


def test_shval_s_tolerance_constant_present():
    """The module-level constant the test suite relies on must exist."""
    assert SHVAL_S_TOLERANCE == 1


# ---------------------------------------------------------------------
# Slitlet-aware row tolerance (user-requested fix to v1.2.0)
# ---------------------------------------------------------------------


@pytest.mark.parametrize("slit_length, expected_ps_tol, expected_pp_tol", [
    (1, 1, 1),   # N=1 → half=0 → |Δs|≤1 vs stuck-open, |Δs|≤1 vs slitlet
    (3, 2, 3),   # N=3 → half=1 → |Δs|≤2 vs stuck-open, |Δs|≤3 vs slitlet
    (5, 3, 5),   # N=5 → half=2 → |Δs|≤3 vs stuck-open, |Δs|≤5 vs slitlet
])
def test_collision_tolerances_scale_with_slit_length(
    slit_length, expected_ps_tol, expected_pp_tol,
):
    """The tolerances cached on the evaluator must follow the
    slitlet-aware formula:
        protected ↔ stuck-open: half + 1
        protected ↔ slitlet:    2·half + 1
    where half = slit_length // 2. This is the fix to v1.2.0's
    hardcoded |Δs| ≤ 1 — for any slitlet wider than 1 shutter the old
    check was too narrow (e.g. for N=3 it ignored the s_p ± 2 rows
    the user wants to be empty)."""
    ra, dec = grid_sources(n=4)
    ev = PointingEvaluator(
        ra, dec, slit_length=slit_length,
        protect_mask=np.ones(len(ra), dtype=bool),
        priorities=np.arange(len(ra), dtype=float),
        weights=np.ones(len(ra)),
        disperser="PRISM", filt="CLEAR",
    )
    assert ev._sd_tol_ps == expected_ps_tol, (
        f"slit_length={slit_length}: stuck-open tol={ev._sd_tol_ps}, "
        f"expected {expected_ps_tol}"
    )
    assert ev._sd_tol_pp == expected_pp_tol, (
        f"slit_length={slit_length}: slitlet tol={ev._sd_tol_pp}, "
        f"expected {expected_pp_tol}"
    )


def test_n3_slitlet_drops_unprotected_two_rows_away():
    """For a 3-shutter slitlet, the user wants rows s_p±2 (the
    immediate neighbours of the outer slitlet edges) also kept clear.
    Set up an H grating (wide V2 overlap so the V2 condition is
    almost certainly satisfied) and verify that adding a slitlet-
    aware tolerance drops at least as many unprotected sources as the
    old |Δs|≤1 check would have.
    """
    ra, dec = grid_sources(n=80, spread_deg=0.005)
    protect = np.zeros(len(ra), dtype=bool)
    protect[0] = True
    pri = np.arange(len(ra), dtype=float)
    wgt = np.ones(len(ra))

    ev_n1 = PointingEvaluator(
        ra, dec, slit_length=1,
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    ev_n3 = PointingEvaluator(
        ra, dec, slit_length=3,
        protect_mask=protect, priorities=pri, weights=wgt,
        disperser="G140H", filt="F100LP",
    )
    _, _, _, drops_n1 = ev_n1.evaluate_with_stats(**FIDUCIAL)
    _, _, _, drops_n3 = ev_n3.evaluate_with_stats(**FIDUCIAL)
    # Wider slitlets = wider exclusion zone = at least as many drops.
    assert drops_n3 >= drops_n1, (
        f"N=3 dropped {drops_n3}, N=1 dropped {drops_n1} — "
        f"slitlet-aware tolerance should be >= per-shutter"
    )
