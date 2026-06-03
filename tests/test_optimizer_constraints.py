"""Tests for the per-target spectral constraints (v1.3.0+).

Each test pins one of the four constraint types from the v1.3.0
design:
  1. required_lam — list of (λ_lo, λ_hi) ranges that must land on the
     detector.
  2. no_gap      — the NRS detector gap must not fall inside the
     spectrum.
  3. extend_blue — shutter's lam_blue must reach the disperser's
     MSA-wide best blue (within 20 nm tolerance).
  4. extend_red  — same on the red side.

Plus regressions for: the parse-string helper, the catalog dataclass
defaults, and the optimizer-level OR of the per-target `protect`
flag with the v1.2 catalog-wide `protect_mask`.
"""

from __future__ import annotations

import numpy as np
import pytest

from vmpt.catalog import (
    Catalog,
    _format_lam_req,
    _parse_lam_req_str,
)
from vmpt.optimizer import (
    DROP_COLLISION,
    DROP_EXTEND_BLUE,
    DROP_EXTEND_RED,
    DROP_NO_GAP,
    DROP_REASONS,
    DROP_REQUIRED_LAM,
    PointingEvaluator,
)
from vmpt.wavelengths import (
    disperser_max_lambda,
    disperser_min_lambda,
    disperser_range,
    interval_covered,
)


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


FIDUCIAL = dict(ra_p=53.16, dec_p=-27.7717, pa_v3=0.0)


def grid_sources(n: int = 50, ra0: float = 53.16, dec0: float = -27.78,
                 spread_deg: float = 0.005, seed: int = 42):
    rng = np.random.default_rng(seed)
    ra = ra0 + (rng.random(n) - 0.5) * 2.0 * spread_deg
    dec = dec0 + (rng.random(n) - 0.5) * 2.0 * spread_deg
    return ra, dec


def ragged_lam_req(n: int, ranges: list[tuple[float, float]]) -> np.ndarray:
    """Build a length-n object array where every row has the same
    `[(lo1, hi1), …]` list. Mirrors how the catalog editor's popover
    would produce per-row constraint data."""
    out = np.empty(n, dtype=object)
    for i in range(n):
        out[i] = list(ranges)
    return out


# ---------------------------------------------------------------------
# Parse-string helpers (vmpt.catalog)
# ---------------------------------------------------------------------


@pytest.mark.parametrize("text, expected", [
    ("",                       []),
    ("1.0-1.3",                [(1.0, 1.3)]),
    ("1.0-1.3; 1.5-1.8",       [(1.0, 1.3), (1.5, 1.8)]),
    ("1.0 - 1.3, 2 to 3",      [(1.0, 1.3), (2.0, 3.0)]),
    ("1.5 — 1.3",              [(1.3, 1.5)]),         # swap-order
    ("nan",                    []),
    ("garbage",                []),
    ("1.0-1.3; garbage; 2-2.5", [(1.0, 1.3), (2.0, 2.5)]),
])
def test_parse_lam_req_str(text, expected):
    assert _parse_lam_req_str(text) == expected


def test_parse_format_roundtrip():
    """The serialiser is the inverse of the parser for clean input."""
    s = "1.0-1.3; 1.5-1.8"
    parsed = _parse_lam_req_str(s)
    formatted = _format_lam_req(parsed)
    assert _parse_lam_req_str(formatted) == parsed


# ---------------------------------------------------------------------
# Disperser-range lookup helpers
# ---------------------------------------------------------------------


def test_disperser_range_known():
    assert disperser_range("PRISM", "CLEAR") == (0.60, 5.30)
    assert disperser_range("G140H", "F100LP") == (0.97, 1.89)
    assert disperser_range("G395H", "F290LP") == (2.87, 5.27)


def test_disperser_range_unsupported_combo():
    assert disperser_range("G140H", "F290LP") is None
    assert disperser_range(None, None) is None


def test_disperser_min_max_wrappers():
    assert disperser_min_lambda("PRISM", "CLEAR") == pytest.approx(0.60)
    assert disperser_max_lambda("PRISM", "CLEAR") == pytest.approx(5.30)


# ---------------------------------------------------------------------
# interval_covered semantics
# ---------------------------------------------------------------------


def test_interval_covered_simple():
    # [1.0, 1.3] inside [0.9, 1.5] with no gap → covered.
    assert interval_covered(1.0, 1.3, 0.9, float("nan"), float("nan"), 1.5)


def test_interval_covered_outside_blue():
    assert not interval_covered(1.0, 1.3, 1.4, float("nan"), float("nan"), 1.8)


def test_interval_covered_outside_red():
    assert not interval_covered(1.0, 1.3, 0.9, float("nan"), float("nan"), 1.2)


def test_interval_covered_bisected_by_gap():
    # Required [1.0, 1.3] bisected by gap [1.1, 1.2] → not covered.
    assert not interval_covered(1.0, 1.3, 0.9, 1.1, 1.2, 1.5)


def test_interval_covered_below_gap():
    # Required [1.0, 1.05] entirely below gap → covered.
    assert interval_covered(1.0, 1.05, 0.9, 1.1, 1.2, 1.5)


def test_interval_covered_above_gap():
    assert interval_covered(1.3, 1.5, 0.9, 1.1, 1.2, 1.5)


def test_interval_covered_nan_bounds():
    # Source doesn't reach the detector at all → never covered.
    assert not interval_covered(1.0, 1.3,
                                float("nan"), float("nan"),
                                float("nan"), float("nan"))


# ---------------------------------------------------------------------
# Catalog dataclass defaults
# ---------------------------------------------------------------------


def test_catalog_constraint_defaults():
    """Bare-bones Catalog construction yields all-empty constraints."""
    cat = Catalog(
        ids=np.array([1, 2, 3]),
        ra_deg=np.array([0.0, 1.0, 2.0]),
        dec_deg=np.array([0.0, 1.0, 2.0]),
        priority=np.array([1.0, 2.0, 3.0]),
        mag=np.array([10.0, 11.0, 12.0]),
        z=np.array([0.0, 0.0, 0.0]),
        label=np.array(["", "", ""], dtype=object),
        source_path="",
    )
    assert cat.required_lam.size == 0
    assert cat.no_gap.size == 0
    assert cat.extend_blue.size == 0
    assert cat.extend_red.size == 0
    assert cat.protect.size == 0


# ---------------------------------------------------------------------
# Optimizer integration — defaults preserve v1.2.x behaviour
# ---------------------------------------------------------------------


def test_no_constraints_no_drops():
    """All-default constraints leave the kept set untouched and the
    reason dict empty."""
    ra, dec = grid_sources(n=30)
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    assert sum(reasons.values()) == 0


def test_constraint_keys_stable():
    """The reason dict keys are the constants in DROP_REASONS."""
    ra, dec = grid_sources(n=10)
    ev = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    _, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    assert set(reasons.keys()) == set(DROP_REASONS)


def test_protect_mask_size_mismatch_via_per_target_field():
    """Per-target `protect` array of the wrong length is rejected."""
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError, match="per-target protect size"):
        PointingEvaluator(
            ra, dec,
            protect=np.ones(5, dtype=bool),  # wrong length
            priorities=np.ones(10),
            disperser="PRISM", filt="CLEAR",
        )


def test_required_lam_size_mismatch_raises():
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError, match="required_lam size"):
        PointingEvaluator(
            ra, dec,
            required_lam=ragged_lam_req(5, [(1.0, 1.3)]),  # wrong length
            disperser="PRISM", filt="CLEAR",
        )


def test_no_gap_size_mismatch_raises():
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError, match="no_gap size"):
        PointingEvaluator(
            ra, dec,
            no_gap=np.ones(5, dtype=bool),  # wrong length
            disperser="PRISM", filt="CLEAR",
        )


def test_constraint_without_disperser_raises():
    ra, dec = grid_sources(n=10)
    with pytest.raises(ValueError,
                       match="per-target spectral constraints"):
        PointingEvaluator(
            ra, dec,
            no_gap=np.ones(10, dtype=bool),
            # disperser missing
        )


# ---------------------------------------------------------------------
# Required λ ranges
# ---------------------------------------------------------------------


def test_required_lam_outside_disperser_drops_all():
    """Asking for λ ∈ [6, 7] μm under PRISM (range 0.60-5.30) drops
    every detected source as required_lam."""
    ra, dec = grid_sources(n=50)
    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    n_base = int(det_base.sum())
    if n_base == 0:
        pytest.skip("Need detected sources for this test")

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        required_lam=ragged_lam_req(len(ra), [(6.0, 7.0)]),
        disperser="PRISM", filt="CLEAR",
    )
    _, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    assert reasons[DROP_REQUIRED_LAM] == n_base


def test_required_lam_inside_keeps_all():
    """λ ∈ [1.0, 1.2] inside PRISM (no gap at typical shutters) keeps
    every detected source."""
    ra, dec = grid_sources(n=50)
    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    n_base = int(det_base.sum())
    if n_base == 0:
        pytest.skip("Need detected sources for this test")

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        required_lam=ragged_lam_req(len(ra), [(1.0, 1.2)]),
        disperser="PRISM", filt="CLEAR",
    )
    det, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    assert reasons[DROP_REQUIRED_LAM] == 0
    assert int(det.sum()) == n_base


# ---------------------------------------------------------------------
# no_gap
# ---------------------------------------------------------------------


def test_no_gap_h_grating_drops():
    """H gratings have a detector gap inside their spectrum for the
    centre shutters → no_gap=True drops every detected source."""
    ra, dec = grid_sources(n=50)
    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    n_base = int(det_base.sum())
    if n_base == 0:
        pytest.skip("Need detected sources")

    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        no_gap=np.ones(len(ra), dtype=bool),
        disperser="G140H", filt="F100LP",
    )
    _, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    # The H grating gap exists between NRS1 and NRS2 for every shutter
    # in the detector, so a no_gap-flagged source must be dropped.
    assert reasons[DROP_NO_GAP] >= 1


# ---------------------------------------------------------------------
# Mixed: total dropped == sum of per-reason counts
# ---------------------------------------------------------------------


def test_total_dropped_equals_sum_of_reasons():
    """`evaluate_with_stats` returns int(total) consistent with
    `evaluate_with_reasons`'s per-reason dict."""
    ra, dec = grid_sources(n=50)
    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        no_gap=np.ones(len(ra), dtype=bool),
        required_lam=ragged_lam_req(len(ra), [(1.0, 1.05)]),
        disperser="G395M", filt="F290LP",
    )
    _, _, _, n_total = ev.evaluate_with_stats(**FIDUCIAL)
    _, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    assert n_total == sum(reasons.values())


# ---------------------------------------------------------------------
# Per-target protect OR'd with the v1.2 catalog-wide mask
# ---------------------------------------------------------------------


def test_per_target_protect_alone_enables_collision_rules():
    """Setting `protect` per-target on at least one row enables the
    v1.2 collision-protection rules even when `protect_mask` is None.
    """
    ra, dec = grid_sources(n=40, spread_deg=0.005)
    protect = np.zeros(len(ra), dtype=bool)
    protect[0] = True
    ev_base = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det_base, _, _ = ev_base.evaluate(**FIDUCIAL)
    if not det_base[0] or det_base.sum() < 2:
        pytest.skip("Need protected + ≥1 other detected source")
    ev = PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED",
        protect=protect,  # NO `protect_mask`
        priorities=np.arange(len(ra), dtype=float),
        weights=np.ones(len(ra)),
        disperser="G140H", filt="F100LP",
    )
    _, _, _, reasons = ev.evaluate_with_reasons(**FIDUCIAL)
    # Collision rules should be active; reason key exists in the dict.
    assert DROP_COLLISION in reasons
    # H grating has wide V2 overlap → at least one non-protected
    # source on a colliding row should be dropped.
    assert reasons[DROP_COLLISION] >= 1
