"""Tests for app.wavelengths: per-grating dispersion model."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.coords import MSA_V2_REF, MSA_V3_REF
from app.wavelengths import GRATING_RANGES, V2_DISP_EXTENT, cutoffs


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
    Path(__file__).resolve().parent.parent / "data" / "dispersion_cutoffs.npz"
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
