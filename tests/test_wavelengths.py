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
    ("G140M", "F100LP", 0.97, 1.84),
    ("G235M", "F170LP", 1.66, 3.07),
    ("G395M", "F290LP", 2.87, 5.14),
    ("G140H", "F070LP", 0.81, 1.27),
    ("G140H", "F100LP", 0.97, 1.84),
    ("G235H", "F170LP", 1.66, 3.07),
    ("G395H", "F290LP", 2.87, 5.14),
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
