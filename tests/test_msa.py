"""Tests for app.msa: shutter grid loader and operability parser."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.msa import (
    load_msa_grid,
    load_operability,
    shutter_center_v2v3,
    shutters_in_bbox,
)


def test_grid_shape():
    v2, v3 = load_msa_grid()
    assert v2.shape == (4, 171, 365)
    assert v3.shape == (4, 171, 365)


def test_grid_v2_range_plausible():
    v2, v3 = load_msa_grid()
    # V2 ~ 226–534 across the four quadrants.
    assert 200 < v2.min() < 350
    assert 500 < v2.max() < 600
    # V3 ~ -584 to -276
    assert -600 < v3.min() < -500
    assert -350 < v3.max() < -200


def test_shutter_center_consistency():
    v2, v3 = load_msa_grid()
    cv2, cv3 = shutter_center_v2v3(2, 86, 200)
    assert cv2 == float(v2[1, 85, 199])
    assert cv3 == float(v3[1, 85, 199])


def test_shutters_in_bbox_picks_known_quadrant():
    v2, v3 = load_msa_grid()
    # Sub-region well inside Q2's bbox (V2 [316, 448], V3 [-406, -276]).
    out = shutters_in_bbox(350.0, 380.0, -350.0, -320.0)
    assert out.shape[1] == 3
    assert out.shape[0] > 0
    # All returned shutters should be from quadrant 2 (sub-bbox is Q2-only).
    qs = out[:, 0]
    assert (qs == 2).all(), f"expected only Q2, got {set(qs.tolist())}"


def test_shutters_in_bbox_round_trip():
    v2, v3 = load_msa_grid()
    cv2, cv3 = float(v2[2, 50, 100]), float(v3[2, 50, 100])  # q=3, s=51, d=101
    out = shutters_in_bbox(cv2 - 0.01, cv2 + 0.01, cv3 - 0.01, cv3 + 0.01)
    # The chosen shutter must appear.
    assert any((row[0] == 3 and row[1] == 51 and row[2] == 101) for row in out)


def test_operability_loader_fallback_when_no_file():
    # Point at a path that has no msaoper json.
    operable, reason = load_operability(crds_path="/nonexistent/path/that/has/no/crds")
    # The function also probes CRDS_PATH env and ~/crds_cache; if a real file
    # is present we just verify shape and dtype.
    assert operable.shape == (4, 171, 365)
    assert reason.shape == (4, 171, 365)
    assert operable.dtype == bool
    assert reason.dtype == np.int8
