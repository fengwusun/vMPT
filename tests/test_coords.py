"""Tests for app.coords: shutter polygon construction and V2/V3 -> RA/Dec."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.coords import (
    MSA_V2_REF,
    MSA_V3_REF,
    rot_matrix,
    shutter_corners_v2v3,
    v2v3_to_radec,
)
from app.msa import load_msa_grid


def test_rot_matrix_identity():
    R = rot_matrix(0.0)
    assert np.allclose(R, np.eye(2))


def test_rot_matrix_90():
    R = rot_matrix(90.0)
    assert np.allclose(R, [[0, -1], [1, 0]])


def test_msa_refs_close_to_known():
    assert abs(MSA_V2_REF - 378.563) < 0.01
    assert abs(MSA_V3_REF - (-428.403)) < 0.01


def test_shutter_corners_shape_and_size():
    v2, v3 = load_msa_grid()
    v2c, v3c = float(v2[1, 85, 199]), float(v3[1, 85, 199])  # q=2, s=86, d=200
    corners = shutter_corners_v2v3(v2c, v3c)
    assert corners.shape == (4, 2)
    # Two adjacent edges should have lengths ~0.20 and ~0.46 arcsec.
    d01 = np.linalg.norm(corners[1] - corners[0])
    d12 = np.linalg.norm(corners[2] - corners[1])
    assert abs(d01 - 0.20) < 1e-9
    assert abs(d12 - 0.46) < 1e-9


def test_shutter_corners_sky_separation():
    v2, v3 = load_msa_grid()
    v2c, v3c = float(v2[1, 85, 199]), float(v3[1, 85, 199])
    corners = shutter_corners_v2v3(v2c, v3c)
    coord_c = SkyCoord(53.14, -27.79, unit=u.deg)
    sky = v2v3_to_radec(coord_c, 321.0, corners)
    assert len(sky) == 4
    # Edge separations on sky should match V2/V3 edge lengths to <1 mas.
    sep01 = sky[0].separation(sky[1]).to(u.arcsec).value
    sep12 = sky[1].separation(sky[2]).to(u.arcsec).value
    assert abs(sep01 - 0.20) < 1e-3
    assert abs(sep12 - 0.46) < 1e-3
