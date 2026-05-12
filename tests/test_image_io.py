"""Tests for app.image_io."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.image_io import load_fits, load_jpg_with_sidecar, stretch_for_display

A370 = "/Users/sunfengwu/nirspec/example_a370/a370_f182m_f200w_f210m.fits"
R0600_JPG = "/Users/sunfengwu/nirspec/example_r0600/JWST_F090W_F200W_F444W.jpg"
R0600_WCS = "/Users/sunfengwu/nirspec/example_r0600/wcs.fits"


def test_load_fits_a370():
    img = load_fits(A370)
    assert img.shape == (2500, 2200)
    assert img.data.dtype == np.float32
    assert not np.any(np.isnan(img.data))
    assert img.mode == "fits"
    crval = img.wcs.wcs.crval
    assert abs(crval[0] - 40.008) < 0.01
    assert abs(crval[1] - (-1.6)) < 0.01


def test_load_fits_explicit_hdu():
    img = load_fits(A370, hdu=0)
    assert img.shape == (2500, 2200)


def test_stretch_for_display_grayscale():
    img = load_fits(A370)
    rgba = stretch_for_display(img.data, stretch="asinh")
    assert rgba.shape == img.data.shape
    assert rgba.dtype == np.uint32


def test_stretch_for_display_rgb():
    arr = np.zeros((10, 12, 3), dtype=np.uint8)
    arr[..., 0] = 255
    rgba = stretch_for_display(arr, stretch="linear")
    assert rgba.shape == (10, 12)
    assert rgba.dtype == np.uint32


def test_stretch_modes():
    arr = np.linspace(0, 100, 64, dtype=np.float32).reshape(8, 8)
    for s in ("linear", "sqrt", "asinh", "log"):
        out = stretch_for_display(arr, stretch=s)
        assert out.shape == (8, 8)
        assert out.dtype == np.uint32


@pytest.mark.slow
def test_load_jpg_with_sidecar_downsampled():
    img = load_jpg_with_sidecar(R0600_JPG, R0600_WCS, max_dim=4000)
    # Either 4000x4000 (factor=4) or close
    assert max(img.shape) <= 4500
    assert img.data.ndim == 3
    assert img.data.shape[2] == 3
    assert img.mode == "jpg+sidecar"

    crval = img.wcs.wcs.crval
    assert abs(crval[0] - 90.047) < 0.01
    assert abs(crval[1] - (-20.133)) < 0.01

    # Pixel scale ~30 mas after scaling factor: 30 mas * factor at the orig resolution
    # Validate that downsampled CRPIX maps back to CRVAL within 0.1"
    crpix = img.wcs.wcs.crpix
    sky = img.wcs.pixel_to_world(crpix[0] - 1, crpix[1] - 1)
    # Separation from CRVAL should be tiny
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ref = SkyCoord(crval[0] * u.deg, crval[1] * u.deg)
    sep = sky.separation(ref).to(u.arcsec).value
    assert sep < 0.1, f"sep={sep:.4f} arcsec"
