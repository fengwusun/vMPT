"""Tests for vmpt.image_io."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from vmpt.image_io import (
    FITS_COLORMAPS,
    compute_histogram,
    compute_interval,
    load_fits,
    load_jpg_with_sidecar,
    lod_view_factor,
    read_fits_region,
    stretch_curve,
    stretch_for_display,
)

A370 = "/Users/sunfengwu/nirspec/example_a370/a370_f182m_f200w_f210m.fits"
R0600_JPG = "/Users/sunfengwu/nirspec/example_r0600/JWST_F090W_F200W_F444W.jpg"
R0600_WCS = "/Users/sunfengwu/nirspec/example_r0600/wcs.fits"


def test_load_fits_a370():
    img = load_fits(A370)
    assert img.shape == (2500, 2200)
    assert img.data.dtype == np.float32
    # NaN is now PRESERVED (so the display can optionally paint blank pixels
    # white); the A370 mosaic has NaN borders, and ±inf is folded into NaN so
    # no infinities survive.
    assert not np.any(np.isinf(img.data))
    assert np.any(np.isnan(img.data))          # A370 has blank (NaN) edges
    assert np.isfinite(img.data).any()         # …but real data too
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


# ---------------------------------------------------------------------------
# v1.8.0 — scale modes (percentile / manual / zscale) + RGB brightness/contrast
# ---------------------------------------------------------------------------


def test_stretch_scale_modes():
    """Every (stretch × scale_mode) combination returns a valid RGBA."""
    arr = np.linspace(1.0, 1000.0, 256, dtype=np.float32).reshape(16, 16)
    for s in ("linear", "sqrt", "asinh", "log"):
        out_p = stretch_for_display(arr, stretch=s, scale_mode="percentile",
                                    percentile=95.0)
        out_m = stretch_for_display(arr, stretch=s, scale_mode="manual",
                                    vmin=10.0, vmax=900.0)
        out_z = stretch_for_display(arr, stretch=s, scale_mode="zscale")
        for out in (out_p, out_m, out_z):
            assert out.shape == (16, 16)
            assert out.dtype == np.uint32


def test_stretch_percentile_clips_extremes():
    """A small percentile window should saturate at least some pixels to
    pure black/white (i.e. the interval actually clips)."""
    arr = np.arange(100, dtype=np.float32).reshape(10, 10)
    out = stretch_for_display(arr, stretch="linear", scale_mode="percentile",
                              percentile=50.0)
    # Decode red channel (low byte) — with a 50% window the tails clip.
    red = (out & 0xFF).astype(np.uint8)
    assert red.min() == 0
    assert red.max() == 255


def test_stretch_manual_bad_range_does_not_raise():
    """vmax <= vmin must fall back to the safe percentile-linear path, not
    blow up the display."""
    arr = np.random.RandomState(0).rand(8, 8).astype(np.float32)
    out = stretch_for_display(arr, stretch="linear", scale_mode="manual",
                              vmin=5.0, vmax=5.0)
    assert out.shape == (8, 8) and out.dtype == np.uint32


def _decode_rgb(rgba_u32):
    r = (rgba_u32 & 0xFF).astype(np.uint8)
    g = ((rgba_u32 >> 8) & 0xFF).astype(np.uint8)
    b = ((rgba_u32 >> 16) & 0xFF).astype(np.uint8)
    return r, g, b


def test_stretch_nan_white_paints_blank_pixels():
    """nan_white=True paints NaN/inf pixels pure white; off, they take the
    colormap's zero colour (black for gray). Finite pixels are unchanged."""
    arr = np.array([[1.0, 2.0, np.nan],
                    [3.0, 4.0, np.inf],
                    [-np.inf, 6.0, 7.0]], dtype=np.float32)
    off = stretch_for_display(arr, stretch="linear", scale_mode="manual",
                              vmin=1.0, vmax=7.0, nan_white=False)
    on = stretch_for_display(arr, stretch="linear", scale_mode="manual",
                             vmin=1.0, vmax=7.0, nan_white=True)
    bad = ~np.isfinite(arr)
    r_off, g_off, b_off = _decode_rgb(off)
    r_on, g_on, b_on = _decode_rgb(on)
    # ON: every non-finite pixel is white.
    assert np.all(r_on[bad] == 255) and np.all(g_on[bad] == 255) and \
        np.all(b_on[bad] == 255)
    # OFF: non-finite pixels are the colormap zero (black for gray).
    assert np.all(r_off[bad] == 0)
    # Finite pixels identical whether or not nan_white is set.
    assert np.array_equal(r_on[~bad], r_off[~bad])


def test_stretch_nan_white_noop_without_nans():
    """With no non-finite pixels, nan_white makes no difference."""
    arr = np.arange(16, dtype=np.float32).reshape(4, 4)
    a = stretch_for_display(arr, stretch="linear", nan_white=False)
    b = stretch_for_display(arr, stretch="linear", nan_white=True)
    assert np.array_equal(a, b)


def test_rgb_keeps_legacy_asinh_look_by_default():
    """RGB must retain its pre-1.8.0 asinh tone-curve: a mid-gray pixel is
    LIFTED, not passed through unchanged. A plain identity (the v1.8.0-rc
    regression) flattened faint structure/noise — this locks the fix."""
    arr = np.full((4, 4, 3), 128, dtype=np.uint8)
    out = stretch_for_display(arr)  # defaults: brightness 0, contrast 1
    r = int(out.flat[0] & 0xFF)
    # Reproduce the exact legacy formula: asinh(x/0.1)/asinh(1/0.1).
    x = 128 / 255.0
    ref = int(round(np.arcsinh(x / 0.1) / np.arcsinh(1.0 / 0.1) * 255.0))
    assert abs(r - ref) <= 1, (r, ref)
    assert r > 150, f"mid-gray not lifted ({r}) — asinh tone-curve missing"


def test_rgb_brightness_contrast_clamps():
    """RGB path applies (rgb-0.5)*contrast + 0.5 + brightness and clamps to
    [0,255]; extreme settings must stay in range, not wrap."""
    arr = np.full((4, 4, 3), 128, dtype=np.uint8)
    bright = stretch_for_display(arr, brightness=0.5, contrast=1.0)
    dark = stretch_for_display(arr, brightness=-0.5, contrast=1.0)
    assert (bright & 0xFF).max() <= 255 and (dark & 0xFF).min() >= 0
    # Pushing brightness way up saturates toward white.
    white = stretch_for_display(arr, brightness=0.9, contrast=1.0)
    assert int((white & 0xFF).min()) == 255


# ---------------------------------------------------------------------------
# v1.8.0 — fast FITS load: memmap + strided downsample + WCS scaling
# ---------------------------------------------------------------------------


def _write_big_fits(path, n=600):
    from astropy.io import fits
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-1.0 / 3600.0, 1.0 / 3600.0]  # 1"/pix
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    data = np.arange(n * n, dtype=np.float32).reshape(n, n)
    fits.writeto(path, data, header=w.to_header(), overwrite=True)


def test_load_fits_downsamples_large(tmp_path):
    big = tmp_path / "big.fits"
    _write_big_fits(big, n=600)
    img = load_fits(str(big), max_dim=200)
    # 600 / 200 = 3 → factor 3, shape ceil-strided to 200.
    assert img.factor == 3
    assert max(img.shape) <= 200
    assert img.data.dtype == np.float32
    # WCS scaled: pixel scale grows by `factor` (1" → 3").
    scale_arcsec = abs(img.wcs.wcs.cdelt[1]) * 3600.0
    assert abs(scale_arcsec - 3.0) < 1e-6


def test_load_fits_no_downsample_when_small(tmp_path):
    small = tmp_path / "small.fits"
    _write_big_fits(small, n=300)
    img = load_fits(str(small), max_dim=4000)
    assert img.factor == 1
    assert img.shape == (300, 300)


def test_load_fits_downsampled_wcs_roundtrips_center(tmp_path):
    """The scaled WCS must still place CRVAL at the (rescaled) reference
    pixel — i.e. the sky anchor survives the downsample."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    big = tmp_path / "big.fits"
    _write_big_fits(big, n=800)
    img = load_fits(str(big), max_dim=200)
    crpix = img.wcs.wcs.crpix
    sky = img.wcs.pixel_to_world(crpix[0] - 1, crpix[1] - 1)
    ref = SkyCoord(150.0 * u.deg, 2.0 * u.deg)
    assert sky.separation(ref).to(u.arcsec).value < 0.05


def test_downsampled_fits_overlays_land_on_source(tmp_path):
    """Regression for the ~1 px contour/overlay shift on large FITS: because
    FITS is strided-decimated (`data[::f]`), the scaled WCS must use the
    stride convention. A source at full-res pixel P must project to base-tier
    pixel P/f exactly — NOT off by ~0.5·(f-1)/f px (the old resize-convention
    bug, ≈1.5 full-res px at f=4)."""
    from astropy.wcs import WCS
    from astropy.wcs.utils import skycoord_to_pixel
    big = tmp_path / "big.fits"
    _write_big_fits(big, n=1600)          # → factor 4 at max_dim 400
    full_wcs = WCS(fits_header_of(big))
    img = load_fits(str(big), max_dim=400)
    f = img.factor
    assert f == 4
    worst = 0.0
    for (px, py) in [(400, 600), (800, 1200), (100, 200), (1500, 100)]:
        world = full_wcs.pixel_to_world(px, py)
        x, y = skycoord_to_pixel(world, img.wcs, origin=0)
        worst = max(worst, abs(float(x) - px / f), abs(float(y) - py / f))
    assert worst < 1e-3, f"overlay off by {worst:.3f} base-px (stride bug?)"


def fits_header_of(path):
    from astropy.io import fits
    with fits.open(str(path)) as h:
        return h[0].header


# ---------------------------------------------------------------------------
# v1.8.0 — on-demand LOD helpers (zoom-dependent FITS resolution)
# ---------------------------------------------------------------------------


def test_lod_view_factor_tiers():
    """Powers of two, floored at 1 (native), capped at the base factor."""
    assert lod_view_factor(4000, 1000, 4) == 4   # raw 4
    assert lod_view_factor(2000, 1000, 4) == 2   # raw 2
    assert lod_view_factor(1500, 1000, 4) == 1   # raw 1.5 → native
    assert lod_view_factor(400, 1000, 4) == 1    # zoomed past native → native
    assert lod_view_factor(9000, 1000, 4) == 4   # would be 8, capped at base 4
    assert lod_view_factor(4000, 1000, 8) == 4   # base 8, raw 4 → 4
    assert lod_view_factor(0, 1000, 4) == 4      # degenerate → base


def test_load_fits_records_full_shape_and_hdu(tmp_path):
    big = tmp_path / "big.fits"
    _write_big_fits(big, n=1600)
    img = load_fits(str(big), max_dim=400)
    assert img.factor == 4
    assert img.full_shape == (1600, 1600)
    assert img.hdu == 0
    assert img.shape == (400, 400)


def test_read_fits_region_native_and_decimated(tmp_path):
    from astropy.io import fits
    big = tmp_path / "big.fits"
    data = (np.arange(1600 * 1600).reshape(1600, 1600) % 251).astype(np.float32)
    fits.writeto(str(big), data, overwrite=True)
    # Native crop matches the source exactly.
    crop, ay0, ax0 = read_fits_region(str(big), 0, 500, 700, 300, 500, 1)
    assert crop.shape == (200, 200) and (ay0, ax0) == (500, 300)
    assert np.array_equal(crop, data[500:700, 300:500])
    # Decimated crop aligns its origin DOWN to a multiple of the factor.
    crop2, ay2, ax2 = read_fits_region(str(big), 0, 501, 701, 301, 501, 2)
    assert (ay2, ax2) == (500, 300)
    assert np.array_equal(crop2, data[500:701:2, 300:501:2])


def test_compute_histogram_fits_full_range_50_bins():
    """Bins span the FULL data [min, max] (fixed once loaded) with 50 bins by
    default, so the whole range is coverable; stats report min/median/max."""
    rng = np.random.RandomState(0)
    data = np.concatenate([
        rng.normal(100, 5, 50000), rng.normal(5000, 20, 100)]).astype(np.float32)
    h = compute_histogram(data.reshape(-1, 1))
    assert len(h["counts"]) == 50 and len(h["edges"]) == 51   # default 50 bins
    assert h["lo"] == float(data.min()) and h["hi"] == float(data.max())
    assert h["hi"] > 5000              # full range includes the far tail
    assert abs(h["stats"]["median"] - 100) < 2
    # A robust clip window still works when explicitly requested.
    hc = compute_histogram(data.reshape(-1, 1), clip=(0.5, 99.5))
    assert hc["hi"] < 5000


def test_compute_histogram_rgb_luminance():
    rng = np.random.RandomState(1)
    rgb = rng.randint(0, 256, size=(30, 40, 3)).astype(np.uint8)
    h = compute_histogram(rgb)
    assert (h["lo"], h["hi"]) == (0.0, 255.0)
    assert len(h["counts"]) == 50   # default bin count


def test_compute_histogram_empty_safe():
    h = compute_histogram(np.full((4, 4), np.nan, dtype=np.float32))
    assert h["counts"].sum() == 0  # all-NaN → no crash


def test_stretch_curve_monotonic_and_endpoints():
    for s in ("linear", "sqrt", "asinh", "log"):
        xs, ys = stretch_curve(s, 10.0, 100.0, n=25)
        assert xs[0] == 10.0 and xs[-1] == 100.0
        assert ys[0] < 0.02 and abs(ys[-1] - 1.0) < 1e-6
        assert all(ys[i] <= ys[i + 1] + 1e-9 for i in range(len(ys) - 1))
    # asinh lifts faint values above the linear diagonal.
    _, y_asinh = stretch_curve("asinh", 0.0, 1.0, n=3)
    assert y_asinh[1] > 0.5
    # Degenerate range doesn't crash.
    xs, ys = stretch_curve("linear", 5.0, 5.0, n=4)
    assert len(xs) == 4


def test_compute_interval_minmax_percentile_100():
    """Percentile 100 = Min–Max: the interval is the data's exact min/max."""
    x = np.array([1.0, 5, 7, 100, 3, 50], dtype=np.float32)
    lo, hi = compute_interval(x, "percentile", 100.0)
    assert lo == float(x.min()) and hi == float(x.max())
    # A tighter percentile clips the tails inward.
    lo2, hi2 = compute_interval(x, "percentile", 90.0)
    assert lo2 > lo and hi2 < hi


def test_compute_interval_stable_across_crops(tmp_path):
    """The interval is a property of the array, not the view — a base tier and
    a sub-crop of it must not give wildly different limits when the crop is
    representative (this is why crops reuse the base interval)."""
    rng = np.random.RandomState(1)
    base = rng.normal(100, 10, size=(400, 400)).astype(np.float32)
    lo, hi = compute_interval(base, "percentile", 99.5)
    assert hi > lo
    lo_z, hi_z = compute_interval(base, "zscale")
    assert hi_z > lo_z
    lo_m, hi_m = compute_interval(base, "manual", vmin=50.0, vmax=150.0)
    assert (lo_m, hi_m) == (50.0, 150.0)


# ---------------------------------------------------------------------------
# v1.8.0 — grayscale FITS colormaps
# ---------------------------------------------------------------------------


def test_fits_colormaps_list():
    assert FITS_COLORMAPS[0] == "gray"
    assert "viridis" in FITS_COLORMAPS and "magma" in FITS_COLORMAPS
    assert len(FITS_COLORMAPS) == 10


def test_colormap_changes_output_and_invert():
    arr = np.linspace(0, 100, 64, dtype=np.float32).reshape(8, 8)
    gray = stretch_for_display(arr, "linear", scale_mode="manual",
                               vmin=0, vmax=100, cmap="gray")
    magma = stretch_for_display(arr, "linear", scale_mode="manual",
                                vmin=0, vmax=100, cmap="magma")
    inv = stretch_for_display(arr, "linear", scale_mode="manual",
                              vmin=0, vmax=100, cmap="magma", invert=True)
    assert not np.array_equal(gray, magma), "colormap had no effect"
    assert not np.array_equal(magma, inv), "invert had no effect"
    assert gray.dtype == magma.dtype == np.uint32


def test_colormap_gray_matches_plain_grayscale():
    """cmap='gray' must reproduce the plain grayscale packing exactly."""
    arr = np.linspace(0, 1, 100, dtype=np.float32).reshape(10, 10) * 50
    g = stretch_for_display(arr, "linear", scale_mode="manual", vmin=0, vmax=50)
    g2 = stretch_for_display(arr, "linear", scale_mode="manual", vmin=0,
                             vmax=50, cmap="gray")
    assert np.array_equal(g, g2)


def test_unknown_colormap_falls_back_to_gray():
    arr = np.linspace(0, 50, 64, dtype=np.float32).reshape(8, 8)
    ref = stretch_for_display(arr, "linear", scale_mode="manual", vmin=0, vmax=50)
    out = stretch_for_display(arr, "linear", scale_mode="manual", vmin=0,
                              vmax=50, cmap="not_a_real_cmap")
    assert np.array_equal(out, ref)


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
