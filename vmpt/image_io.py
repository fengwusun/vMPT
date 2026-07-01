"""Image loaders (FITS / JPG+sidecar) and display stretching."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from astropy.io import fits
from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    ManualInterval,
    PercentileInterval,
    SqrtStretch,
    ZScaleInterval,
)
from astropy.wcs import WCS
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# Default on-screen resolution cap for FITS. Larger images are strided-
# decimated on load so GB-scale FITS read fast and don't overflow the browser.
DEFAULT_FITS_MAX_DIM = 4000


@dataclass
class LoadedImage:
    data: np.ndarray
    wcs: WCS
    shape: tuple
    source_path: str
    mode: str
    wcs_sidecar_path: Optional[str] = None  # set only for jpg+sidecar mode
    factor: int = 1  # downsample factor applied on load (1 = full resolution)
    full_shape: Optional[tuple] = None  # (H, W) of the on-disk FITS (LOD source)
    hdu: int = 0  # HDU index the pixel data came from (for on-demand crops)


def _first_image_hdu(hdul: fits.HDUList) -> int:
    for i, h in enumerate(hdul):
        d = h.data
        if d is not None and getattr(d, "ndim", 0) == 2:
            return i
    raise ValueError("No 2D image HDU found")


def load_fits(
    path: str,
    hdu: Optional[int] = None,
    max_dim: int = DEFAULT_FITS_MAX_DIM,
) -> LoadedImage:
    """Load a 2D FITS image for display.

    Uses a memory-mapped read and, when the image is larger than ``max_dim``
    on a side, strided-decimates it on load (``data[::f, ::f]``) so GB-scale
    FITS read quickly (~1/f² of pages touched) and produce a display-sized
    array instead of a browser-overflowing full-res RGBA. The WCS is scaled to
    match, and the downsample ``factor`` is recorded for image-coordinate
    overlays.
    """
    with fits.open(path, memmap=True) as hdul:
        idx = _first_image_hdu(hdul) if hdu is None else hdu
        header = hdul[idx].header
        full = hdul[idx].data
        h0, w0 = full.shape[:2]
        factor = 1
        if max_dim and max(h0, w0) > max_dim:
            factor = int(np.ceil(max(h0, w0) / float(max_dim)))
            data = np.array(full[::factor, ::factor], dtype=np.float32)
        else:
            data = np.array(full, dtype=np.float32)
        # Keep NaN so the display can optionally render blank pixels white;
        # fold ±inf into NaN too so interval / percentile / histogram math
        # (all of which drop non-finite values) stays well-behaved.
        data[np.isinf(data)] = np.nan
        wcs0 = WCS(header)
        wcs = wcs0.celestial if wcs0.has_celestial else wcs0
        if factor > 1:
            wcs = _scale_wcs(wcs, factor, convention="stride")
    return LoadedImage(data=data, wcs=wcs, shape=data.shape,
                       source_path=path, mode="fits", factor=factor,
                       full_shape=(int(h0), int(w0)), hdu=int(idx))


def compute_interval(
    arr: np.ndarray,
    scale_mode: str = "percentile",
    percentile: float = 99.5,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> tuple:
    """Return the ``(lo, hi)`` intensity limits for a FITS array under the given
    scaling mode. Computed once from the base tier and reused for every
    on-demand zoom crop, so the image brightness stays stable as you zoom
    (a percentile taken per-crop would flicker)."""
    x = np.asarray(arr, dtype=np.float32)
    try:
        if scale_mode == "manual" and vmin is not None and vmax is not None:
            return float(vmin), float(vmax)
        if scale_mode == "zscale":
            lo, hi = ZScaleInterval().get_limits(x)
        else:
            lo, hi = PercentileInterval(float(percentile)).get_limits(x)
        lo, hi = float(lo), float(hi)
    except Exception:  # noqa: BLE001
        finite = np.isfinite(x)
        lo, hi = (np.percentile(x[finite], [1.0, 99.5])
                  if finite.any() else (0.0, 1.0))
        lo, hi = float(lo), float(hi)
    return lo, (hi if hi > lo else lo + 1.0)


def compute_histogram(arr: np.ndarray, nbins: int = 50,
                      clip: Optional[tuple] = None) -> dict:
    """Pixel-value histogram for the Image display dialog.

    Bins the finite values over the FULL data range ``[min, max]`` (so the
    binning is fixed once an image is loaded and the whole range is
    draggable), split into ``nbins`` equal-width bins. An RGB image bins
    luminance over ``[0, 255]``. Pass ``clip=(lo_pct, hi_pct)`` to bin over a
    robust percentile window instead. Values are sub-sampled above ~1M points
    to keep it fast. Returns a dict with ``edges`` (nbins+1), ``counts``
    (nbins), ``lo``/``hi`` (the binned range), and ``stats`` (min/median/max)."""
    a = np.asarray(arr)
    if a.ndim == 3 and a.shape[2] >= 3:
        rgb = a[..., :3].astype(np.float32)
        vals = (0.299 * rgb[..., 0] + 0.587 * rgb[..., 1]
                + 0.114 * rgb[..., 2]).ravel()
        lo, hi = 0.0, 255.0
    else:
        vals = np.asarray(a, dtype=np.float32).ravel()
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return {"edges": np.array([0.0, 1.0]), "counts": np.array([0]),
                    "lo": 0.0, "hi": 1.0, "stats": {}}
        if clip is not None:
            lo, hi = (float(x) for x in np.percentile(vals, list(clip)))
        else:
            lo, hi = float(vals.min()), float(vals.max())
        if hi <= lo:
            hi = lo + 1.0
    if vals.size > 1_000_000:                    # sub-sample huge arrays
        vals = vals[:: int(vals.size // 1_000_000)]
    counts, edges = np.histogram(vals, bins=int(nbins), range=(lo, hi))
    stats = {"min": float(vals.min()), "max": float(vals.max()),
             "median": float(np.median(vals))}
    return {"edges": edges, "counts": counts, "lo": float(lo),
            "hi": float(hi), "stats": stats}


def stretch_curve(stretch: str, vmin: float, vmax: float,
                  n: int = 48) -> tuple:
    """Trace how data values in ``[vmin, vmax]`` map to display brightness
    ``[0, 1]`` under a tone curve — for overlaying on the pixel histogram so
    the user sees the stretch they're looking at. Returns ``(xs, ys)``."""
    stretch_cls = _STRETCHES.get(str(stretch).lower(), AsinhStretch)
    xs = np.linspace(float(vmin), float(vmax), int(n))
    if vmax <= vmin:
        return xs.tolist(), [0.0] * len(xs)
    norm = (xs - float(vmin)) / (float(vmax) - float(vmin))
    try:
        ys = np.clip(np.nan_to_num(np.asarray(stretch_cls()(norm))), 0.0, 1.0)
    except Exception:  # noqa: BLE001
        ys = np.clip(norm, 0.0, 1.0)
    return xs.tolist(), np.asarray(ys).tolist()


def lod_view_factor(visible_full_px: float, target_px: float,
                    base_factor: int) -> int:
    """Choose a decimation factor for an on-demand FITS zoom crop.

    ``visible_full_px`` is the longer side of the visible region measured in
    ORIGINAL (full-res) pixels; ``target_px`` is the on-screen pixel budget
    (~canvas frame size). Returns a power-of-two factor in ``[1, base_factor]``
    so the rendered crop is ≳ ``target_px`` — capped at the base tier (never
    coarser) and floored at 1 (native resolution). Powers of two give a small,
    stable set of tiers (e.g. base_factor 4 → {1, 2, 4}) so we re-render only
    when crossing a tier."""
    if target_px <= 0 or visible_full_px <= 0:
        return int(max(1, base_factor))
    raw = float(visible_full_px) / float(target_px)
    f = 1
    while f * 2 <= raw:
        f *= 2
    return int(min(max(1, f), max(1, int(base_factor))))


def read_fits_region(
    path: str,
    hdu: int,
    y0: int,
    y1: int,
    x0: int,
    x1: int,
    view_factor: int,
) -> tuple:
    """Read a decimated crop ``data[y0:y1:f, x0:x1:f]`` (full-res pixel bounds)
    from the memory-mapped FITS. The read origin is aligned DOWN to a multiple
    of ``f`` so the sampled pixel grid is stable across pans (no shimmer).
    Returns ``(crop_float32, ay0, ax0)`` where ``(ay0, ax0)`` is the aligned
    top-left origin actually used."""
    f = int(max(1, view_factor))
    ay0 = int(y0) - (int(y0) % f)
    ax0 = int(x0) - (int(x0) % f)
    with fits.open(path, memmap=True) as hdul:
        full = hdul[int(hdu)].data
        crop = np.array(full[ay0:int(y1):f, ax0:int(x1):f], dtype=np.float32)
    crop[np.isinf(crop)] = np.nan   # keep NaN (see load_fits); neutralise ±inf
    return crop, ay0, ax0


def _scale_wcs(wcs: WCS, factor: int, convention: str = "resize") -> WCS:
    """Scale a WCS to a downsampled image.

    The CRPIX transform depends on HOW the image was downsampled — the two
    conventions differ by a half-pixel-ish shift that otherwise mis-places
    every overlay:

    - ``"stride"`` — strided decimation ``data[::f, ::f]`` (FITS path). The new
      pixel *j* IS old pixel *j·f*, so ``crpix_new = (crpix-1)/f + 1`` (FITS
      1-based). This keeps pixel *centers* aligned.
    - ``"resize"`` — area/bilinear ``Image.resize`` (JPG path), where the new
      pixel samples the block *center*, so ``crpix_new = (crpix-0.5)/f + 0.5``.

    Using the wrong one shifts overlays by ``0.5·(f-1)/f`` downsampled px
    (≈1.5 full-res px at f=4) — the cause of contours sitting off-source on
    large FITS before v1.8.0's fix.
    """
    w = wcs.deepcopy()
    c = np.asarray(wcs.wcs.crpix)
    if str(convention).lower() == "stride":
        w.wcs.crpix = (c - 1.0) / factor + 1.0
    else:
        w.wcs.crpix = (c - 0.5) / factor + 0.5
    if wcs.wcs.has_cd():
        w.wcs.cd = wcs.wcs.cd * factor
    else:
        w.wcs.cdelt = np.asarray(wcs.wcs.cdelt) * factor
    return w


def load_jpg_with_sidecar(
    jpg_path: str,
    sidecar_fits_path: str,
    max_dim: int = 8000,
) -> LoadedImage:
    im = Image.open(jpg_path)
    jpg_w, jpg_h = im.size

    with fits.open(sidecar_fits_path) as hdul:
        header = hdul[0].header.copy()

    naxis1 = header.get("NAXIS1")
    naxis2 = header.get("NAXIS2")
    crpix1 = header.get("CRPIX1")
    crpix2 = header.get("CRPIX2")

    if naxis1 is None or naxis2 is None:
        if crpix1 is not None and crpix2 is not None:
            implied_w = 2 * (crpix1 - 0.5)
            implied_h = 2 * (crpix2 - 0.5)
            if abs(implied_w - jpg_w) / max(implied_w, 1) > 0.1 or abs(implied_h - jpg_h) / max(implied_h, 1) > 0.1:
                print(
                    f"WARNING: JPG dims ({jpg_w}x{jpg_h}) disagree with CRPIX-implied "
                    f"sidecar dims ({implied_w:.0f}x{implied_h:.0f}); using JPG dims."
                )
        header["NAXIS"] = 2
        header["NAXIS1"] = jpg_w
        header["NAXIS2"] = jpg_h

    wcs = WCS(header)
    if wcs.has_celestial:
        wcs = wcs.celestial

    factor = 1
    if max(jpg_w, jpg_h) > max_dim:
        factor = int(np.ceil(max(jpg_w, jpg_h) / max_dim))
        new_w = jpg_w // factor
        new_h = jpg_h // factor
        im = im.resize((new_w, new_h), Image.BILINEAR)
        wcs = _scale_wcs(wcs, factor)

    arr = np.asarray(im)
    if arr.ndim == 2:
        shape = arr.shape
    else:
        shape = arr.shape[:2]

    return LoadedImage(
        data=arr,
        wcs=wcs,
        shape=shape,
        source_path=jpg_path,
        mode="jpg+sidecar",
        wcs_sidecar_path=sidecar_fits_path,
    )


_STRETCHES = {
    "linear": LinearStretch,
    "sqrt": SqrtStretch,
    "asinh": AsinhStretch,
    "log": LogStretch,
}


def _pack_rgba(rgb_u8: np.ndarray) -> np.ndarray:
    """(H, W, 3) uint8 → (H, W) uint32 RGBA (opaque), Bokeh image_rgba format."""
    h, w = rgb_u8.shape[:2]
    rgba = np.empty((h, w, 4), dtype=np.uint8)
    rgba[..., :3] = rgb_u8
    rgba[..., 3] = 255
    return rgba.view(np.uint32).reshape(h, w)


def _rgb_tone_asinh(x01: np.ndarray, a: float = 0.1) -> np.ndarray:
    """vMPT's legacy RGB tone-curve: an asinh lift that brings up faint
    structure and noise. This was applied to every RGB image before v1.8.0,
    so re-applying it keeps a loaded JPG/PNG looking exactly as it did
    (brightness/contrast then adjust on top). ``x01`` in [0, 1]."""
    return np.arcsinh(x01 / a) / np.arcsinh(1.0 / a)


# Colormaps offered for grayscale FITS (Image display dialog). All are
# matplotlib names; "gray" uses a fast pure-numpy path (no matplotlib import).
FITS_COLORMAPS = [
    "gray", "viridis", "magma", "inferno", "plasma", "cividis",
    "cubehelix", "hot", "afmhot", "turbo",
]


def _apply_colormap(g01: np.ndarray, cmap: str = "gray",
                    invert: bool = False) -> np.ndarray:
    """Map a normalized [0, 1] grayscale array to (H, W, 3) uint8 RGB through a
    matplotlib colormap. ``invert`` flips the mapping (dark↔bright). ``gray``
    without invert takes a fast pure-numpy path so the common case never
    imports matplotlib; an unknown name or a missing matplotlib falls back to
    grayscale."""
    g = 1.0 - g01 if invert else g01
    if str(cmap).lower() in ("gray", "grey", ""):
        u8 = np.clip(g * 255.0, 0, 255).astype(np.uint8)
        return np.repeat(u8[..., None], 3, axis=2)
    try:
        from matplotlib import colormaps
        rgb = colormaps[str(cmap)](g)[..., :3]        # (H, W, 3) float [0,1]
        return np.clip(rgb * 255.0, 0, 255).astype(np.uint8)
    except Exception:  # noqa: BLE001 — unknown cmap / no matplotlib → grayscale
        u8 = np.clip(g * 255.0, 0, 255).astype(np.uint8)
        return np.repeat(u8[..., None], 3, axis=2)


def stretch_for_display(
    arr: np.ndarray,
    stretch: str = "asinh",
    *,
    scale_mode: str = "percentile",
    percentile: float = 99.5,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    brightness: float = 0.0,
    contrast: float = 1.0,
    cmap: str = "gray",
    invert: bool = False,
    nan_white: bool = False,
) -> np.ndarray:
    """Convert raw pixels to a uint32 RGBA array for Bokeh's ``image_rgba``.

    RGB images (JPG/PNG, ndim==3): the legacy asinh tone-curve is re-applied
    (so the image looks exactly as it did pre-v1.8.0), then ``brightness`` and
    ``contrast`` adjust on top (fitsmap-style) — ``out = (tone(rgb)-0.5) *
    contrast + 0.5 + brightness``. Defaults (0, 1) reproduce the original look.

    Grayscale images (FITS, ndim==2): normalise with
    ``astropy.visualization.ImageNormalize`` using ``scale_mode``
    ('percentile' → central ``percentile`` %, 'manual' → ``[vmin, vmax]``,
    'zscale') and a ``stretch`` tone curve (linear/sqrt/asinh/log), then apply
    a ``cmap`` colormap (``invert`` flips dark↔bright). When ``nan_white`` is
    set, non-finite (NaN/blank/±inf) pixels are painted white instead of
    taking the colormap's zero colour.
    """
    # ---- RGB (JPG/PNG): legacy asinh tone-curve + brightness/contrast -----
    if arr.ndim == 3 and arr.shape[2] >= 3:
        rgb = arr[..., :3].astype(np.float32) / 255.0
        rgb = _rgb_tone_asinh(rgb)   # restore the pre-1.8.0 faint-structure look
        rgb = (rgb - 0.5) * float(contrast) + 0.5 + float(brightness)
        rgb = np.clip(rgb * 255.0, 0, 255).astype(np.uint8)
        return _pack_rgba(rgb)

    # ---- Grayscale / FITS: astropy ImageNormalize ------------------------
    x = np.asarray(arr, dtype=np.float32)
    bad = ~np.isfinite(x)          # NaN / blank / ±inf pixels
    if bad.any():
        # Treat ±inf exactly like NaN so a bad pixel never clips to the
        # interval max (white) on its own — it's a "blank", handled uniformly
        # by nan_white below.
        x = x.copy()
        x[bad] = np.nan
    stretch_cls = _STRETCHES.get(str(stretch).lower(), AsinhStretch)
    try:
        if scale_mode == "manual" and vmin is not None and vmax is not None:
            interval = ManualInterval(float(vmin), float(vmax))
        elif scale_mode == "zscale":
            interval = ZScaleInterval()
        else:  # "percentile"
            interval = PercentileInterval(float(percentile))
        norm = ImageNormalize(x, interval=interval, stretch=stretch_cls(),
                              clip=True)
        g01 = np.nan_to_num(np.clip(np.asarray(norm(x)), 0.0, 1.0), nan=0.0)
    except Exception:  # noqa: BLE001 — never let a bad range break display
        finite = ~bad
        lo, hi = (np.percentile(x[finite], [1.0, 99.5])
                  if finite.any() else (0.0, 1.0))
        if hi <= lo:
            hi = lo + 1.0
        g01 = np.nan_to_num(np.clip((x - lo) / (hi - lo), 0.0, 1.0), nan=0.0)
    rgb_u8 = _apply_colormap(g01, cmap=cmap, invert=invert)
    if nan_white and bad.any():
        rgb_u8 = np.array(rgb_u8, copy=True)
        rgb_u8[bad] = 255          # blank pixels → white
    return _pack_rgba(rgb_u8)
