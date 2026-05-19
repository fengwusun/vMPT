"""Image loaders (FITS / JPG+sidecar) and display stretching."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


@dataclass
class LoadedImage:
    data: np.ndarray
    wcs: WCS
    shape: tuple
    source_path: str
    mode: str
    wcs_sidecar_path: Optional[str] = None  # set only for jpg+sidecar mode


def _first_image_hdu(hdul: fits.HDUList) -> int:
    for i, h in enumerate(hdul):
        d = h.data
        if d is not None and getattr(d, "ndim", 0) == 2:
            return i
    raise ValueError("No 2D image HDU found")


def load_fits(path: str, hdu: Optional[int] = None) -> LoadedImage:
    with fits.open(path) as hdul:
        idx = _first_image_hdu(hdul) if hdu is None else hdu
        header = hdul[idx].header
        data = np.asarray(hdul[idx].data, dtype=np.float32)
        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
        wcs = WCS(header).celestial if WCS(header).has_celestial else WCS(header)
    return LoadedImage(data=data, wcs=wcs, shape=data.shape, source_path=path, mode="fits")


def _scale_wcs(wcs: WCS, factor: int) -> WCS:
    w = wcs.deepcopy()
    w.wcs.crpix = (np.asarray(wcs.wcs.crpix) - 0.5) / factor + 0.5
    if wcs.wcs.has_cd():
        w.wcs.cd = wcs.wcs.cd * factor
    elif wcs.wcs.has_pc():
        w.wcs.cdelt = np.asarray(wcs.wcs.cdelt) * factor
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


def _apply_stretch(x: np.ndarray, stretch: str) -> np.ndarray:
    if stretch == "linear":
        return x
    if stretch == "sqrt":
        return np.sqrt(np.clip(x, 0, 1))
    if stretch == "asinh":
        a = 0.1
        return np.arcsinh(x / a) / np.arcsinh(1.0 / a)
    if stretch == "log":
        a = 1000.0
        return np.log1p(a * np.clip(x, 0, 1)) / np.log1p(a)
    raise ValueError(f"Unknown stretch: {stretch}")


def stretch_for_display(
    arr: np.ndarray,
    stretch: str = "asinh",
    percentile_lo: float = 1.0,
    percentile_hi: float = 99.5,
) -> np.ndarray:
    if arr.ndim == 3 and arr.shape[2] >= 3:
        rgb = arr[..., :3].astype(np.float32) / 255.0
        stretched = _apply_stretch(rgb, stretch)
        stretched = np.clip(stretched * 255.0, 0, 255).astype(np.uint8)
        h, w = stretched.shape[:2]
        rgba = np.empty((h, w, 4), dtype=np.uint8)
        rgba[..., :3] = stretched
        rgba[..., 3] = 255
        return rgba.view(np.uint32).reshape(h, w)

    x = arr.astype(np.float32)
    finite = np.isfinite(x)
    if finite.any():
        lo, hi = np.percentile(x[finite], [percentile_lo, percentile_hi])
    else:
        lo, hi = 0.0, 1.0
    if hi <= lo:
        hi = lo + 1.0
    norm = np.clip((x - lo) / (hi - lo), 0.0, 1.0)
    s = _apply_stretch(norm, stretch)
    g = np.clip(s * 255.0, 0, 255).astype(np.uint8)
    h, w = g.shape
    rgba = np.empty((h, w, 4), dtype=np.uint8)
    rgba[..., 0] = g
    rgba[..., 1] = g
    rgba[..., 2] = g
    rgba[..., 3] = 255
    return rgba.view(np.uint32).reshape(h, w)
