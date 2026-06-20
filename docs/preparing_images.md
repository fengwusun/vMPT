# Preparing your own image (RGB + WCS sidecar)

vMPT draws the MSA overlay **on an image of your field**. It accepts two
kinds of image:

1. **A FITS image with a WCS in its header** — load it directly (Input tab
   → *Load image…* → *FITS image*, or `--fits`). Nothing to prepare.
2. **A colour JPG/PNG + a small WCS "sidecar" FITS** — best for a nice
   multi-band RGB. This page shows how to make both from drizzled
   single-band mosaics (e.g. JWST/NIRCam `*_i2d.fits` / `*_drz.fits`).

```{important}
The RGB image and its WCS sidecar **must describe the same pixel grid** —
same width × height, same WCS. The simplest way to guarantee that is to
drizzle every band onto one common mosaic grid, then build the RGB and the
sidecar from those aligned mosaics.
```

## 1. Write the WCS sidecar

Take the WCS from any one of your aligned mosaics and write it to a tiny
FITS that carries only the WCS keywords plus the image size:

```python
from astropy.io import fits
from astropy.wcs import WCS

src = "F200W_drz.fits"          # any mosaic on your common grid
out = "wcs.fits"

with fits.open(src) as hdul:
    hdr = hdul[1].header        # JWST i2d/drz: WCS is in the SCI extension (1);
                                # use hdul[0].header if your WCS is in the primary HDU
    sidecar = WCS(hdr).to_fits()             # HDUList holding just the WCS
    for key in ("NAXIS", "NAXIS1", "NAXIS2"):
        sidecar[0].header[key] = hdr[key]    # record the grid size vMPT needs
    sidecar.writeto(out, overwrite=True)
```

`wcs.fits` is only a header — a few kB. vMPT reads `NAXIS1/NAXIS2` from it
to map the JPG's pixels onto the sky.

## 2. Build the RGB JPG

Pick three bands (bluest → reddest), stretch each, and stack them into an
RGB. Any three filters work — choose whatever shows your field best.

```python
import numpy as np
from astropy.io import fits
from PIL import Image, ImageEnhance

# Three single-band mosaics on the SAME pixel grid as the sidecar above,
# ordered bluest -> reddest.
blue_path  = "F090W_drz.fits"
green_path = "F200W_drz.fits"
red_path   = "F444W_drz.fits"

def stretch(path, floor=-1.5, ceil=0.5):
    """Log-stretch a surface-brightness image into [0, 1].
    Lower `floor` to bring up fainter features; raise `ceil` to tame
    bright cores. Tune per dataset."""
    data = fits.getdata(path).astype(float)         # SCI data of the mosaic
    x = np.log10(data + 10.0 ** floor)
    return np.nan_to_num((np.clip(x, floor, ceil) - floor) / (ceil - floor))

rgb = np.dstack([stretch(red_path),                 # R
                 stretch(green_path),               # G
                 stretch(blue_path)])               # B

# FITS arrays have their origin at the bottom-left; image files at the
# top-left. Flip vertically so the saved JPG lines up with the WCS sidecar
# (vMPT applies the complementary flip when it loads a JPG+sidecar).
rgb = np.flipud(rgb)

im = Image.fromarray(np.uint8(np.clip(rgb, 0.0, 1.0) * 255), "RGB")
im = ImageEnhance.Color(im).enhance(2.0)            # optional: punch up saturation
im.save("field_rgb.jpg", format="JPEG", subsampling=0, quality=100)
```

Tips:

- **Per-band balance** — if one band dominates, divide its stretched array
  by a factor before stacking (e.g. `stretch(green_path) / 1.8`) to even
  out the colours.
- **Faint vs bright** — the `floor`/`ceil` of `stretch()` set the displayed
  dynamic range in log10(flux). Start at `(-1.5, 0.5)` and adjust.
- **Format** — PNG works too; JPEG with `quality=100, subsampling=0` keeps
  the image sharp while staying small.

## 3. Load into vMPT

Input tab → **Load image…** → **JPG/PNG + WCS sidecar**, point at
`field_rgb.jpg` and `wcs.fits`. Or from the command line:

```bash
./run.sh --jpg field_rgb.jpg --wcs wcs.fits --catalog targets.csv
```

The MSA overlay should sit correctly on the field. If the image comes out
**vertically mirrored** relative to your catalog circles, remove (or add)
the `np.flipud` in step 2 — that's the only orientation knob, and it
depends on how your mosaics were written.
```{note}
The bundled `example_r0600` JPG + `wcs.fits` were made exactly this way, so
they're a working reference if you want to compare.
```
