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

## 4. Display controls — stretch, scale & large FITS

Once an image is loaded, click **Settings tab → 🎨 Image display…** to open the
image-display dialog and re-stretch it live — no reload needed. Changing any
control **keeps your current pan & zoom** (only the pixels are re-rendered). The
controls shown depend on the image type:

**FITS / grayscale images**

- **Stretch** — the tone curve applied to the data: **Linear**, **Square root**,
  **Asinh** (good default for high-dynamic-range fields), or **Log**.
- **Scale mode** — how the black/white limits are chosen:
  - **Percentile** — clip to the central *N*% of pixel values; presets
    **99.5 / 99 / 95 / 90** (lower = more aggressive, brings up faint features)
    plus **Min–Max** (100% — the full data range, no clipping).
  - **Manual** — type explicit `vmin` / `vmax` data values.
  - **ZScale** — the DS9/IRAF ZScale algorithm; robust for most astronomical
    images without tuning.

These map directly onto `astropy.visualization.ImageNormalize`, so the on-screen
result matches what you'd get in a notebook with the same interval + stretch.

- **Colormap** — apply a matplotlib colormap to the normalized image, chosen
  from a **rich dropdown** where every row is a **0→1 gradient bar** of that map
  (Gray, Viridis, Magma, Inferno, Plasma, Cividis, Cubehelix, Hot, Afmhot,
  Turbo), plus an **Invert colormap** checkbox to flip dark↔bright. The default
  is **gray** on every launch (the colormap is a per-session view tweak and is
  not persisted).
- **Render NaN as white** — a checkbox that paints **blank pixels** (NaN, and
  any ±inf) pure white instead of the colormap's darkest colour. Handy when a
  mosaic has NaN borders or gaps you'd rather see as white than as black. It's
  off by default and **persisted**, so once ticked it stays on across launches.

**RGB images (JPG/PNG)**

- **Brightness** and **Contrast** sliders (fitsmap-style). vMPT keeps the RGB
  image's original look — an asinh tone-curve that lifts faint structure — and
  applies brightness/contrast on top: `out = (tone(rgb) − 0.5) · contrast +
  0.5 + brightness`. At the defaults (0, 1) the image looks exactly as when
  loaded; the sliders just nudge the overall levels.

**Pixel histogram** — the dialog shows a histogram of the image's pixel values
over the **full data range** (50 bins, fixed once the image loads, with one bin
of padding each side so the edge handles are easy to grab), with a **log-scale
count** axis and a linear value axis. The blue band marks the **[vmin, vmax]**
currently shown, a red curve traces the **stretch** (how values map to
brightness), and a line underneath reports the data min / median / max and the
exact vmin / vmax. **Drag the handles** to set vmin / vmax right on the
histogram — the **△ (vmin)** sits at the bottom and the **▽ (vmax)** at the top,
at opposite ends. Dragging switches scaling to *manual*; for fine control, the
manual **vmin / vmax** boxes take exact numbers. For an RGB image it shows the
luminance distribution. Everything updates live as you change the stretch,
scaling, or vmin / vmax.

The same dialog also holds the **canvas size** (on-screen width/height of the
figure). All display settings are saved to `~/.vmpt/preferences.json` and
restored on the next launch. Slider re-stretches are throttled, so dragging stays
smooth even on large images. (Overlay-layer visibility and per-layer alpha /
stroke live in a separate **Settings → 🗂 Layers…** dialog.)

### Large (GB-scale) FITS — automatic zoom-in resolution

vMPT memory-maps FITS files and, when an image is larger than a display cap
(~4000 px on the long axis), shows a downsampled base view that opens in a
fraction of a second. **As you zoom in, it automatically re-reads just the
visible region at higher resolution — up to the file's exact native pixels at
full zoom** — reading on demand from the memory-mapped file (no pyramid to
precompute, no cache to manage). The refresh is debounced, so it sharpens a
moment after you stop panning/zooming. The WCS is rescaled so catalog targets
and overlays stay pixel-accurate at every zoom level. A 15k × 11k, 700 MB mosaic
opens in ~0.7 s and reaches native detail wherever you look.

## 5. DS9 region & contour overlays

Click **Input tab → 🧩 Load Add-on…** to open the add-on dialog, then
**Add files…** and pick any mix of DS9 **region** (`.reg`) and **contour**
(`.ctr` / `.con`) files — handy for showing apertures, masks, footprints, or
SNR/flux contours next to your targets.

- **One picker, both types, many files.** Select as many files as you like at
  once; vMPT decides per file whether it's a region or a contour (by extension,
  then by content), and for a contour it reads the coordinate frame (sky vs
  image) from the file itself. No pre-selection needed.
- **Per-file row in the sidebar.** Each loaded file gets a compact row under
  **Loaded add-ons** in the Input tab (right below *Loaded catalogs*) — a square
  **colour** swatch, an **on/off** checkbox with the filename, and **✕** to
  remove — so you toggle / recolour / delete a file without reopening the dialog.
  **Click the swatch** to open a small popover with a **colour picker** and a
  **fill-opacity** slider. **Clear all overlays** (in the dialog) removes them
  all. The on/off state, colour, and fill alpha all round-trip in the saved
  session.
- **From the command line.** `--addon` (repeatable) pre-loads overlay files:
  `./run.sh --fits img.fits --addon sources.reg --addon snr.ctr`.
- **Regions** — circle, ellipse, box, polygon, line, point, in any DS9 frame
  (`fk5`, `icrs`, `galactic`, `image`, …), parsed with the astropy-affiliated
  [`regions`](https://astropy-regions.readthedocs.io/) package.
- **Contours** — the vertex file DS9 writes from *Analysis → Contours → … →
  Apply* then *File → Save Contours*. Both the modern DS9 `.ctr` layout (a
  `level=N` marker per level with each contour wrapped in `( … )`) and a plain
  `.con` (whitespace `x y` per line, blank lines between segments) are accepted.

Overlays are projected through the image WCS, so they stay aligned through pan,
zoom, and re-pointing.

**Colour & fill are per file.** Because you can load several catalogs and DS9
files at once, each one's **colour** and **fill opacity** are set right where
it's loaded, not in a shared layer panel:

- **DS9 regions / contours** — each loaded file's row under **Loaded add-ons**
  (Input tab) has a square colour swatch; **click it** for a popover with a
  colour picker and a fill-opacity slider. Fill defaults to `0` (no face colour,
  like fitsmap); raise it toward `1` to shade the interior of closed regions /
  contours in that file's colour.
- **Catalog sources** — each catalog in the **catalog list** (Input tab) has its
  own square colour swatch that recolours just that catalog's markers.

The **Settings → 🗂 Layers** dialog keeps the layer-wide controls: show / hide
each layer (*Show DS9 regions*, *Show contours*, …) and tune its outline
**Alpha** (opacity) and **Stroke** (line width).

Loaded overlays (with per-file on/off state) are saved in the session JSON, so a
shared session reopens with the same overlays (as long as the files still exist
at their recorded paths).
