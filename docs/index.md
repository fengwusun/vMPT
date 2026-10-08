# vMPT — visual MSA Planning Tool

```{toctree}
:hidden:
:maxdepth: 2

installation
quickstart
preparing_images
optimizer
multiconfig
catalogs
spec_overlap
exporting
troubleshooting
api
changelog
```

## vMPT in 20 seconds

```{raw} html
<video controls autoplay muted loop playsinline preload="metadata"
       poster="_static/vmpt_intro_poster.jpg"
       style="width:100%;max-width:960px;display:block;margin:0 auto 8px;border-radius:12px;background:#0b0e21">
  <source src="_static/vmpt_intro.mp4" type="video/mp4">
  Your browser does not play MP4 — <a href="_static/vmpt_intro.mp4">download the video</a>.
</video>
<p style="text-align:center;margin:0 0 1.4em;font-size:0.95em;opacity:0.85">
  Point the micro-shutter mask at a field and let the optimizer find the pointing that catches the most
  targets; pick shutters by hand on a zoomed grid; spectra spread sideways, masked shutters turn orange and
  a true conflict pulses violet. <a href="_static/vmpt_intro.html">Open the interactive version</a>
  (pause, replay, scrub).
</p>
```

vMPT is an interactive Bokeh app for planning JWST/NIRSpec micro-shutter
(MSA) observations **directly on an image of your field**: point the mask,
let the optimizer find the best pointing, pick shutters by hand, see
spectral overlaps as they happen, and export a bundle that loads straight
into APT.

**What it does**

- **Finds the best pointing** — searches RA, Dec and roll angle for the
  position that puts the most targets into operable shutters, scored by
  count, by weight, or by strict priority tiers. Multi-configuration plans
  (up to five pointings) come out of a single run. The search is inspired
  by [hMPT](https://github.com/zihaowu-astro/hMPT) (Z. Wu et al.,
  CfA / Harvard) and ESA's eMPT (Bonaventura et al. 2023); the MSA
  geometry, coordinate mapping and constraint machinery are vMPT's own.
- **Hand-picking with live conflict feedback** — click a shutter to open
  an N-shutter slitlet and see where every open and stuck-open shutter's
  spectrum lands: **orange** = masked by a spectrum, **pink** = masked by
  a stuck-open shutter, **purple** = two slitlets in conflict, matching
  APT's MPT. Undo at will; Space toggles a single shutter; W A S D pans.
- **Collision protection** — protect high-priority targets so that no
  other spectrum, including from stuck-open shutters, can overlap theirs.
- **Works on your image** — FITS with a WCS, or JPG/PNG with a WCS
  sidecar. GB-scale mosaics open instantly and sharpen as you zoom; live
  stretch, colormap and histogram controls; DS9 regions and contours as
  overlays, each with its own colour.
- **Round-trips with APT and eMPT** — export an APT-importable catalog and
  MPT plan (with a README of the import steps) plus the eMPT pipeline
  inputs; import APT plans, shutter masks or `.aptx` archives; save and
  share whole sessions as JSON.

---

## Where to go next

:::{list-table}
:header-rows: 0
:widths: 30 70

* - **[Installation](installation.md)**
  - `pip install jwst-vmpt`, or from source if you're a developer.
* - **[Quick start](quickstart.md)**
  - Load an example field, aim the MSA, pick shutters, export.
* - **[MSA pointing optimizer](optimizer.md)**
  - Democracy / Meritocracy / Hierarchy + shutter-collision
    protection.
* - **[Multiple configurations](multiconfig.md)**
  - Plan up to five APT-style MPT configs, the auto-all optimizer, the
    max-configs cap, the MPT catalog viewer, and multi-source shutters.
* - **[Catalogs](catalogs.md)**
  - CSV / ASCII / FITS columns, multi-catalog stacking, weights,
    priorities.
* - **[Spec-overlap colours](spec_overlap.md)**
  - Pink / orange / purple: what each colour means, how the
    detector-pixel collision check works, per-disperser behaviour.
* - **[Exporting & sharing](exporting.md)**
  - APT MPT plan, eMPT bundle, session JSON.
* - **[Troubleshooting](troubleshooting.md)**
  - Common install + runtime issues with their one-line fixes.
* - **[API reference](api.md)**
  - `vmpt.optimizer`, `vmpt.msa`, `vmpt.wavelengths`, …
* - **[Changelog](changelog.md)**
  - Release notes for every version.
:::

## Citation

If you use vMPT in a paper, please cite the underlying algorithms:

- **hMPT** — Wu, Z. et al. (in prep); CfA/Harvard. Lightweight
  Python script for optimizing MSA pointing and roll angles,
  inspired by ESA's eMPT (and the inspiration for vMPT's own
  optimizer).
  <https://github.com/zihaowu-astro/hMPT>
- **eMPT** — Bonaventura, N. et al. (2023), *A&A* 672, A40.
  ESA's reference MSA Planning Tool pipeline.

vMPT itself is a free reimplementation that adds the visual
hand-pick layer + shutter-collision protection on top.

## Links

- 📦 PyPI: <https://pypi.org/project/jwst-vmpt/>
- 🐙 GitHub: <https://github.com/fengwusun/vMPT>
- 🐞 Issues: <https://github.com/fengwusun/vMPT/issues>
- 📝 Changelog: [in this site](changelog.md) or
  [on GitHub](https://github.com/fengwusun/vMPT/blob/main/CHANGELOG.md)
