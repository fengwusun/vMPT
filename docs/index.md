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

Interactive Bokeh app for planning JWST/NIRSpec MSA observations
**directly on an image of the target field**. vMPT combines:

- **Automated MSA pointing optimization** —
  `vmpt/optimizer.py` is a lightweight Python module that
  searches (RA, Dec, V3 PA) for the best MSA pointing and roll
  angle. Inspired by
  [hMPT](https://github.com/zihaowu-astro/hMPT) (Z. Wu et al.,
  CfA / Harvard), which is itself inspired by ESA's eMPT
  (Bonaventura et al. 2023).
  The MSA shutter geometry, V2/V3 ↔ (s, d) coordinate mapping,
  and constraint machinery are independently implemented; the
  search algorithm (grid → differential-evolution refine) is a
  simpler version than hMPT's. Three modes: Democracy (count),
  Meritocracy (Σ weight), Hierarchy (strict priority tiers).
- **Shutter-collision protection** — mark high-priority targets
  whose spectra must not overlap any other source on the
  detector under the current Disperser / Filter. Slitlet-aware
  row buffer (`v1.2.1+`) reserves one row above and below each
  protected slitlet, including against stuck-open shutters.
- **Hand-picking + live MPT-faithful conflict feedback** — click
  any shutter to open an N-shutter slitlet, watch the three-colour
  spec-overlap layer light up in real time. **Pink** (Mask Stuck)
  warns where a stuck-open's spectrum would land on an operable
  shutter; **orange** (Masked) warns where a user-pick's spectrum
  lands; **purple** (Mask Conflict) fires when two open slitlets
  touch (boundary rows adjacent, no operable row between) — bounded
  to the `±2`-row window where they crowd (not the whole band),
  matching APT MPT.
- **Image display + overlays** (`v1.8.0`) — a **🎨 Image display**
  dialog with FITS stretch / scaling (percentile · Min–Max · manual
  vmin–vmax · zscale) + a gradient-swatch **colormap** picker and a
  live **pixel histogram** (data range, [vmin, vmax], stretch curve);
  RGB brightness / contrast; and a **pixel-brightness histogram**.
  Large FITS open at a preview resolution and **auto-sharpen to native
  pixels as you zoom in**. Load DS9 **region** (`.reg`) and **contour**
  (`.ctr`/`.con`) files as sky-fixed **🧩 Add-on** layers — any mix at
  once, per-file on/off, recolour + fill, or `--addon` from the shell.
- **APT- / eMPT-ready export** — write a bundle (target-prefixed
  `<catalog>_MPT_plan.json` + `<catalog>_APT_catalog.cat`, with a
  generated `README.md` of the import steps) that loads straight into
  APT, plus the files that feed the
  [eMPT pipeline](https://github.com/esdc-esac-esa-int/eMPT_v1).
- **Sharing** — save the whole session as a JSON file, send it
  to a collaborator, and they pick up exactly where you left off.

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
