# vMPT — visual MSA Planning Tool

```{toctree}
:hidden:
:maxdepth: 2

installation
quickstart
optimizer
catalogs
exporting
api
changelog
```

Interactive Bokeh app for planning JWST/NIRSpec MSA observations
**directly on an image of the target field**. vMPT combines:

- **Automated MSA pointing optimization** — derived from
  [hMPT](https://github.com/zihaowu-astro/hMPT) (Z. Wu et al.,
  CfA / Harvard), itself a Python re-implementation of ESA's
  eMPT (Bonaventura et al. 2023). Grid + differential-evolution
  search over (RA, Dec, V3 PA) in three modes: Democracy
  (count), Meritocracy (Σ weight), Hierarchy (strict priority
  tiers).
- **Shutter-collision protection** — mark high-priority targets
  whose spectra must not overlap any other source on the
  detector under the current Disperser / Filter. Slitlet-aware
  row buffer (`v1.2.1+`) reserves one row above and below each
  protected slitlet, including against stuck-open shutters.
- **Hand-picking + live conflict feedback** — click any shutter
  to open an N-shutter slitlet, watch the orange spec-overlap
  layer light up in real time, undo / redo at will.
- **APT- / eMPT-ready export** — write an MPT_plan.json + .cat
  bundle that loads straight into APT, plus the three CSVs that
  feed the [eMPT pipeline](https://github.com/esdc-esac-esa-int/eMPT_v1).
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
* - **[Catalogs](catalogs.md)**
  - CSV / ASCII / FITS columns, multi-catalog stacking, weights,
    priorities.
* - **[Exporting & sharing](exporting.md)**
  - APT MPT plan, eMPT bundle, session JSON.
* - **[API reference](api.md)**
  - `vmpt.optimizer`, `vmpt.msa`, `vmpt.wavelengths`, …
* - **[Changelog](changelog.md)**
  - Release notes for every version.
:::

## Citation

If you use vMPT in a paper, please cite the underlying algorithms:

- **hMPT** — Wu, Z. et al. (in prep); CfA/Harvard. Python re-
  implementation of ESA's eMPT. <https://github.com/zihaowu-astro/hMPT>
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
