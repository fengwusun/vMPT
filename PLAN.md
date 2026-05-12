# NIRSpec MSA Overlay & Planner — Implementation Plan

See `CONTEXT.md` for data formats, coordinate transforms, and the MPT-style JSON schema. This file is the milestone-level plan.

The end-state is an **MPT/eMPT-style planner**: load image → drop targets → pick pointing & APA_V3 → hand-pick shutters → export a config that round-trips with APT.

---

## Stack decision (settled 2026-05-11)

**Bokeh server from M1 onward.** Run with `bokeh serve app/` for local use.

- Hand-picking shutters is fundamentally "click on a polygon glyph" + hover tooltips on hundreds of thousands of polygons. Bokeh's `MultiPolygons` + `TapTool` + `HoverTool` is the natural fit.
- Avoids a mid-project rewrite (the earlier "Streamlit for M1–M3, switch later" plan was wasted work).
- Math/data layers (`coords.py`, `catalog.py`, `operability.py`, `wavelengths.py`, `empt_io.py`) are UI-agnostic Python — Bokeh only sits on top.

If file-upload UX in Bokeh gets painful, fallback is Streamlit + `streamlit-bokeh-events` keeping the Bokeh canvas — same canvas, different shell.

---

## Milestones

### M0 — Scaffolding (0.5 day)

- `app/`, `data/`, `notebooks/`.
- Copy `nirspec_msa_v2v3.npz` into `data/`.
- `app/coords.py` — pure functions ported from `footprint_emerald.ipynb`:
  - `rot_matrix(deg)`
  - `shutter_corners_v2v3(v2c, v3c, w=0.20, h=0.46)`
  - `v2v3_to_radec(coord_c, pa_v3, corners_v2v3) -> SkyCoord`
  - `load_msa_grid()` → `(v2_msa, v3_msa)` plus operable placeholder.
- `notebooks/00_sanity.ipynb`: reproduce one figure from `footprint_emerald.ipynb` on a JADES cutout to confirm orientation and sign conventions before anything else lands.

### M1 — Image + overlay (Bokeh server) (1.5 days)

Sidebar inputs (Bokeh widgets):
- FITS upload + HDU index (default 1).
- **JPG + sidecar-FITS** upload: JPG supplies pixels, sidecar FITS supplies WCS via its header. Validate that the JPG's `(NAXIS1, NAXIS2)` match the sidecar's image dimensions; if mismatch, refuse.
- Fiducial RA/Dec (decimal or sexagesimal).
- `apa_v3` slider + numeric.
- Display toggles: MSA outline / all shutters / quadrant filter.
- Stretch controls (linear/sqrt/log + percentile clip).

Main pane: Bokeh `figure` with image as `image_rgba` glyph + MSA outline (`MultiPolygons`) + compass + scale bar. Pan/zoom on by default.

Deliverable: drop a JADES F356W mosaic, see the MSA at any APA_V3. No interaction beyond pan/zoom yet.

### M2 — Target catalog overlay (0.5 day)

- `app/catalog.py`: load CSV / ASCII / FITS-table with `ID, RA, DEC` (+ optional cols per `CONTEXT.md`).
- Render markers on the image; size/color by `priority` / `mag` if present.
- Sidebar: catalog upload, marker style controls, label-visibility toggle.
- Table view below the figure listing the targets inside the current field, sortable.

### M3 — Operability mask (0.5 day)

- Fetch the latest `jwst_nirspec_msaoper_*.json` from CRDS once; store in `data/`.
- `app/operability.py`: `(operable, reason)` arrays of shape (4, 171, 365).
- Toggle: "Apply operability". Failed-open in red, failed-closed in gray, operable in default style.

### M4 — Hand-picking UI (2 days)

This is the planner core. (Already in Bokeh; no migration needed.)

- Two shutter `MultiPolygons` layers sharing a `ColumnDataSource` each:
  - **Background layer** (light gray, low alpha): all operable shutters in view.
  - **Open layer** (saturated, opaque): user-selected shutters. Appending/removing from this CDS is the core interaction.
- Click handlers via `TapTool` + a custom `CustomJS` or Python callback on the relevant CDS:
  - **Click a target marker** → propose the nearest operable shutter as the slitlet center, add a 3-shutter slitlet (s−1, s, s+1 in column d of quadrant q) to the open set, link it to the target ID. Slitlet height configurable (1 / 3 / 5).
  - **Click an empty shutter** → add it standalone (no target link; role = `manual`).
  - **Click an already-open shutter** → remove it (and its slitlet siblings if part of a group).
- Sidebar:
  - Slitlet height (1/3/5).
  - "Snap target to nearest operable" toggle.
  - Counter: "N open shutters covering M targets".
  - "Undo last action" and "Clear all".
- Conflict detection (M4 stretch, can defer to M5): no two open shutters in the same `(q, d)` column unless they belong to the same contiguous slitlet (because they'd disperse onto the same detector rows).

### M5 — Wavelength cutoffs (analytic model, no CRDS) (1 day) — **shipped with placeholder, needs replacement**

**Known issue (flagged by reviewer 2026-05-12)**: `V2_DISP_EXTENT = 180″` in `app/wavelengths.py` is a placeholder. The per-shutter wavelength shift reaches ±4 μm for PRISM and ±2 μm for G395M across the full MSA — far above the ±0.05 μm accuracy target. Tooltip wavelengths are correct at the fiducial only. Replace with per-grating `dλ/dV2` constants sourced from JDox before trusting the tooltip for non-fiducial shutters.



Per the user (2026-05-11): skip `jwst.assign_wcs`. Use an analytic per-grating dispersion model keyed off the shutter's V2/V3 (already in hand from `nirspec_msa_v2v3.npz`).

**Model.** NIRSpec disperses along V2. For each (disperser, filter):

- `λ_fid` = wavelength at the MSA fiducial (V2_ref, V3_ref).
- `dλ/dV2` ≈ const (μm / arcsec) — published per grating in JDox / NIRSpec dispersion solution tables.
- Detector boundaries in V2 (NRS1, gap, NRS2) are fixed offsets from the fiducial.
- Filter blue cutoff (e.g. F070LP > 0.70 μm, F100LP > 0.97 μm, F170LP > 1.66 μm, F290LP > 2.87 μm, CLEAR > 0.6 μm) is a hard lower bound.

For shutter at `(V2_sh, V3_sh)`:

```
λ_center  = λ_fid + (V2_sh − V2_ref) * dλ/dV2
λ_blue    = max(filter_blue_cutoff, λ_center − Δλ_to_NRS1_left)
λ_red     = λ_center + Δλ_to_NRS2_right
λ_gap_lo, λ_gap_hi = λ at fixed V2 positions of the NRS1/NRS2 boundary
```

A small per-grating table of `(λ_fid, dλ/dV2, NRS1_x_range, NRS2_x_range, filter_cutoff)` is enough.

**M5a — `app/wavelengths.py`** (in repo, no notebook needed):
- `GRATING_PARAMS` dict literal with the constants per (disperser, filter). Sourced from JDox (jwst-docs.stsci.edu → NIRSpec → MSA → Dispersers and Filters) and NIRSpec dispersion-solution PDFs.
- `cutoffs(v2_arcsec, v3_arcsec, disperser, filter) → {lam_blue, lam_gap_lo, lam_gap_hi, lam_red}` returns microns; returns `None` for fields where the shutter falls off-detector.
- Validation: cross-check 3–5 shutter positions per grating against the published JDox wavelength-range plots. Target accuracy ±0.05 μm.

**M5b — UI**:
- Sidebar dropdowns for disperser & filter.
- Hover tooltip on each shutter shows λ_blue / gap / λ_red.
- **Spectral-trace overlay**: for each open slitlet, draw a thin strip on the sky representing where the dispersed spectrum lands (~0.2″ wide × the V2-extent of NRS1+NRS2 in arcsec). Two open slitlets whose strips overlap → flag in red; that's a spectral conflict the user must resolve.

### M6 — Export & re-import (1.5 days)

Two distinct export paths (see `CONTEXT.md` → "Export formats").

**A. APT-loadable triple** (matches eMPT's `m_pick_output/`):

- `app/empt_io.py`:
  - `write_observed_targets_cat(path, targets_in_config)` — eMPT whitespace-separated catalog (`No, No_sub, No_cat, Pr, RA[deg], Dec[deg]`).
  - `write_pointing_summary_txt(path, pointing)` — free-form text with `Nod 0` central RA/Dec, `PA_AP`, `PA_V3`. (`Nod 1`/`Nod 2` echo `Nod 0` for our single-PA case.)
  - `write_shutter_mask_csv(path, open_shutters, operability)` — the 730 × 365 (TBC) grid of `{x, s, 1, 0}` cells in eMPT's exact tiling. **Validate by writing one and diffing against `trial_00/m_pick_output/pointing_100/shutter_mask.csv` from the eMPT repo** before declaring done. If the cell alphabet, tile order, or line endings don't match byte-for-byte, APT will silently mis-place shutters.
- A help panel in the UI showing the 3-step APT import workflow from `CONTEXT.md`.

**B. Internal session JSON** (round-trips inside our tool, not for APT):

- `app/session_io.py`:
  - `export_session_json(...)` writing the schema at the bottom of `CONTEXT.md`.
  - `import_session_json(text)` rebuilding open-shutter set + pointing + image/catalog paths.
  - Test: export, clear app state, import, assert the rendered open-shutter set is byte-identical.

Sidebar gets "Export → APT bundle", "Export → session JSON", "Import session JSON" buttons.

### M7 — Polish (open-ended)

- Save/load full session (image path, catalog, pointing, picks) as a workspace file.
- Export the figure as PDF/PNG with scale bar and metadata.
- Optional: PA sweep animation (vary APA_V3 over a range, render frames) to help pick a good visit-PA — even though each *run* is a single PA, exploring before committing is valuable.

---

## Open decisions resolved with user (2026-05-11)

1. "NIRSpec PA" = APA_V3. ✓
2. JPG mode = sidecar FITS provides WCS. ✓
3. One PA per session, no multi-PA stacking. ✓
4. End-state is MPT/eMPT-style hand-picking with APT-loadable export. ✓
5. **Stack**: Bokeh server from M1 onward. ✓
6. **APT ingestion**: emit eMPT's three files (`observed_targets.cat`, `pointing_summary.txt`, `shutter_mask.csv`). User has used eMPT before; this is the proven path. ✓

## Open decisions still to settle

- **Slitlet height default**: 3 (standard NIRSpec MOS for nod-and-shuffle) or user-set per pick? Default proposal: 3, with a slider to change.
- **Wavelength-cutoff source of truth**: `jwst.assign_wcs` (accurate, needs CRDS) vs parameterized JDox figures (approximate, no deps)? Recommend the former; only fall back if `stenv` can't reach CRDS.
- **`PA_AP` vs `PA_V3` in `pointing_summary.txt`**: APT distinguishes aperture PA from V3 PA. Our session-level `apa_v3_deg` maps to one or the other depending on which aperture APT thinks is in use. Resolve at M6 by reading what eMPT writes for both fields.

---

## Risks

- **PA sign / orientation**: the 138.5° MSA tilt and `spherical_offsets_by` order are subtle. Keep the notebook's exact arithmetic. Sanity-check M1 against an APT-exported visualization before declaring done.
- **WCS via sidecar FITS**: scaling/flipping discrepancies between the JPG and the sidecar are easy to introduce (e.g. JPG saved with origin-top-left vs FITS origin-bottom-left). Validate by overlaying a known-position marker (a Gaia star) before trusting an alignment.
- **APT round-trip**: until we verify which format APT accepts, "exports a config that loads in APT" is aspirational. Build the internal JSON first; add adapter(s) once the target format is confirmed.
- **Spectral-trace correctness**: a wrong trace overlay will silently let users pick spectrally overlapping slitlets. Validate M5b traces against a real MOS simulation (`jwst` pipeline `assign_wcs` output) at a handful of shutters before trusting the conflict detector.
- **Performance**: ~250 k shutters total. Cull aggressively to those in the current view bbox + 5″ margin. Bokeh handles ~10 k polygons interactively without trouble; matplotlib chokes above ~30 k. Don't render quadrants the user has filtered out.
- **pysiaf PRD mismatch** (`stenv` on PRDOPSSOC-068 vs online 072): not blocking. Upgrade if absolute positions need to be quoted at <0.05″.
