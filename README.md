# NIRSpec MSA hand-picking planner

A Bokeh-based planner for laying out NIRSpec MSA shutter configurations
on an image of the target field. Mirrors the MPT/eMPT workflow but lets
you pick shutters by hand instead of by automated optimization, and
exports an APT-loadable bundle.

## Run

```bash
cd /Users/sunfengwu/nirspec
./run.sh
# or manually:
conda activate stenv
bokeh serve app/ --websocket-max-message-size 524288000 --show
```

The browser opens at `http://localhost:5006/app`.

## Sidebar workflow

### 1. Load an image (Image section)

Two paths supported. Always prefer the **path** inputs over the upload
widgets — uploads are limited to ~20 MB by Bokeh's WebSocket protocol
(the `run.sh` wrapper bumps this to 500 MB but typing a path is still
faster).

- **FITS**: paste an absolute path to a 2D-image FITS. The first HDU
  with image data is auto-selected. The WCS is read from that HDU's
  header. NaNs are zeroed for display.
- **JPG + sidecar FITS**: paste paths to a JPG (or PNG) and a FITS file
  whose header carries the WCS. The sidecar may be header-only
  (`NAXIS=0`); pixel dimensions are inferred from the JPG. Images
  larger than 6000 px on a side are downsampled with WCS scaling.

### 2. Optional: target catalog (Catalog section)

Accepts CSV, whitespace-ASCII, or FITS-table. Required columns
(case-insensitive): `ID`, `RA`, `DEC`. Optional: `priority`/`Pr`,
`mag`/`F444W_mag`, `z`/`zspec`/`zphot`, `label`/`name`. Targets are
drawn as yellow circles, clipped to the visible image.

### 3. Set the pointing (Pointing section)

- **Pointing RA / Pointing Dec** — fiducial sky position the MSA centers
  on. Auto-filled with the image center when an image loads.
- **V3 PA** — the position angle of the JWST V3 axis on sky. This is
  what drives the V2/V3 → RA/Dec transform.
- **NIRSpec APA** — the *aperture* PA of NIRSpec, equal to
  `V3PA + V3IdlYAngle (mod 360)` where `V3IdlYAngle ≈ 138.575°` for
  `NRS_FULL_MSA`. This is what APT/MPT calls "NIRSpec PA". Editing one
  PA field syncs the other automatically.
- **Visibility date** + **Compute allowed V3 PA** — queries
  `jwst_gtvt`'s ephemeris and reports the allowed V3PA window for the
  requested date. First call takes ~5–8 s.

The lime cross on the figure is the **pointing center marker**.
**Shift-click anywhere on the image** to move the pointing center to
that sky location; RA/Dec inputs update and the MSA overlay follows.

### 4. Display layers

Toggle independently:
- MSA outline (dodgerblue quadrant rectangles)
- Operable shutters (faint white grid)
- Apply operability (hides failed shutters in the bg layer; stuck-open
  are always shown regardless)
- Targets (yellow circles)

### 5. Hand-pick shutters

The planner uses a "snap-to-nearest" model so you don't have to land
precisely on a tiny polygon:

- **Single click** on the image → snaps to the nearest shutter and
  toggles it open. If you click near a yellow target, opens a
  3-shutter slitlet centered on the nearest operable shutter to that
  target.
- **Double click** on the image → toggles a cyan highlight on the
  nearest shutter. Highlight is a visual flag only; it is *not*
  exported.

Layer colors:

| Color | Meaning |
|---|---|
| Faint white | Operable shutter |
| Red edge, faint red fill | Stuck-open shutter (always visible) |
| Red filled | Open (user-selected) shutter |
| Cyan edge | Highlighted shutter |
| **Orange tint** | Shutter that would have **spectral overlap** with an open shutter (shares its `s` row in the same quadrant) — open it only if you accept the conflict |
| Gold | NIRSpec fixed slits (S200A1/A2, S400A1, S1600A1, S200B1) |
| Lime cross | Draggable pointing handle |

Slitlet height (1, 3, or 5) and "snap target to nearest operable"
toggle are in the Display section.

**Undo last** / **Clear open** revert recent changes.

### 6. Disperser / filter (Instrument section)

Drives the wavelength tooltip on open shutters (λ_blue / λ_red in μm).
**Caveat**: the analytic wavelength model is calibrated at the MSA
fiducial only; tooltip values for shutters far from center are
approximate (off by up to several μm for PRISM).

### 7. Export (Export section)

Click **Export eMPT bundle** to write three files into a timestamped
subdirectory of the export dir:

1. `observed_targets.cat` — whitespace-separated source catalog. Load
   into APT via *Form Editor → Targets → Import MSA Source Catalog*.
2. `pointing_summary.txt` — RA/Dec/PA values to copy-paste into APT's
   *MSA Planner → Search Grid* panel (set search box width/height to 0
   and Number of configurations to 1, then Generate Plan).
3. `shutter_mask.csv` — the 730 × 342 MSA shutter grid. For each nod
   row in the generated plan, *Edit Configuration → Edit → Import CSV*
   to overwrite APT's auto-generated mask with this one.

Format is byte-compatible with eMPT's outputs; the writer was reverse
engineered from `reference_files/shutter_routines_new.f90` and
round-trips identical against `trial_00_ref/m_pick_output/pointing_100/`
in the eMPT repo.

## Toolbar interactions

- **Wheel**: zoom both axes equally (locked aspect, scrolling on an
  axis does not zoom that axis alone).
- **Drag**: pan the view.
- **Box zoom**: select box-zoom icon, then drag a rectangle.
- **Reset**: toolbar reset icon.
- **Click** (no modifier): snap to nearest shutter / target (see above).
- **Shift-click**: move pointing center to the click location.
- **Double-click**: toggle cyan highlight on the nearest shutter.

## Architecture

- `app/main.py` — Bokeh server entry, all UI wiring.
- `app/coords.py` — V2/V3 ↔ sky transforms ported verbatim from
  `footprint_emerald.ipynb`. The 138.5° MSA tilt and pysiaf
  `NRS_FULL_MSA` parameters are load-bearing — don't reorder.
- `app/msa.py` — loads the 4×171×365 shutter grid and CRDS operability.
- `app/wavelengths.py` — analytic per-grating dispersion model.
- `app/image_io.py` — FITS + JPG-with-sidecar loaders.
- `app/catalog.py` — CSV/ASCII/FITS catalog reader.
- `app/empt_io.py` — three eMPT writers.

## Known limitations / TODOs

- `V2_DISP_EXTENT = 180″` in `app/wavelengths.py` is a placeholder.
  Per-shutter λ tooltips are correct at the MSA fiducial but drift
  by several μm at the edges; needs JDox-sourced `dλ/dV2` constants.
- Bokeh sessions share the module-level `state` dict — opening two
  browser tabs lets picks bleed across them. Fine for single-user/single-tab.
- `load_jpg_with_sidecar` warns but doesn't refuse on JPG/sidecar
  dimension mismatch >10%.
- The wavelength model assumes uniform pixel scales in the WCS;
  uneven `cdelt1`/`cdelt2` would mis-place overlays.

## Tests

`pytest tests/` — 46 unit and end-to-end tests covering coord
transforms, shutter-mask CSV format (byte-diff against eMPT reference),
wavelength model, image loaders, catalog parser, and overlay rendering
on both example FITS and JPG fields.

## References

- JWST PA conventions: [JDox PA reference](https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-position-angles-ranges-and-offsets).
- MSA operability: STScI CRDS `jwst_nirspec_msaoper_*.json`.
- eMPT (the inspiration for our export format): Bonaventura et al.
  2023, A&A 672, A40 (arXiv:2302.10957); code at
  https://github.com/esdc-esac-esa-int/eMPT_v1.
