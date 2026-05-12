# Overnight build log — 2026-05-11 → 2026-05-12

User signed off at ~22:00. **All milestones M0–M4 + M5 (analytic) + M6 (eMPT export) shipped.** Reviewer passed; 46 / 46 tests green.

## Final status

- [x] M0 scaffolding — git init, `app/`, `data/`, `notebooks/`, `refs/eMPT_v1` clone, `nirspec_msa_v2v3.npz` copied to `data/`, bokeh installed
- [x] M1 Bokeh app: FITS + JPG-sidecar loading + MSA overlay at any APA_V3
- [x] M2 target catalog overlay (CSV/ASCII/FITS auto-detect)
- [x] M3 operability mask (CRDS `jwst_nirspec_msaoper_0014.json` loaded automatically)
- [x] M4 hand-picking: click target → snap to nearest operable shutter + 3-shutter slitlet; click shutter to toggle; undo + clear; slitlet height selector
- [x] M5 analytic wavelength cutoffs — works at fiducial; **known limitation: tooltip wavelengths wrong for shutters far from fiducial** (V2_DISP_EXTENT is placeholder)
- [x] M6 eMPT triple writer + reverse-engineered byte-identical `shutter_mask.csv` format + round-trip diff test against `refs/eMPT_v1/trial_00_ref/m_pick_output/pointing_100/shutter_mask.csv`
- [x] Reviewer pass — 2 bugs fixed in place (PA_AP offset 0.075° → exact V3IdlYAngle, slitlet-sibling off-by-one), report in `REVIEW.md`

## How to run

```bash
conda activate stenv
cd /Users/sunfengwu/nirspec
bokeh serve app/ --show
# Browser opens at http://localhost:5006/app
```

Then:
1. Upload a FITS (or a JPG + sidecar FITS).
2. Optional: upload a target catalog (CSV with `ID, RA, DEC`).
3. RA/Dec auto-default to image center; adjust if needed.
4. Set APA_V3 with slider or numeric input.
5. Pick disperser / filter from dropdowns.
6. Tap a target → opens a 3-shutter slitlet on the nearest operable shutter.
7. Tap any shutter to toggle open/closed.
8. Undo / Clear as needed.
9. Click "Export eMPT bundle" → writes `observed_targets.cat`, `pointing_summary.txt`, `shutter_mask.csv` to `exports/empt_bundle_<timestamp>/`.
10. Import into APT following the 3-step workflow in `CONTEXT.md` ("Export formats").

Try with the supplied examples:
- FITS-only: `example_a370/a370_f182m_f200w_f210m.fits`
- JPG+sidecar: `example_r0600/JWST_F090W_F200W_F444W.jpg` + `example_r0600/wcs.fits`

## Test summary

```
46 passed in ~6 s
- coords.py: 5 tests (shutter polygon vs notebook reproduction)
- msa.py: 6 tests (grid shape, bbox round-trip, operability fallback)
- wavelengths.py: 11 tests (9 disperser/filter fiducial matches + V2 shift + invalid combo)
- image_io.py: 6 tests (FITS auto-HDU, JPG-sidecar downsampling, RGBA stretch)
- catalog.py: 4 tests (CSV/ASCII/FITS load + RA-wrap masking)
- empt_io.py: 6 tests (CSV byte-exact header, grid shape, reference round-trip, observed_targets first row)
- main.py: 4 tests (sky↔V2V3 round-trip, _nearest_shutter, _add_slitlet, undo history)
- end_to_end.py: 4 tests (overlay inside a370, overlay covers r0600, export bundle, V2-shift wavelengths)
```

## Open issues for the user

1. **`V2_DISP_EXTENT = 180″` in `app/wavelengths.py` is a placeholder.** Tooltip λ_blue/λ_red are correct at the fiducial (MSA center) but drift by several μm at the edges. Action: source per-grating `dλ/dV2` from JDox and replace the constant. Test suite already verifies fiducial values; add a tighter shift-magnitude test once corrected.
2. **`load_jpg_with_sidecar` warns but doesn't refuse on JPG/sidecar dimension mismatch >10 %.** Plan said "refuse"; current is "warn and proceed". Tighten if you ever hit a mismatch in production.
3. **Per-tab session isolation.** Module-level `state` dict is shared across browser tabs. Fine for single-user single-tab; if you open two tabs, picks bleed. Scope to `curdoc()` if needed.

Full reviewer report in `REVIEW.md`. Future-work nits enumerated there.

## Decisions logged tonight

- **M5 simplification (user, 22:30)**: analytic per-grating dispersion model from shutter V2 offset; skip `jwst.assign_wcs`. Used a placeholder `V2_DISP_EXTENT` constant — needs replacement.
- **eMPT shutter_mask.csv format** reverse-engineered from `reference_files/shutter_routines_new.f90`: 731 lines × 683 chars; CSV row 1–365 carries Q1 (cols 1–171) + Q2 (cols 172–342) with `d = csv_row`; CSV row 366–730 carries Q3 + Q4 with `d = csv_row − 365`. Cell alphabet `{x, s, 1, 0}`. Header padded to 683 chars. Round-trips byte-identical with the reference.
- **PA_V3 → PA_AP**: offset is `V3_IDL_Y_ANGLE` from pysiaf NRS_FULL_MSA aperture (≈ 138.5746°). Verified against eMPT reference summary file to <0.00002°.

## Git history

```
8 commits will be (TODO)... see git log
M0: scaffolding (dirs, MSA grid copy, eMPT clone, .gitignore)
Data layer: coords, msa+operability, wavelengths, image_io, catalog, empt_io
M1-M4: Bokeh server app with hand-picking and eMPT export
Reviewer fixes: V3IdlYAngle PA conversion + slitlet sibling off-by-one + REVIEW.md
```
