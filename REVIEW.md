# Overnight review — 2026-05-12

## Test results

- **46 / 46 tests pass** (no failures, no skips reached, ~6 s).
- **Bokeh server**: started on port 5009, `curl` → HTTP 200. `/tmp/bokeh-review.log` shows clean boot: server up, CRDS msaoper JSON loaded, only the cosmetic `PRDOPSSOC-068 vs 072` pysiaf warning (already noted in CONTEXT.md as non-blocking). Server stopped after review.

## Correctness findings

### Critical (must fix before user uses the tool)

- **None.** All the load-bearing pieces (V2/V3 ↔ sky transforms, shutter-mask tiling, JPG/FITS orientation) verified correct.

### Important (should fix but not blocking)

- **`V2_DISP_EXTENT = 180″` placeholder gives wildly wrong wavelength tooltips for shutters away from the fiducial.** With the current constant and the actual V2 extent of the MSA (226″ → 534″, span 308″), the per-shutter wavelength shift in `app/wavelengths.py:cutoffs` reaches ±4 μm for PRISM/CLEAR, ±2 μm for G395M/F290LP, and ±0.4 μm even for the modest G140H/F070LP — far above the >0.5 μm threshold the brief asks to flag. The tooltip values are useful only at the fiducial; everywhere else they are nonsense (often clipped to the filter cutoff). The user already labelled this as placeholder in `PLAN.md`. Needs a real `dλ/dV2` per grating sourced from JDox before the tooltip is trusted. Module-level constant in `app/wavelengths.py:27`.
- **`pa_ap` formula was off by 0.075° (now fixed).** Was `(pa_v3 + 138.5) % 360`; the eMPT reference `pointing_summary.txt` shows PA_AP = 99.579041 for PA_V3 = 321.004456, which requires the offset to be `V3IdlYAngle = 138.5745697°` (from pysiaf `NRS_FULL_MSA`), not 138.5°. APT may or may not be tolerant of 0.075°; safer to be exact.
- **Off-by-one in `on_open_shutter_tap` sibling-removal radius (now fixed).** Slitlet-siblings filter used `abs(other.s - s) <= slitlet_height // 2 + 1`, which for height=3 catches s±2 (i.e. an adjacent unrelated slitlet sharing the same `target_id`). Tightened to `<= slitlet_height // 2`. No test broke; the existing tests don't exercise this branch.

### Nits / future work

- `_add_slitlet` silently skips failed shutters (`OPERABLE[q-1,s-1,d-1]` test) and just returns a smaller count. Status message says "{n} shutters opened" which is honest, but the user has no obvious cue that the slitlet was truncated. Future: surface "skipped k failed".
- `_nearest_shutter(..., require_operable=True)` returns the nearest operable shutter even if it's arcminutes away. If the target is off the MSA entirely the user will silently get a faraway shutter. Future: also threshold by separation and return `None` if > ~30″.
- `on_open_shutter_tap` calls `_push_history()` *before* the in-set check; if the tap fires for a stale index that's no longer in `state["open_shutters"]`, history grows by an empty snapshot. Minor — only matters if undo behavior surprises.
- `on_pa_slider` ↔ `on_pa_text` are protected from infinite loops by Bokeh's "only fire on real change" semantics, but the slider value is formatted to `.2f` which forces 0.01° quantization on round-trip. Fine for hand-picking, not for precision use.
- `load_jpg_with_sidecar` warns but does not refuse when JPG dims disagree with sidecar `NAXIS{1,2}` by >10 %. `PLAN.md` says "refuse"; current behavior is "warn and proceed". A real-world mismatch would silently mis-place overlays.
- `observed_targets.cat` rows come out in insertion order (i.e. click order). The eMPT reference is sorted by `No_cat`. APT may not care, but worth sorting for tidy diffs.
- `state` is a single module-level dict shared across all sessions of the Bokeh server. For one user / one tab this is fine and matches Bokeh's per-session callback threading. If the user ever opens two tabs, picks bleed across them. Document or scope to `curdoc()` if multi-session is ever expected.
- `_nearest_shutter` brute-forces over all 249,660 shutters every tap (~1 ms on this machine, fine), but allocates two full grid-shaped arrays each call. Could be cached. Not urgent.
- Operability JSON: `msa.py` interprets `"open"` as failed-open and `"closed"` as failed-closed via substring matching. CRDS schema seems to use phrases like "stuck open"/"stuck closed", which are caught. Anything containing only "operable" wouldn't match — would silently leave them operable, which is the safe default.

## Fixes applied in this review

- `app/main.py:677` — slitlet-sibling removal radius `<= slitlet_height // 2 + 1` → `<= slitlet_height // 2`. Off-by-one bug that would, on un-pick, also delete shutters 2 rows away if they happened to share `target_id`.
- `app/coords.py:14-17` — added `V3_IDL_Y_ANGLE` from pysiaf (`= 138.5745697°` on PRDOPSSOC-068) so PA_V3 → PA_AP conversion in the export is correct, not just close.
- `app/main.py:42-49` and `app/main.py:752-754` — import `V3_IDL_Y_ANGLE` and use it in `on_export` instead of the literal `138.5`. Verified against the eMPT reference: `(321.004456 + 138.5745697) % 360 = 99.579025` vs reference `99.579041` (agreement to <0.00002°, limited by sigfigs in PRD).

Full test suite re-run after fixes: 46 / 46 pass.

## Verified to match the notebook

- `app/coords.py:rot_matrix` ← `footprint_emerald.ipynb` `rot_matrix` (cell containing `def rot_matrix(rotation = 30)`): identical `[[c,-s],[s,c]]`.
- `app/coords.py:shutter_corners_v2v3` ← notebook `shutter_corners_v2v3`: identical, including the `rot_matrix(138.5)` MSA tilt and the LL/LR/UR/UL corner order.
- `app/coords.py:v2v3_to_radec` ← notebook `v2v3_to_radec`: identical — `corners - [V2Ref,V3Ref]`, then `@ rot_matrix(pa_v3)`, then `spherical_offsets_by(Δlon=east, Δlat=north)`. Sign convention preserved.
- MSA reference: `MSA_V2_REF = 378.563202″, MSA_V3_REF = -428.402832″` matches CONTEXT.md’s "≈ (378.563, −428.403)".

## Verified vs. the eMPT reference

- Header line of our `shutter_mask.csv` matches `refs/eMPT_v1/trial_00_ref/m_pick_output/pointing_100/shutter_mask.csv` byte-for-byte (683 chars, exact text).
- All-operable line count = 731 = ref. All line widths = 683 = ref.
- Parsed the reference CSV, fed its `(operable, reason, opens)` back into the writer, diffed re-written vs reference: **0 / 731 lines differ.** Round-trip is byte-identical.
- Per-quadrant failed-shutter counts from the reference (15520/14517/15637/18115 'x' and 6/3/12/1 's') are exactly the values pinned in `test_operability_roundtrip_against_reference`.
- Cross-checked the Fortran tiling at `refs/eMPT_v1/reference_files/shutter_routines_new.f90:581-672`: top half iterates `ir=1..365` writing `kk=1` then `kk=2`, inner `jj=1..171`; bottom half `ir=366..730` writing `kk=3` then `kk=4`. Our writer iterates `d=1..365` with `cells[0,:,d-1]` then `cells[1,:,d-1]` etc. — matches exactly.
- `_sky_to_v2v3` inverse algebra: forward is `(Δlon, Δlat) = (Δv2·c + Δv3·s, -Δv2·s + Δv3·c)`; inverse must give `Δv2 = Δlon·c - Δlat·s, Δv3 = Δlon·s + Δlat·c`. Code computes exactly that via `rot_matrix(-pa_v3)`. Existing round-trip test in `tests/test_main_handpicking.py:test_sky_to_v2v3_round_trip` is a real check: it picks a shutter, forward-transforms, inverse-transforms, and asserts both V2 and V3 recovered to 1e-3 arcsec — not tautological.

## Sign-off

The tool is in a usable state for hand-picking and for producing an APT-loadable export bundle. The transforms are verbatim from the notebook, the shutter-mask CSV is byte-identical to eMPT on a round-trip, and Bokeh boots cleanly. The biggest caveat is `V2_DISP_EXTENT = 180″`: wavelength tooltips are correct at the MSA fiducial only — for shutters off-fiducial they are off by up to several μm and should not be quoted to anyone. Two small bugs (sibling-removal off-by-one, PA_AP offset) were fixed in place; all 46 tests still pass.
