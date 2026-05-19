# vMPT — Project Context

Cold-start reference for the interactive NIRSpec MSA planning tool.
Captures the decisions an agent (human or otherwise) needs to know
before touching the code — invariants, file roles, the bundle layout,
key formulas, and the gotchas we've already paid for in bug-hunts.

For end-user docs (install, two-minute tour, troubleshooting) see
[`README.md`](README.md).

---

## What the tool does

1. Loads a JWST-field image (FITS with WCS, or JPG + sidecar FITS).
2. Optionally loads a target catalog (CSV / ASCII / FITS) and overlays
   marker circles on the image.
3. Overlays the NIRSpec MSA at a user-chosen pointing (RA, Dec, V3 PA),
   the 5 fixed slits, and the CRDS operability mask.
4. Lets the user hand-pick **N-shutter slitlets** (N ∈ {1, 2, 3, 5})
   by clicking the image. The clicked shutter snaps to the nearest
   operable shutter; clicking an open one closes it + siblings.
5. Computes live **spectral-overlap warnings** based on the chosen
   disperser/filter and the open + stuck-open shutter set.
6. Exports a bundle (`Save session` / `Export eMPT bundle`) that
   round-trips back into vMPT AND loads into APT MPT.
7. Imports plans from APT MPT JSON, shutter-mask CSVs, or `.aptx`
   archives (local or fetched directly from STScI by program ID).

---

## Architecture (live, current state)

```
app/
├── main.py            Bokeh server entry; UI wiring; refresh_overlays
├── coords.py          V2/V3 ↔ RA/Dec transforms (pysiaf-backed)
├── msa.py             MSA shutter grid + CRDS operability loader
├── wavelengths.py     Per-grating dispersion model + cutoffs
├── image_io.py        FITS + JPG-with-sidecar loaders (LoadedImage)
├── catalog.py         Catalog reader → Catalog dataclass
├── empt_io.py         eMPT-format writers + MPT-importable .cat writer
├── session_io.py      Bundle save/load (MPT plan + workspace sidecar)
├── mpt_io.py          APT MPT JSON parser + .aptx archive reader
├── static/favicon.svg vMPT favicon (MSA grid + pointing cross)
└── templates/index.html  Injects favicon via data URI

data/
└── nirspec_msa_v2v3.npz   (4 × 171 × 365) shutter V2/V3 arcsec

tests/                 pytest suite (60+ tests; ~7 s)
example_a370/          Abell 370 FITS (44 MB)
example_r0600/         RXCJ0600 JPG + sidecar (240 MB)
exports/               default output dir
```

---

## Coordinate plumbing (do not re-derive)

### Per-shutter V2/V3 grid

`data/nirspec_msa_v2v3.npz` carries two `(4, 171, 365)` arrays in
arcsec. Indexing convention (matches APT MPT JSON `slitlets`):

- axis 0: **quadrant q − 1**  (q ∈ {1,2,3,4})
- axis 1: **shutter row s − 1**  (s ∈ {1…171})
- axis 2: **shutter column d − 1**  (d ∈ {1…365})

Quadrant V2/V3 bounding boxes (sanity check):

| Q | V2 range | V3 range |
|---|---|---|
| Q1 | [+399.55, +533.73] | [−503.04, −370.62] |
| Q2 | [+316.00, +448.20] | [−406.62, −275.99] |
| Q3 | [+309.20, +442.30] | [−583.70, −450.22] |
| Q4 | [+226.02, +357.41] | [−486.74, −355.33] |

Detector pairing (used by spec-overlap calc):
- **NRS1** images **Q1 + Q3**
- **NRS2** images **Q2 + Q4**
- Cross-quadrant overlap only happens within these pairs.

Shutter pitch: **0.20″ along V2** (within a row), **~0.40″ along V3**.
The open aperture per shutter is **0.20″ × 0.46″**. The MSA is
rotated **138.5°** within the V2/V3 frame
(`coords.V3_IDL_Y_ANGLE ≈ 138.5746°`).

### Critical relations

- **APA = V3 PA + V3IdlYAngle (mod 360)**  — the NIRSpec aperture PA
  vs. the V3 axis PA. Used everywhere we cross between APT inputs
  and vMPT's V3-PA state. Stored as `state["pa_v3"]`.
- **`v2v3_to_radec`** (`coords.py:33`): rotates V2/V3 offsets into the
  sky tangent plane via `rot_matrix(pa_v3)`, then applies
  `SkyCoord.spherical_offsets_by`.
- **`_sky_to_v2v3`** (`main.py:1331`): inverse, using
  `SkyCoord.spherical_offsets_to` which has the `cos(dec)` correction
  built in. We fixed a +1/cos(dec) bug here for Dec≠0 fields
  (Dec=−20° was off by 6 % = ~140 px).

### Loading the WCS

- FITS path → `astropy.wcs.WCS(header).celestial`.
- JPG + sidecar → `WCS` from the sidecar header, then rescaled to the
  downsampled JPG dimensions. The sidecar lives next to the JPG;
  vMPT auto-discovers it when loading a session that points only at
  the JPG.

---

## Operability mask

Loaded once at startup from CRDS via `app/msa.py`. Exposes:

- `OPERABLE` — `(4, 171, 365)` bool. True if the shutter is commandable.
- `REASON`  — `(4, 171, 365)` int8. `0`=operable, `1`=failed-closed,
              `2`=failed-open (stuck open).

Flattened views (`_FLAT_REASON`) are used in `refresh_overlays` for
fast vectorised masking.

There are **22 stuck-open shutters** in the current CRDS reference
(`jwst_nirspec_msaoper_0014.json` or similar). Stuck-opens always
disperse light — they contribute to the spec-overlap calculation
even if the user hasn't opened them.

---

## State (single global `state` dict in `main.py`)

| Key | Meaning |
|---|---|
| `image` | `LoadedImage` or `None`. Carries `data`, `wcs`, `shape`, `source_path`, `mode`, and (JPG-only) `wcs_sidecar_path`. |
| `catalog` | `Catalog` or `None`. RA/Dec/IDs/priority/mag/z arrays + `source_path`. |
| `ra_deg`, `dec_deg` | Pointing center. |
| `pa_v3` | V3 PA in degrees (mod 360). APA = `pa_v3 + V3_IDL_Y_ANGLE`. |
| `disperser`, `filter` | e.g. `"PRISM"`, `"CLEAR"`. |
| `slitlet_height` | N ∈ {1, 2, 3, 5}. Determines `_slitlet_offsets` and toggle-off siblings. |
| `open_shutters` | `dict[(q, s, d) → OpenShutter]`. The user's picks. |
| `highlighted` | `set[(q, s, d)]` — cyan-edge visual flag, not exported. |
| `history` | Undo stack (capped at 50 snapshots of `open_shutters`). |
| `shutter_to_catids` | `dict[(q, s, d) → [source_id, …]]` — catalog sources whose footprint lands in each shutter. Rebuilt on pointing / PA / catalog change. |
| `snap_to_operable` | When a target click misses, snap to the nearest operable shutter. |

---

## Render pipeline (`refresh_overlays`)

Heavy path, runs on every pointing / PA / disperser / open-shutter
change. Order matters — later layers read masks computed earlier:

1. **MSA outline** (4 dodgerblue quadrant rectangles).
2. **In-view mask** — `_in_view_mask`: projects all 249,660 shutter
   centres to image pixels and masks to the figure's current
   `x_range` × `y_range` bbox (post-zoom). Reused by all per-shutter
   layers.
3. **Stuck-open** (`REASON == 2`, always visible) — thick dark-red
   outline (`#b30000`, 2.5 px), light red fill (α=0.15).
4. **Spectral-overlap** — see below.
5. **Operable, unaffected** — silver edge (α=0.20), no fill. Filtered
   to `REASON == 0 ∧ in_view ∧ NOT in open_shutters ∧ NOT in
   overlap_idx`. Capped at `MAX_OPERABLE_RENDER = 10000`; above the
   cap the layer is blanked rather than stride-sampled.
6. **Open shutters** — red fill (`#ff8888`, α=0.35) + thicker red edge.
7. **Highlighted** — cyan-edge overlay.
8. **Fixed slits** — gold polygons.
9. **Pointing handle** — lime cross.
10. **Catalog targets** — yellow circles (if layer toggled on).

### Spec-overlap calculation (the subtle one)

For each dispersion source (every user-open + every stuck-open
shutter), the set of overlap-affected operable shutters is:

```
(s ∈ open ± SHVAL_S_TOLERANCE)   AND
(q ∈ same_detector_half)         AND  // NRS1={Q1,Q3}; NRS2={Q2,Q4}
(|ΔV2| < v2_overlap_distance(disperser, filter))
```

Current values:
- `SHVAL_S_TOLERANCE = 1` — only the open shutter's row + immediate
  neighbours above/below disperse onto overlapping detector pixels
  (matches eMPT's `shval ≈ s` exactly).
- `v2_overlap_distance`: PRISM = 35″, M-gratings = 200″,
  H-gratings = 500″ (`wavelengths.SPECTRUM_V2_HALFEXTENT`).
- Detector pairing prevents H-gratings (500″ window covers most of
  the MSA in V2) from spuriously lighting up Q1↔Q2 or Q3↔Q4.

The overlap polygons are drawn fill-only (orange, α=0.10, **no
edge**) so multiple dispersion sources stack and the colour
intensifies where many spectra overlap.

---

## Slitlet picking semantics

`_slitlet_offsets(N)` returns the relative-s offsets from the
clicked shutter:

| N | offsets | role layout |
|---|---|---|
| 1 | `[0]` | clicked = target |
| 2 | `[-1, 0]` | s=clicked-1 → sky; s=clicked → target |
| 3 | `[-1, 0, +1]` | sky / **target** / sky |
| 5 | `[-2,-1,0,+1,+2]` | sky / sky / **target** / sky / sky |

The clicked shutter is always at offset 0 ("target" role).
`_add_slitlet` skips offsets that fall outside `s ∈ [1, 171]` or land
on a failed-closed shutter.

If the slitlet has no caller-supplied `target_id`, `_add_slitlet`
looks up `state["shutter_to_catids"]` for each opened shutter and
adopts the first catalog source ID it finds. All shutters in the
slitlet share that ID. The status bar surfaces the auto-match.

---

## Bundle output (six files)

Naming convention: prefix by role.

```
<export_dir>/
├── MPT_plan.json                   ← APT MPT plan
├── <catalog_stem>.cat              ← APT-importable Target List
├── vMPT_workspace.json             ← vMPT-only state (paths, target_id, roles)
├── eMPT_observed_targets.cat       ← eMPT pipeline input
├── eMPT_pointing_summary.txt       ← eMPT pipeline input
└── eMPT_shutter_mask.csv           ← eMPT pipeline input
```

Filenames are constants in [`app/session_io.py`](app/session_io.py)
(`MPT_PLAN_FILENAME`, `MPT_CATALOG_FILENAME`, `WORKSPACE_FILENAME`,
`EMPT_*`). Pre-1.4 names (`session_MPT_plan.json`,
`vmpt_workspace.json`) are still recognised on load.

### `MPT_plan.json` — APT MPT plan

Mirrors the reference plan structure field-for-field
(verified against `exports/empt_bundle_*/G395H_F290LP.json` etc.):

```
instrument, name, aperturePA, theta, catalog, referencePointing,
configs[{name, version, info, masterBackground, slitlets, exposures,
         primaryIds, fillerIds}],
stats, errors,
plannerSpecification{gratingSpecification, planName, planAngle,
                     theta, candidates, slitSpecification,
                     searchParameters{spectralOverlapShutterOffsetMap, …},
                     slitSearchSpecification, maskingSpecification,
                     pointingSpecification, searchGridSpecification,
                     wavelengthRangeSpecification}
```

Key encoded values:

- `aperturePA` = `pa_v3 + V3_IDL_Y_ANGLE` (mod 360).
- `configs[0].slitlets` = `_group_into_slitlets(open_shutters)` —
  consecutive-s runs at fixed `(q, d)` collapsed into `{q, d, s, h}`.
- `configs[0].primaryIds` = list of catalog target_ids, in slitlet
  order (slitlets WITH a target_id come first; slitlets without are
  appended after — so positional `primaryIds[j] ↔ slitlets[j]`
  alignment holds for the targeted prefix).
- `configs[0].exposures[0].gratingFilter` = `"{disperser}_{filter}"`,
  e.g. `"PRISM_CLEAR"`.
- `configs[0].exposures[0].msaSlitlet` = `"ONE_SHUTTER"` /
  `"TWO_SHUTTER"` / `"THREE_SHUTTER"` / `"FIVE_SHUTTER"`.
- `catalog.name` = `catalog.primariesName` =
  `plannerSpecification.candidates.{primaries, catalog}` =
  **the stem of the `.cat` file in the bundle**. So importing the
  `.cat` under its default name on the APT side, then loading the
  plan, binds them automatically.
- `plannerSpecification.searchParameters.spectralOverlapShutterOffsetMap`
  = `"JWST_NIRSPEC_<DISPERSER>"`, e.g. `"JWST_NIRSPEC_G395H"`.

### `<catalog_stem>.cat` — MPT-importable target list

ASCII, **tab-separated**, `#`-prefixed header. Column names match
the [JDox MPT Catalogs spec](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-catalogs)
exactly:

| Column | JDox role | vMPT use |
|---|---|---|
| `ID` | integer source id | `No_cat` |
| `RA` | decimal degrees | (no `[deg]` suffix — APT matches the bare token) |
| `DEC` | decimal degrees | |
| `Weight` | numeric priority weight | from input catalog's `priority`/`Pr` (default 1; synth = 5) |
| `Primary` | Number-typed flag | 1 = primary, 0 = filler — user can edit to split the list |
| `Label` | free text (recognized label) | **`real`** for input-catalog matches, **`vMPT_synth`** for entries we synthesized at slitlet centres |

```
# ID	RA	DEC	Weight	Primary	Label
1	39.9826125000	-1.5916444000	1	1	real
2	39.9870166600	-1.5891333000	5	1	vMPT_synth
…
```

Downstream the user can filter on `Label == "vMPT_synth"` inside APT
(or in any text editor) to see which rows weren't in their input
catalog and decide whether to keep them.

### Shutter ↔ source matching ("Unconstrained")

`_rebuild_shutter_catalog_index` matches each catalog source to the
nearest operable shutter, allowing the source centre to sit **anywhere
inside the full MSA shutter pitch** (≈0.27″ × 0.53″) — not just inside
the narrower open aperture (0.20″ × 0.46″). This mirrors APT's
[*Unconstrained* Source Centering Constraint](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-planner):
a source whose centre falls behind a bar still matches the
neighbouring shutter, because the geometry forces it to be inside one
shutter's Voronoi cell.

Implementation: half-pitch box `(SHUTTER_HALF_PITCH_V2,
SHUTTER_HALF_PITCH_V3) = (0.135″, 0.265″)`. Sources outside that box
must be entirely off the MSA grid; they're discarded.

Filename basename: if the user loaded a catalog file the stem is
`Path(catalog.source_path).stem`; otherwise `MPT_catalog`. Either
way it matches `catalog.name` in `MPT_plan.json`.

The rows include:
1. **Real catalog sources** tied to user picks (look up via
   `state["shutter_to_catids"]` at export time).
2. **Fake entries** synthesised at the centre of any open-shutter
   slitlet with no real source. Priority `5`; ID is the smallest
   unused integer. Stuck-open shutters never get faked entries.

### `vMPT_workspace.json` — vMPT-only state

Carries everything `MPT_plan.json` can't: per-shutter `target_id` +
`role`, highlighted set, image / sidecar / catalog paths, slitlet
height, exact `pa_v3_deg` (the plan only carries `aperturePA`).
Lossless round-trip for the picking session.

### `eMPT_*.cat / .txt / .csv`

Same writers vMPT has shipped since M5 — formats reverse-engineered
from `refs/eMPT_v1/reference_files/shutter_routines_new.f90` and
matching the `trial_00/` example byte-for-byte for `shutter_mask.csv`.

### Loading the bundle

`import_session_json` accepts EITHER `MPT_plan.json` or
`vMPT_workspace.json` — the sibling is auto-discovered. Also
supports two legacy schemas:
- pre-1.4 filenames `session_MPT_plan.json` + `vmpt_workspace.json`
- pre-1.1 single-file format (flat top-level `open_shutters`,
  `pointing` block).

Image is routed by extension: `.fits/.fit/.fts` → `load_fits`;
`.jpg/.jpeg/.png` → `load_jpg_with_sidecar` (auto-finds the WCS
FITS in the same directory if the workspace doesn't name one).
If the recorded `image_path` is no longer on disk, the picks /
pointing still restore — the status bar tells the user to load an
image manually.

---

## Wavelength dispersion model

`app/wavelengths.py`. Linear-shift model per grating:

```
λ(V2) = λ_min + (V2 - MSA_V2_REF) × (λ_max - λ_min) / V2_DISP_EXTENT
```

with `V2_DISP_EXTENT = 180″` and per-grating `(λ_min, λ_max)` from
JDox. Endpoints are **clamped** to the grating's intrinsic range —
prevents the old "PRISM shows λ > 5.3 µm" bug.

`v2_overlap_distance(disperser, filter)` returns the V2 half-extent
used by the spec-overlap calculation (35″ PRISM, 200″ M, 500″ H).

---

## Importing APT plans (MPT tab)

`app/mpt_io.py`:

- `parse_mpt_json(path)` — top-level → list of `MPTPlan`. Picks the
  first **dispersed** exposure (skips `gratingFilter: null` target
  acquisition / imaging steps) for grating + RA/Dec. Falls back to
  `exposures[0]` for plans with no dispersed step (e.g. shutter-mask
  preview configs).
- `list_mpt_plans_in_aptx(aptx_path)` / `parse_mpt_json_in_aptx` —
  treats `.aptx` as a zip archive; finds embedded JSONs with both
  `configs` and `aperturePA` keys.
- `download_apt_program(program_id)` — fetches
  `https://www.stsci.edu/jwst-program-info/download/jwst/apt/<pid>/`,
  verifies zip magic bytes, writes to a tempfile. Some programs
  return 404 (unreleased).
- `parse_shutter_csv(path)` — open-mask CSV (730×342) → flat list of
  `OpenShutter`. Tiling: rows 1–365 carry Q1+Q2, rows 366–730 carry
  Q3+Q4; cols 1–171 carry Q1/Q3 s-indices, cols 172–342 carry Q2/Q4.

`_apply_plan` in `main.py` writes the parsed plan into the live UI:
RA/Dec/PA/disperser-filter/open_shutters. If no image is loaded yet,
emits a warning ("Plan loaded; load an image to see overlay") but
the state still applies — and a subsequently-loaded image keeps the
plan's pointing instead of recentering.

---

## Loading-overlay UX

`_show_loading(msg)` / `_hide_loading()` toggle a full-page spinner
overlay. The widget is a zero-size Bokeh Div whose inner HTML uses
`position: fixed` to escape the layout and cover the whole viewport.
84-px amber ring on a 45 % black backdrop with 2 px blur; CSS
`@keyframes vmpt-spin` rotates the ring every 0.9 s; 180 ms fade-in.

Triggered by: image loads, JPG+sidecar loads, catalog loads,
.aptx fetches, jwst_gtvt queries, **and** every pointing / V3 PA /
APA / disperser change. The change handlers defer the actual
recompute to the next document tick (`_deferred(_do)`) so the
spinner paints before the heavy work starts; a `finally:
_hide_loading()` clause guarantees it disappears even on error.
60 s safety timeout backstops broken callsites.

---

## Tests

```
pytest tests/    # 60+ tests, ~7 s
```

Notable coverage:
- `tests/test_session_io.py` — round-trip, legacy-schema fallback,
  workspace sidecar pairing, MPT-side load (`parse_mpt_json` on
  exported `MPT_plan.json`), "no file paths in MPT plan".
- `tests/test_mpt_io.py` — APT MPT JSON parse, `.aptx` round-trip,
  shutter CSV matches JSON unfolding.
- `tests/test_empt_io.py` — eMPT shutter-mask byte compatibility,
  observed-targets formatting, MPT-catalog writer header/data.
- `tests/test_end_to_end.py` — full bundle write + parse cycle.
- `tests/test_wavelengths.py` — fiducial wavelengths match published
  values, gap behaviour, clamping.

---

## Environment

- Activate with `conda activate stenv` (path: `/Users/sunfengwu/anaconda3`).
  `stenv` ships `astropy`, `pysiaf`, `jwst`, `numpy`. Add `bokeh` +
  `jwst_gtvt` with pip.
- `pysiaf` PRD note: stenv may report "PRDOPSSOC-068 doesn't match
  online PRDOPSSOC-072". The 138.5° MSA tilt and V2/V3 reference
  values are stable across these PRDs; safe to ignore unless
  precision <0.05″ matters.

---

## Bug-history (don't pay these costs twice)

- **`+1/cos(dec)` Jacobian bug** (Dec≠0 fields off by 6 % at Dec=−20°).
  Fixed by using `SkyCoord.spherical_offsets_to` instead of
  hand-rolled offsets in `_sky_to_v2v3`. Commit `4682571`.
- **PRISM showing λ > 5.3 µm**: linear shift model didn't clamp to
  the grating range. Fixed in `b96f126` — `wavelengths.cutoffs`
  clamps to `[lam_min, lam_max]`.
- **TapTool fading non-selected open shutters to 20 % alpha**:
  Bokeh's default nonselection rendering. Fix: drop the `TapTool`
  entirely, keep `fig.on_event(Tap, on_tap)`.
- **Spurious cross-quadrant orange tints** with H-gratings: the
  500″ V2 window passed for any pair of quadrants regardless of
  which detector they image. Fix: enforce NRS1={Q1,Q3} /
  NRS2={Q2,Q4} pairing in `refresh_overlays`.
- **`parse_mpt_json` missed grating when `exposures[0]` was target
  acquisition**: configs like `step1 copy` had `gratingFilter: null`
  on the first 1–2 exposures, the dispersed step came later. Fix:
  scan `exposures[]` for the first non-null `gratingFilter`.
- **`_set_image_and_recenter` overwrote APT-plan pointing**: loading
  an image after a plan clobbered the plan's RA/Dec. Fix: preserve
  RA/Dec if `ra_input.value` and `dec_input.value` are both already
  set.
- **APT couldn't load our session.json**: extra top-level `vmpt` key
  + missing `stats/errors/plannerSpecification` blocks. Fix: move
  vMPT-only data to a sidecar `vMPT_workspace.json`; fill all the
  nested null blocks with the reference shape.

---

## External references

- **eMPT** (export format inspiration): Bonaventura et al. 2023,
  A&A 672 A40 — [arXiv:2302.10957](https://arxiv.org/abs/2302.10957)
  / [GitHub](https://github.com/esdc-esac-esa-int/eMPT_v1).
- **JDox MPT Catalogs**: <https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-catalogs>
- **JDox JWST PA reference**: <https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-position-angles-ranges-and-offsets>
- **STScI APT downloader**: `https://www.stsci.edu/jwst-program-info/download/jwst/apt/<program_id>/`
- **MSA operability** (CRDS): `jwst_nirspec_msaoper_*.json`
- **Source of coordinate code**: `/Users/sunfengwu/jwst_cycle4/footprint_emerald.ipynb`
  (cells around In[54], In[151], In[164], In[166], In[125]).
