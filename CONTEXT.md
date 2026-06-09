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
├── wavelengths.py     Per-shutter dispersion model + lookup table
├── optimizer.py       MSA pointing optimizer (hMPT-derived, see docstring)
├── image_io.py        FITS + JPG-with-sidecar loaders (LoadedImage)
├── catalog.py         Catalog reader → Catalog dataclass
├── empt_io.py         eMPT-format writers + MPT-importable .cat writer
├── session_io.py      Bundle save/load (MPT plan + workspace sidecar)
├── mpt_io.py          APT MPT JSON parser + .aptx archive reader
├── static/favicon.svg vMPT favicon (MSA grid + pointing cross)
└── templates/index.html  Injects favicon via data URI

data/
├── nirspec_msa_v2v3.npz   (4 × 171 × 365) shutter V2/V3 arcsec
└── dispersion_cutoffs.npz per-shutter wavelength bounds (all 9 disperser
                           × filter combos; built from msaviz)

scripts/
└── precompute_dispersion_cutoffs.py  (re-)generates dispersion_cutoffs.npz

tests/                 pytest suite (110+ tests; ~10 s)
example_a370/          Abell 370 FITS (44 MB)
example_r0600/         RXCJ0600 JPG + sidecar (240 MB)
exports/               default output dir
```

### Sidebar tabs (current names)

- **Input** — image + catalog loading (FITS, JPG+sidecar, multi-catalog).
- **Pointing** — RA/Dec, V3 PA, disperser/filter, visibility window,
  pointing-optimizer panel.
- **Setting** — layers toggle, slitlet size, snap-to-operable, overlay
  appearance picker, undo/clear.
- **MPT** — import APT plans, session save/load, export eMPT bundle.

Older tab names (`Image / Aim / Pick`) were renamed in the May 2026
UI pass. Comments still using the old names should be updated when
touched.

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
| `catalog` | `Catalog` or `None`. Merged-active cache: union of every enabled entry in `catalogs`. Re-derived from scratch on add/remove/toggle so existing readers don't need to know there are multiple. |
| `catalogs` | List of `{name, catalog, enabled, color}` entries — the multi-catalog source of truth. Catalogs render on the canvas at their assigned palette colour; the per-entry checkbox toggles visibility, the × button removes the entry. |
| `catalog_colors` | `np.ndarray[object]` of marker colours parallel to `catalog.ra_deg`. Populated by `_rebuild_merged_catalog` so the overlay can emit per-source `line_color` strings without looking up `catalogs` for every row. |
| `ra_deg`, `dec_deg` | Pointing center. |
| `pa_v3` | V3 PA in degrees (mod 360). APA = `pa_v3 + V3_IDL_Y_ANGLE`. |
| `disperser`, `filter` | e.g. `"PRISM"`, `"CLEAR"`. |
| `slitlet_height` | N ∈ {1, 2, 3, 5}. Determines `_slitlet_offsets` and toggle-off siblings. |
| `open_shutters` | `dict[(q, s, d) → OpenShutter]`. The user's picks. |
| `highlighted` | `set[(q, s, d)]` — cyan-edge visual flag, not exported. |
| `history` | Undo stack (capped at 50 snapshots of `open_shutters`). |
| `shutter_to_catids` | `dict[(q, s, d) → [source_id, …]]` — catalog sources whose footprint lands in each shutter. Rebuilt on pointing / PA / catalog change. |
| `snap_to_operable` | When a target click misses, snap to the nearest operable shutter. |
| `catalog_alphas` | Per-source line_alpha (decayed by z-depth in the catalog stack so earlier-loaded catalogs read more strongly than later ones). |
| `_autoload_active` | Guard flag set while `_autoload_from_args` is driving sequenced loads from `run.sh` args. The path-input on_change handlers no-op while it's True so they don't double-trigger loads. |

Two module-local dicts hold transient work:

- `_opt_run` (in `main.py`) — in-flight optimizer state machine for
  the chunked grid + DE driver. Cleared on Close / when the run
  finishes / on a fresh run.
- `_inverse_cache` (in `app/optimizer.py`) — per-quadrant
  CloughTocher2D interpolators (Axy → fractional shutter indices).
  Lazy-built on first use (~2 s for the Delaunay triangulations).

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
10. **Catalog targets** — coloured circles (if layer toggled on).
    Unpicked sources take the loaded catalog's palette colour
    (`CATALOG_COLOR_PALETTE` in `main.py` — yellow / magenta / pale
    green / coral / lavender / sky-blue / white / salmon, cycled by
    load order). Picked sources (target_id matches an open shutter)
    flip to green (`#2e9b3f`) with a thicker line so matched vs.
    unmatched is obvious at a glance.

### Spec-overlap calculation (the subtle one)

For each dispersion source (every user-open + every stuck-open
shutter), the set of overlap-affected operable shutters is:

```
(s ∈ open ± SHVAL_S_TOLERANCE)                       AND
(primary_detector(q,s,d) == primary_detector(open))  AND
(|ΔV2| < v2_overlap_distance(disperser, filter))
```

Current values:
- `SHVAL_S_TOLERANCE = 1` — only the open shutter's row + immediate
  neighbours above/below disperse onto overlapping detector pixels
  (matches eMPT's `shval ≈ s` exactly).
- `primary_detector(disperser, filt, q, s, d)` is a per-shutter
  lookup into the new per-detector wavelength arrays added to
  `dispersion_cutoffs.npz` (`*_nrs1_lo/_hi`, `*_nrs2_lo/_hi`).
  Returns the detector that carries the larger λ-range for that
  shutter (0=NRS1, 1=NRS2, -1=off-detector). Replaces the v1.0–
  v1.3.x static `NRS1_QUADS={1,3}` / `NRS2_QUADS={2,4}` pairing,
  which was approximately right for PRISM but wrong for grating
  modes where (e.g.) Q4 G395M actually lands on NRS1, not NRS2.
- `v2_overlap_distance(disperser, filter)` is the *full* on-detector
  V2 extent of the spectrum, measured from `slit_frame → detector`
  traces in stenv: detector x-span × ≈ 0.077 ″/V2-px. Two same-row,
  same-primary spectra of length L share a detector pixel iff
  `|ΔV2| < L`.

  | Disperser | Filter | Full extent |
  |---|---|---|
  | PRISM | CLEAR | 32″ |
  | G140M | F070LP | 98″ |
  | G140M | F100LP | 109″ |
  | G235M | F170LP | 110″ |
  | G395M | F290LP | 103″ |
  | G140H | F070LP | 185″ |
  | G140H | F100LP | 307″ |
  | G235H | F170LP | 300″ |
  | G395H | F290LP | 281″ |

  Earlier values were sometimes treated as "half-extent" and
  sometimes as full extent; the table is now unambiguously the full
  V2 extent. The user-reported Q4 d=349 G395M case is correctly
  suppressed by the new `primary_detector` check (Q4 → NRS1, Q2 →
  NRS2) even though their V2 separation is well inside the 103″
  extent — the distance check is just there to catch same-detector
  same-quadrant collisions (Q4 d=349 ↔ Q4 d=50 same s, ΔV2 ≈ 60″,
  real overlap ≈ 557 px on NRS1).

### Three-colour MPT-faithful overlap classification

`refresh_overlays` records two count dicts per source type — `direct`
(candidate's row is in the slitlet's actual row range, with tilt) and
`buffer` (candidate is at the ±1 tolerance edge) — across operable AND
stuck-open candidates. After the contamination pass it derives:

- `hit_sources[i]`: set of `(source_type, source_idx)` tuples that
  produced a hit on candidate `i`. Used for chain-propagation.
- `conflicted_user` / `conflicted_stuck`: sets of slitlet indices
  whose OWN shutters are hit by another source. These are slitlets
  in active touching collision (no operable row between them and
  another open).

Three CDS (`src_spec_overlap_stuck/_user/_both`) feed three
`multi_polygons` glyphs. Per-polygon `fill_alpha` field; alpha is
`min(1, base × n_total)` so contamination from multiple sources
intensifies.

Classification rule applied to every contaminated candidate:

| Case | Colour |
|---|---|
| user-pick with any hit (direct or buffer) from another source | **purple** (Mask Conflict) |
| operable hit by a CONFLICTED source (chain propagation) | **purple** |
| operable hit by user-source(s), none conflicted | **orange** (Masked) |
| operable hit by stuck-source(s) only, none conflicted | **pink** (Mask Stuck) |
| no hit / hit only by an excluded path | no overlay |

The chain rule means: opening a slitlet in a clean (silver-edged)
area only ever produces orange/pink warnings, never purple — even
across many shared dispersion targets. Purple appears only when at
least one of the contributing slitlets is itself touching another
open shutter (so its row ranges are adjacent, no operable row
between them).

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

## Catalog loader (`app/catalog.py`)

The loader is intentionally permissive — catalogs come in from every
corner of the community and we'd rather tolerate weird header
spellings than force users to hand-edit.

### Loose column matching

`_find_col` normalises both the table's column names AND each
candidate before comparing, via `_norm()`:

1. Lowercase.
2. Strip bracketed / parenthesised unit annotations: `[deg]`,
   `(deg)`, `[arcsec]`.
3. Collapse remaining non-alphanumerics: `RA_deg` → `radeg`,
   `R.A.` → `ra`.
4. Peel off trailing unit / epoch tokens (`deg`, `degrees`, `rad`,
   `radian`, `arcsec`, `arcseconds`, `j2000`, `icrs`, `fk5`) in a
   loop, so `radeg` → `ra`, `decJ2000` → `dec`.

This means `RA`, `ra`, `RA[deg]`, `RA(deg)`, `RA_deg`, `RAJ2000`,
`Right Ascension`, `ALPHA_J2000`, `R.A.[deg]` all map to the same
key. The full candidate lists are `_RA_KEYS`, `_DEC_KEYS`,
`_ID_KEYS`, `_PRI_KEYS`, `_MAG_KEYS`, `_Z_KEYS`, `_LABEL_KEYS`.

### ID resolution + mod 10⁷

`_find_id_col` tries `_ID_KEYS` strictly, then falls back to
`_ID_FALLBACK_KEYS` (`name`, `label`, `tag`, `target`, `#`) — but the
fallback path is only accepted when the column's values coerce to
`int`, so a `name` column full of `"NGC-123"` does NOT silently
become the ID column.

If no usable ID column is found, vMPT **synthesises sequential IDs
1..N** so the catalog still loads (the user can refine later).

Numeric IDs are passed through `_coerce_int_ids`, which takes any
value `|id| >= ID_MOD = 10_000_000` mod ID_MOD. This keeps the
integer space compact for APT MPT and eMPT (which both expect
short numeric source numbers); JADES-style 8–9-digit IDs collapse
to 7 digits cleanly. Collision risk after the mod is treated as
acceptable.

String IDs (`"RJ0600-12345678-P0"`) are kept verbatim through the
in-memory model; the integer extraction happens in `main.py`'s
`_to_int_id` at export time, which picks the largest digit run from
the token.

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

`app/wavelengths.py`. Per-shutter lookup table for every (disperser,
filter) combo, derived from spacetelescope/msaviz's numerical
integration of the pipeline dispersion polynomials.

The table lives at `data/dispersion_cutoffs.npz` and is regenerated
by `scripts/precompute_dispersion_cutoffs.py`. For each combo it
stores four (4, 171, 365) float32 arrays under keys
`{DISPERSER}_{FILTER}_blue_edge`, `..._gap_lo`, `..._gap_hi`,
`..._red_edge`. NaN = the shutter's spectrum doesn't reach that
detector (e.g. PRISM Q3/Q4 shutters never span the NRS1/NRS2 gap;
H-grating edge shutters can miss one detector entirely).

`cutoffs(v2, v3, disperser, filt, *, q, s, d)` does O(1) lookup
when shutter indices are supplied; for the existing call sites that
don't have indices yet, a linear V2-shift fallback returns the
fiducial endpoints (used only by the old tests and a safety net for
fresh checkouts that haven't run the precompute).

`v2_overlap_distance(disperser, filter)` returns the V2 half-extent
used by the spec-overlap calculation. Looked up per (disperser, filter)
combo — see the table in the dispersion-row tilt section above
(18″ PRISM, 50–55″ M-gratings, 95–155″ H-gratings — depends on filter).
All values measured from `slit_frame → detector` traces in stenv.

### Regenerating the table

```
git clone https://github.com/spacetelescope/msaviz.git /tmp/msaviz
PYTHONPATH=/tmp/msaviz python scripts/precompute_dispersion_cutoffs.py
```

Takes ~13 minutes (PRISM is the slow one — per-shutter ODE
integration; gratings are batch-integrated per quadrant). vMPT does
NOT depend on msaviz at runtime — only the precompute does.

---

## MSA pointing optimizer (`app/optimizer.py`)

Searches for an (RA, Dec, V3 PA) that maximises the count — or
weighted flux — of catalog sources placed in operable, well-centred
MSA shutters. Re-implemented in vMPT style from
[**hMPT**](https://github.com/zihaowu-astro/hMPT) (Eisenstein,
McCarty, Wu; CfA/Harvard), itself derived from ESA's eMPT
(Bonaventura et al. 2023). Credit in the module docstring.

### Algorithm

1. **`radec_to_axy`** — vectorised gnomonic projection of source
   (RA, Dec) → MSA aperture plane (ax, ay) in arcsec. Includes the
   APT DVA correction (parameter `theta_deg`, default 90 = no shift)
   and the PA + intra-MSA rotation.
2. **`axy_to_shutter`** — per-quadrant `scipy.interpolate.
   CloughTocher2DInterpolator` maps (ax, ay) → fractional shutter
   indices `(quad, s_frac, d_frac)`. Built lazily on first call
   (Delaunay triangulation over ~62 k points per quadrant, ~2-3 s
   one-time). Cached at module level via `_inverse_cache`.
3. **`PointingEvaluator.evaluate`** — combines the above with the
   CRDS operability mask (incl. configurable 3-shutter vertical
   slit), an APT-style centration buffer (`CENTRATION_BUFFERS` dict
   — UNCONSTRAINED → TIGHTLY_CONSTRAINED), and a Gaussian-PSF
   throughput integral through the 0.27″×0.53″ aperture.
4. **`grid_search`** — brute force over a (ΔRA, ΔDec, ΔPA) cube.
   `ΔX = 0` freezes that axis (`n_X` forced to 1; the central
   value is the only sample). Returns the top-ranked candidates.
5. **`refine_top`** — `scipy.optimize.differential_evolution`
   polish of the top-N grid candidates inside a small bound box.
   When ΔX = 0 the corresponding variable is dropped from the
   DE problem (scipy can't take zero-width bounds). After
   refinement the list is sorted then **de-duplicated** within
   tolerance (default 0.3″ in RA / Dec, 0.05° in PA) so the user
   doesn't see N copies of the same pointing when the score
   landscape has a plateau.

### Method (Democracy / Meritocracy / Hierarchy)

The **Method** dropdown determines how source counts are weighted:

- **Democracy** — uniform weight (1) per source. Maximises count.
  Works with any catalog.
- **Meritocracy** — weights = `Catalog.weight` (NaN → 0). Maximises
  the sum of weights of placed sources. Requires a populated weight
  column.
- **Hierarchy** — strict priority-tier lex ordering. Requires a
  populated priority column. Multi-stage:
  1. Grid search with uniform weights (Democracy-style score)
     gives a candidate pool of up to K = 100 pointings.
  2. For each priority tier (ascending priority value = descending
     importance), re-evaluate every surviving candidate with
     `weights = 𝟙(priority == tier)` and keep only those tied with
     the per-tier max. Single tick per tier; fast.
  3. DE-refine the surviving pool with top-tier weights.

  Same result as constructing a single "huge multiplier per tier"
  weights vector, but explicit multi-stage filtering makes the
  intermediate counts inspectable in the progress modal.

Sources with NaN in the required column contribute weight 0 — they're
visible on the canvas but invisible to that mode's optimizer.

### Collision protection (v1.2.0+)

`PointingEvaluator` accepts optional `protect_mask`, `priorities`,
`weights`, `disperser`, `filt`, `reason` kwargs. When `protect_mask`
flags any source, three drop rules apply at every `evaluate()` call,
reusing the **same** physics as the live canvas's spec-overlap layer
(Q1/Q3 → NRS1 vs Q2/Q4 → NRS2 detector halves, V2 separation
< `v2_overlap_distance(disperser, filt)`). Row tolerance is
**slitlet-aware** (v1.2.1+) — `SHVAL_S_TOLERANCE = 1` is the
per-individual-shutter constant used by the live-canvas glyph, but
the optimizer's evaluator pre-computes two slitlet-aware tolerances
in `_init_protection`:

  - `_sd_tol_ps = slit_length // 2 + 1` — protected slitlet against
    a single stuck-open shutter (rule 1).
  - `_sd_tol_pp = 2 * (slit_length // 2) + 1` — protected slitlet
    against another slitlet, both using the optimizer's global
    `slit_length` (rules 2 + 3).

For the default `slit_length=3` (`half=1`) these are 2 and 3 — i.e.
no other shutter at rows `s_p±2`, no other slitlet centered closer
than `s_q − s_p > 3`. v1.2.0 used a flat `|Δs|≤1` between centres
which under-counted collisions for any `slit_length > 1`.

Rules:

1. **Protected ↔ stuck-open** — a protected source on a row colliding
   with any REASON==2 shutter is dropped (unavoidably contaminated).
2. **Protected ↔ protected** — within each colliding cluster, the
   lowest-priority-number source wins. Ties on priority break on
   higher weight; ties on weight break on lower source index
   (stable). Encoded as a per-source `_collision_rank` precomputed in
   `_init_protection`.
3. **Protected (still kept) ↔ unprotected** — every unprotected
   source colliding with any still-kept protected one is dropped.

Dropped protected sources do **not** propagate collision pressure to
rules 2/3 — losing one high-priority spectrum doesn't justify
compounding the loss by also dropping unprotected sources from the
same row. `evaluate_with_stats(ra, dec, pa)` returns the
3-tuple-plus-drop-count; `evaluate(...)` still returns the 3-tuple
but its `detected` is now the kept (post-drop) mask, so existing
scoring code (`np.sum(det * weights)` etc.) automatically reflects
the protection without changes.

The UI builds the `protect_mask` in `_protect_mask_for_catalog`:
priority cutoff `pri ≤ X` OR weight cutoff `wgt ≥ Y` (mutually
exclusive — `opt_protect_mode_radio.active` is 0 or 1).
`opt_protect_status_div` updates live via `_update_protect_status_div`
hooked to the checkbox / radio / threshold input. The results modal
appends `−K` to each Score cell and prefixes 🛡 on protected IDs in
the hover top-10; `n_dropped[i]` is computed in `_opt_de_step` by
calling `evaluate_with_stats` for each refined candidate.

Important: `_rebuild_merged_catalog` had to be extended to propagate
`weight` for the multi-catalog case (was only carrying it through
the single-catalog fast path); without this fix the "By weight ≥"
rule would silently select zero sources in multi-catalog mode.

### Catalog editor (`app/catalog_ops.py`)

The editor's **Compute w from p** and **Compute p from w** buttons
call `compute_weights_from_priorities` / `compute_priorities_from_weights`
from a separate module (`app/catalog_ops.py`) so the formulas can be
unit-tested without importing Bokeh.

The w-from-p formula iterates from the LOWEST priority class
(largest p) upward:
- `w(lowest_p) = 1`
- For each higher class p, find the smallest integer `w(p)` such that
  `w(p) > w(p+1)` and `N(p) * w(p) > N(p+1) * w(p+1)`.

This guarantees that one source at any priority tier strictly
outweighs every source at all lower tiers COMBINED — i.e. equivalent
to lex ordering. The multi-stage Hierarchy filter achieves the same
behaviour without needing huge integers.

### Driver in main.py (`_opt_drive`)

The UI runs the grid + DE asynchronously as a state machine so the
Bokeh IO loop can keep rendering. Phases:

- `phase = "grid"`: process `_OPT_GRID_CHUNK = 400` pointings per
  tick (~0.5 s), then re-schedule via `curdoc().add_next_tick_callback`.
  Progress fraction 0 → 0.85 of the bar.
- `phase = "hierarchy"`: ONLY for the Hierarchy method. Multi-stage
  tier filter, one tier per tick (cheap — ~ms per tier). Progress
  0.85 → ~0.95 (sub-step proportional to tier index).
- `phase = "de"`: one DE refinement per tick. Progress remaining → 1.0.
- `phase = "done"`: results rendered in the pop-up modal.

### Modal UI

Two-state pop-up (`opt_modal_card`, with `opt_modal_backdrop`):

1. **Progress** — spinning ring + status line + animated striped
   progress bar. The bar uses a CSS custom property
   (`--vmpt-pct`) on the wrapper Div so width updates don't
   replace the inner DOM — the stripe / glow animations keep
   running uninterrupted across the dozens of progress updates.
2. **Results** — top-10 distinct solutions as a Bokeh column of
   `row(rank, score, ΔRA, ΔDec, ΔPA, Apply_btn)`. Each Apply
   button sets RA/Dec/V3 PA via `_apply_optimizer_result` and
   closes the modal. **Picks are NOT auto-placed** — the user
   chooses how to fill shutters.

Advanced settings (grid resolution, DE max-iter, objective, σ, θ)
live in a separate pop-up (`opt_advanced_modal_card`). The widgets
retain their values regardless of modal visibility; the optimizer
reads them when Run is clicked.

**v1.2.1 layout** — the Pointing tab no longer holds the optimizer
inputs inline. A single primary `Open optimizer…` button
(`opt_open_btn`) opens `opt_config_modal_card`, which holds the
Method dropdown, ΔRA/ΔDec/ΔPA / N inputs, centration, priority
cutoff, the Protect-spectra group, the Advanced settings opener,
and the actual `Run optimization` button + status line. Both the
config and advanced modals stack above the results modal:

  - results card / backdrop  → `z-index` 1000 / 999
  - optimizer config card / backdrop → `z-index` 1000 / 999
  - advanced settings card / backdrop → `z-index` 1002 / 1001

The advanced bump is so opening `Advanced settings…` *from inside*
the config modal actually draws above it; v1.2.0 had a tie that
left the advanced card hidden behind the (newer-in-DOM) config
modal. Both pop-up cards (catalog editor and config) carry a
top-right `×` dismiss button positioned via inline
`position: absolute` styles (the Bokeh nested wrapper layout
defeats outer-page CSS for that one corner).

### Performance

For a typical run (~500 sources, 20³ = 8 000 grid pointings, top-10
DE refinement) total wall time is **~5–15 s** depending on machine.
Grid dominates; per-pointing eval is ~1 ms after the
CloughTocher2D triangulation is cached.

---

## CLI auto-load (`run.sh --args`)

`run.sh` forwards `--port`, `--fits`, `--jpg`, `--wcs`,
`--catalog` (repeatable) via Bokeh's `--args`. Inside main.py,
`_autoload_from_args` parses `sys.argv`, builds a **sequenced
queue** of (image → catalogs) loaders, and invokes them via an
explicit on_complete chain — each loader's `finally` block
schedules the next step on the next tick.

Sequencing is critical: the catalog overlay's pixel positions
come from the image WCS, so running the catalog load in the same
tick as the image load races against `_set_image_and_recenter`
and breaks the canvas aspect. The guard flag
`state["_autoload_active"]` muzzles the on_change handlers on
the path-input widgets while autoload drives them directly, so
each input value-set doesn't fire its own redundant load.

### 1:1 pixel aspect lock

The figure uses `frame_width` + `frame_height` (the inner canvas
frame's pixel dimensions) — NOT `sizing_mode` + `aspect_ratio`.
`refresh_image_glyph` sets `frame_width` / `frame_height` to
match the image's W:H exactly, with the longer axis pinned to
`FRAME_MAX = 800 px`. Window resizes change only the surrounding
black space; the canvas itself never reflows. Earlier attempts
with `stretch_both + match_aspect` and `scale_both + aspect_ratio`
both ended up distorting the image during reflow.

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
pytest tests/    # 110+ tests, ~10 s
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
- `tests/test_wavelengths.py` — per-shutter table values match
  msaviz fiducial; gap_lo varies > 0.5 µm across MSA for PRISM;
  Q3/Q4 PRISM shutters have no gap.
- `tests/test_optimizer.py` — radec → axy + quadrant lookup
  correctness, hMPT-published-example score range,
  centration-class monotonicity, ΔX=0 freezes axes, refine_top
  dedups, flux objective differs from count.
- `tests/test_catalog.py` — loose column matching (RA[deg],
  RAJ2000, DEJ2000, etc.), ID synth, mod-10⁷, name-as-numeric-ID.

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
- **PRISM gap wrong everywhere except the central shutter**: the
  old model held the NRS1/NRS2 gap at a fixed wavelength pair
  (≈ 2.7–3.2 μm). msaviz integration showed the gap actually varies
  by ±1.5 μm across the MSA because PRISM dispersion is non-linear.
  Fix: precomputed per-shutter table at `data/dispersion_cutoffs.npz`
  (built from msaviz; 9 disperser × filter combos).
- **Canvas aspect distorted on window resize** in run.sh autoload
  mode: `stretch_both` + `match_aspect=True` looked right standing
  still but stretched the image when the user dragged the window
  edges. Even `scale_both` + `aspect_ratio=W/H` wasn't reliable.
  Fix: use `frame_width`/`frame_height` (inner canvas frame pixels)
  pinned to the image's W:H exactly; no `sizing_mode` on the figure.
  Window resizes change only the surrounding letterbox.
- **Catalog overlay raced the image load** in run.sh autoload mode:
  setting fits_path + catalog_path in the same tick fired both
  loads in parallel, with the catalog's pixel-position computation
  using stale WCS. Fix: each loader now takes `on_complete=cb`;
  `_autoload_from_args` chains image → catalogs explicitly via
  these callbacks. The `state["_autoload_active"]` guard muzzles
  the path-input on_change handlers while autoload is driving.
- **Optimizer Apply buttons drifted out of alignment with HTML
  table rows**: separate "HTML table + Bokeh button column" pattern
  accumulated 5 px of column-spacing per row in the buttons
  column. Fix: build the table as `column(header_row, row(cell,
  cell, cell, cell, cell, Apply_btn), ...)` with `spacing=0` —
  one Bokeh row per result so the Apply button is a sibling of
  its own cells.
- **Progress-bar CSS animation restarted on every progress tick**:
  Bokeh's `Div.text = ...` replaces innerHTML, which unmounts the
  bar element and resets its `@keyframes`. Fix: bar HTML is set
  once at construction; only the wrapper's `--vmpt-pct` CSS
  custom property changes per tick (via `Div.styles`), which
  Bokeh applies without touching innerHTML. The stripe + glow
  animations now run continuously across the entire optimization.
- **Optimizer top-10 list filled with identical solutions** when
  ΔX = 0 was used: DE in a zero-width box returned the same
  optimum from every starting point. Fix: `refine_top` drops the
  frozen variable from the DE problem entirely, and dedups
  near-identical results before returning (default tolerance:
  0.3″ RA/Dec, 0.05° PA).

---

## External references

- **eMPT** (export format inspiration): Bonaventura et al. 2023,
  A&A 672 A40 — [arXiv:2302.10957](https://arxiv.org/abs/2302.10957)
  / [GitHub](https://github.com/esdc-esac-esa-int/eMPT_v1).
- **JDox MPT Catalogs**: <https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-catalogs>
- **JDox JWST PA reference**: <https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-position-angles-ranges-and-offsets>
- **STScI APT downloader**: `https://www.stsci.edu/jwst-program-info/download/jwst/apt/<program_id>/`
- **hMPT** (optimizer algorithm): Eisenstein, McCarty, Wu (CfA),
  [GitHub](https://github.com/zihaowu-astro/hMPT). vMPT
  re-implements the algorithm; see `app/optimizer.py` docstring.
- **msaviz** (per-shutter dispersion reference, used by precompute
  script only): <https://github.com/spacetelescope/msaviz>.
- **MSA operability** (CRDS): `jwst_nirspec_msaoper_*.json`
- **Source of coordinate code**: an internal cycle-4 footprint
  notebook (cells around In[54], In[151], In[164], In[166], In[125]).
