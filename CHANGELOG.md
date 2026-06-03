# Changelog

All notable changes to vMPT are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
project adheres to [Semantic Versioning](https://semver.org/).

## [1.1.1] — 2026-06-03

Patch release. Polish + several real bugs in the v1.1.0 optimizer
and catalog editor. The big change is that **Hierarchy mode actually
optimises lower tiers now**, plus a much richer results table.

### Optimizer

- **Slitlet centre is now right under the target.** The optimizer's
  `axy_to_shutter` returns 0-based fractional indices, but
  `_add_slitlet` expects 1-based — the missing `+1` was opening every
  slitlet one row up and one column to the left of the target. Now
  centred correctly.
- **Confirm dialog before Apply.** Clicking **Apply #N** opens a
  browser confirm: "This will CLEAR all previously open shutters and
  replace them with the optimizer's slitlets." OK → clears + applies
  (single Undo step); Cancel → no-op. Wired via `Button.js_on_click`
  → hidden trigger `TextInput` → Python handler. The trigger pattern
  is needed because `CustomJS.args` only accepts Bokeh Model
  instances, not floats — that's why the Apply button silently did
  nothing in v1.1.0; embed per-button scalars via Python f-string
  interpolation into the JS body.
- **Hierarchy mode now genuinely optimises every priority tier.**
  Previously DE refinement used `weights = 1 at top tier, 0
  elsewhere`, so DE happily slid to any pointing that kept the
  top-tier count even if it lost lower-tier sources in the process.
  DE now uses **auto-derived lex weights** (smallest int weights
  such that any higher tier strictly outweighs the sum of all lower
  tiers); their sum is a lex-equivalent scalar that DE maximises
  without violating priority ordering. The grid + multi-stage filter
  phase is unchanged.
- **Results table shows tier breakdown.** For Hierarchy, the Score
  column reads e.g. `P0:4 · P1:12 · P2:30 (46)` — per-tier source
  count + total in parens. For Meritocracy, `Σw 287.0 (46)`. For
  Democracy, just the count `46`.
- **Hover any Score cell** to see the top 10 placed sources at that
  pointing, sorted by priority ascending then weight descending —
  IDs + P + W per line.
- **Modal widened** to 740 px to fit the new columns; Score column
  width is method-specific; cells now `overflow: hidden +
  white-space: nowrap + text-overflow: ellipsis` so a label that
  overruns its column truncates instead of wrapping under the row.

### Catalog editor

- **Numeric sort on Priority + Weight.** Both columns are now stored
  as floats with NaN for missing (was strings → `"10"` < `"2"`
  lexicographically). An HTMLTemplateFormatter renders the cell as
  a rounded integer or blank; cell edits via StringEditor are
  coerced back to float in `_on_cat_edit_data_change` so the column
  stays a sortable numeric.
- **After a header click, the table scrolls to row 1.** Document-
  level click delegate on `.slick-header-column` resets the table's
  `.slick-viewport.scrollTop` to 0 (with an 80 ms delay so the
  re-render finishes first).
- **CSV save** uses `_fmt_int_or_blank` for Priority + Weight so
  the output is `5` not `5.0` and blanks stay blank.
- **`Compute w from p` / `Compute p from w`** write floats to the
  source (was strings) so the new column stays numerically sortable.

### Misc

- The `compute_weights_from_priorities` helper now correctly
  satisfies BOTH `w(p) > w(p+1)` AND `N(p)·w(p) > N(p+1)·w(p+1)`
  using `max(w_prev + 1, n_prev * w_prev // n_q + 1)` as the
  smallest integer that dominates the prior class (regression-tested
  in `tests/test_catalog_ops.py`).
- Loader: empty cells in numeric columns now properly become NaN
  even when the source column was masked-int (previously came
  through as 0).
- A couple of additional patterns added to `.gitignore` so stray
  personal files in the repo root can't accidentally be staged.

### Tests

- 124 passing, 1 skipped (same as v1.1.0; no test regressions).

[1.1.1]: https://github.com/fengwusun/vMPT/releases/tag/v1.1.1

## [1.1.0] — 2026-06-02

Headline feature: a complete **MSA pointing optimizer** with three
methods (Democracy / Meritocracy / Hierarchy), plus an editable,
sortable catalog editor. Several quality-of-life improvements
elsewhere.

### MSA pointing optimizer

- New panel at the bottom of the **Pointing** tab. Searches over
  (ΔRA, ΔDec, ΔPA) for an (RA, Dec, V3 PA) that maximises a
  user-selectable objective. Re-implemented in vMPT style from
  [hMPT](https://github.com/zihaowu-astro/hMPT) (Eisenstein, McCarty,
  Wu; CfA / Harvard); see `app/optimizer.py` for attribution +
  algorithm notes.
- **Three methods**:
  - **Democracy** — raw source count; ignores priority and weight.
  - **Meritocracy** — sum of `weight` of placed sources (MPT-style).
    Requires a populated weight column.
  - **Hierarchy** — strict priority-tier lex ordering (eMPT-style).
    Multi-stage filter: a higher-priority source is never traded for
    any number of lower-priority sources.
- **Pop-up modal** with an animated striped progress bar, a spinning
  ring, and a status line showing the current phase
  (`Grid: 5,200 / 20,000 · 4.2s elapsed · ~12s left`,
  `Hierarchy filter: tier 2 / 4 (p=1) — survivors: 18`,
  `Refining top 10: 3 / 10 · 7.4s elapsed`).
- **Results table** with the top-10 distinct solutions (near-
  duplicates collapsed). Each row pairs score + (ΔRA, ΔDec, ΔPA)
  with an **Apply #N** button.
- **Apply #N** sets the pointing AND opens an N-shutter slitlet
  (N from the Setting tab) at every observable target's shutter,
  auto-tagged with the catalog source ID. One Undo step reverts
  the whole apply.
- **ΔX = 0 freezes the axis** — set ΔPA = 0 to search RA/Dec only
  at the current roll, etc. Both the grid sweep and the DE
  refinement honour the freeze.
- **Advanced settings…** modal exposes grid resolution (n_RA, n_Dec,
  n_PA), DE max iterations, objective (count/flux), source σ, and
  the APT DVA θ.

### Catalog editor

- New **Edit catalog…** button in the Input tab opens a sortable,
  in-cell-editable spreadsheet pop-up.
  - Single-click any cell to edit. Tab / Enter commits; Esc cancels.
  - Drag inside a cell to highlight text; Cmd/Ctrl-C / Cmd/Ctrl-V
    copy / paste — a custom capture-phase keydown handler bypasses
    SlickGrid's column-copy default so only the selected text is
    copied.
  - 🗑️ icon at the end of each row deletes that row.
  - ↶ Undo / ↷ Redo for every edit, delete, derivation, and
    column add (100-step history).
- **Column picker** — toggle which columns are visible. Extras
  columns from the source CSV/FITS (the loader now preserves every
  column it didn't claim) live alongside the standard set and can
  be turned on or off.
- **Add a custom column** via a text input + button. Empty by
  default; useful for `weight`, `reference`, etc. Round-trips
  through Apply changes and Save as CSV.
- **Compute w from p** and **Compute p from w** buttons derive
  one column from the other:
  - `w(lowest p) = 1`; for each higher-priority class, the smallest
    integer `w(p)` satisfying `w(p) > w(p+1)` AND
    `N(p) * w(p) > N(p+1) * w(p+1)`. Guarantees strict-dominance:
    one source at any tier outweighs every source at all lower
    tiers combined.
  - `p` from `w` groups unique weights descending and assigns
    priorities 1, 2, 3, …
- **Save as CSV** with a Browse… file picker.
- **Apply changes & close** commits the working copy to the
  in-memory catalog so the eMPT bundle export reflects edits.

### Catalog model

- `Catalog.weight` is now a first-class field (sibling of
  `priority`). Loader detects `weight` / `w` / `wt` / `weights`
  aliases. Empty cells in numeric columns properly become NaN
  (previously masked integer columns silently became 0).
- Extras columns the loader didn't claim are preserved on the
  `Catalog.extras` dict (object arrays, original column name as
  key) and surfaced through the editor's column picker.

### UI polish

- **`run.sh`** gained `--port N`, `--fits PATH`, `--jpg PATH`,
  `--wcs PATH`, `--catalog PATH` (repeatable) flags. Mutual-
  exclusion rules: `--jpg` and `--wcs` come as a pair; `--fits` is
  exclusive with them.
- **Tabs renamed**: Image → Input, Aim → Pointing, Pick → Setting.
- **Pointing tab** now also hosts Disperser/Filter (was on the
  former Pick tab). RA/Dec inputs share a row; V3 PA/APA share a
  row; Visibility date + button share a row.
- **Canvas pixel aspect locked** — `frame_width` / `frame_height`
  match the loaded image's pixel W:H exactly. Window resizes
  letterbox around the canvas; the image is never stretched.
- **Sequenced autoload**: `run.sh --jpg ... --wcs ... --catalog ...`
  loads the image first and the catalogs strictly after, via an
  `on_complete` callback chain so the catalog overlay never races
  the image's `_set_image_and_recenter`.
- **Status bar** moved out of the scrollable sidebar column and
  pinned to the bottom-left of the viewport (position:fixed) so
  it can't render on top of tab content.
- **Optimizer Advanced settings** moved into a pop-up modal.
- **6 new tips** in the help-panel carousel — `run.sh` args,
  optimizer, catalog editor, multi-catalog, pixel aspect, big-ID
  mod.

### Tests

- 124 passing (was 96 at v1.0.1). New coverage: catalog `weight`
  column, mod-1e7 + empty-cell NaN handling, optimizer correctness
  (radec→Axy, quadrant inverse, centration monotonicity,
  Hierarchy-vs-Democracy divergence, dΔ=0 freezes, dedup), the two
  weight↔priority compute helpers in `app/catalog_ops.py`.

[1.1.0]: https://github.com/fengwusun/vMPT/releases/tag/v1.1.0

## [1.0.1] — 2026-05-21

Patch release. Two large quality-of-life corrections — accurate
per-shutter wavelength values, and a much friendlier catalog
loader — plus polish on the catalog UI and overlay defaults.

### Wavelength accuracy

- **Per-shutter dispersion table for every (disperser, filter)
  combo**, derived from numerical integration of the pipeline
  reference files via [`spacetelescope/msaviz`](https://github.com/spacetelescope/msaviz).
  Lives at `data/dispersion_cutoffs.npz` (19 MB compressed) and is
  regenerated by `scripts/precompute_dispersion_cutoffs.py`. Replaces
  the old linear V2-shift approximation that was wrong in two ways
  for PRISM:
  - Previous PRISM gap was held at 2.7–3.2 μm everywhere. Real
    gap location varies dramatically across the MSA (5–95 %
    spread: gap_lo 0.65–3.59 μm, gap_hi 3.03–5.02 μm) because
    PRISM dispersion is highly non-linear. The new lookup gives
    msaviz-accurate values per shutter.
  - PRISM endpoints used to drift with V2; in reality they're
    essentially constant (msaviz spread is ~0.01 μm).
- **Q3 / Q4 PRISM shutters** correctly report "no gap on this
  spectrum" — their spectra fall entirely on one detector.
- **Grating endpoints** updated to match the pipeline-reference
  `sci_range` instead of the slightly narrower JDox "useful range":
  - G140M/F100LP: 0.97–**1.89** (was 0.97–1.84)
  - G140H/F100LP: 0.97–**1.89** (was 0.97–1.84)
  - G235M/F170LP: 1.66–**3.17** (was 1.66–3.07)
  - G235H/F170LP: 1.66–**3.17** (was 1.66–3.07)
  - G395M/F290LP: 2.87–**5.27** (was 2.87–5.14)
  - G395H/F290LP: 2.87–**5.27** (was 2.87–5.14)
  - G140H/F070LP: **0.70**–1.27 (was 0.81–1.27)
- **vMPT does NOT depend on msaviz at runtime** — only the
  precompute script does. The shipped npz is everything the app
  needs.

### Catalog loader: looser column matching

- **`_norm()`** lowercases, strips bracketed/parenthesised unit
  annotations (`[deg]`, `(deg)`), collapses non-alphanumerics, and
  peels trailing unit/epoch tokens (`deg`, `degrees`, `rad`,
  `arcsec`, `J2000`, `ICRS`, `FK5`). All of these now match RA:
  `RA`, `ra`, `RA[deg]`, `RA(deg)`, `RA_deg`, `RAJ2000`,
  `Right Ascension`, `ALPHA_J2000`, `R.A.[deg]`. Same for Dec
  (including Vizier's `DEJ2000`).
- **ID resolution** accepts the usual aliases (`id`, `no`,
  `source_id`, `objid`, `srcid`, …) plus permissive fallbacks
  (`name`, `label`, `tag`, `target`, `#`) — fallbacks honoured only
  when values coerce to integer.
- **Missing ID column → synthesised** sequential IDs 1..N so the
  catalog still loads.
- **Numeric IDs ≥ 10⁷ are taken mod 10⁷** (`ID_MOD = 10_000_000`)
  so JADES-style 8–9-digit IDs collapse to APT's compact space.
- **Priority class strings** (`P0`, `P1`, …) and masked numeric
  cells now flow through cleanly — the old loader threw
  `ValueError` on a `P0` priority cell and produced `0.0` for
  masked `mag` / `z` instead of `NaN`.

### Multi-catalog

- **Load multiple catalogs at once**. Each gets a colour chip in
  the sidebar list. Toggle visibility per-catalog with a checkbox;
  × to remove; ▲ / ▼ to reorder the visual stack.
- **Per-catalog marker colours** cycle through an 8-entry palette
  (yellow / magenta / pale green / coral / lavender / sky-blue /
  white / salmon), picked to read clearly on dark fields and avoid
  the other overlay colours.
- **Z-order by list order** (earlier-loaded catalogs draw on top)
  with **alpha decay** by depth (1.0 → 0.35 floor). Matched-shutter
  targets always render fully opaque so a "picked" marker is never
  visually demoted.
- Sessions serialise the list (`catalog_paths`) — workspace JSON
  remembers each path and its enabled flag.

### MPT-importable catalog (export)

- **Output is now a superset of the input catalog** — every input
  source is included, plus any synthesised entries for slitlets
  without a real match. The `Label` column carries `real` or
  `vMPT_synth` so downstream tools can tell which is which.
- **Integer IDs only**, extracted as the largest digit run from
  the original token (so `RJ0600-10274-P0` → 10274). The original
  string token is preserved in the `Label` column for traceability.

### Overlay defaults

- Operable-shutter stroke: 0.75 px → **1.0 px**.
- Spectral-overlap fill alpha: 0.10 → **0.20**.
- Spectral-overlap edge colour now explicitly **orange (#d97a00)**
  — when you reveal the edge via the stroke slider it now matches
  the orange fill instead of Bokeh's default blue-grey.

### Tests

- 96 passing (was 63 at v1.0.0). Coverage growth concentrated on
  the catalog loader (column aliases, ID synth, mod-10⁷, string
  IDs, name-as-numeric-ID) and the wavelength model (every
  disperser × filter combo verified against the msaviz table).

[1.0.1]: https://github.com/fengwusun/vMPT/releases/tag/v1.0.1

## [1.0.0] — 2026-05-20

First public release. The tool is feature-complete for hand-picking
JWST/NIRSpec MSA shutter configurations on a target field and
exporting a bundle that loads into APT MPT and the eMPT pipeline.

### Highlights

- **Interactive shutter picker** with N-shutter slitlets (N ∈ {1, 2, 3, 5}),
  snap-to-nearest-operable, undo / clear, double-click highlights,
  shift-click to move the pointing, wheel-zoom and pan.
- **Live overlays** — MSA outline, operable shutters (silver edge),
  stuck-open (dark-red outline), user picks (red fill), spectral
  conflicts (orange fill, stackable), 5 fixed slits (gold), catalog
  targets (yellow / green when matched), lime pointing cross.
- **APT-ready bundle export** — 6 files per export, with role-prefixed
  filenames (`MPT_*`, `vMPT_*`, `eMPT_*`). The MPT plan JSON matches
  APT's reference schema field-for-field; the `<catalog>.cat` uses
  JDox-recognized column names (ID, RA, DEC, Weight, Primary, Label).
  Labels distinguish `real` catalog rows from `vMPT_synth` synthesised
  entries.
- **APT plan importer** — load any `MPT_plan.json`, shutter mask CSV,
  local `.aptx` archive, or fetch by JWST program ID directly from
  STScI. Reads multi-plan archives (e.g. program 1208 with 40+ plans).
- **Bundle round-trip** — Save session → load session restores
  pointing, V3 PA, disperser/filter, every open shutter with its
  `target_id` + `role`, the highlighted set, and the image + sidecar
  paths. Point at either `MPT_plan.json` OR `vMPT_workspace.json` —
  the sibling auto-loads.
- **Responsive layout** — canvas stretches to fill the browser
  window; sidebar / help panel scroll on overflow; left-sidebar
  fixed at 340 px, right help panel at 340 px.
- **Rotating tip card** in the help panel (13 hand-written tips,
  15-second rotation with CSS fade-in).
- **GitHub version-check on startup** — non-blocking background
  thread compares the local HEAD to `origin/main`; shows a
  dismissible amber notification if the local copy is behind.
- **Custom favicon** (4 MSA quadrants + lime pointing cross).
- **One-page summary slide** generator (`build_vmpt_slide.js`,
  pptxgenjs-based).

### Science correctness

- **MSA geometry** sourced from `pysiaf` (`NRS_FULL_MSA`); 138.575°
  intra-MSA rotation, V2/V3 reference at (378.563, −428.403).
- **APA = V3 PA + V3IdlYAngle (mod 360)**; both quantities are
  surfaced in the status bar and editable from the Aim tab.
- **Operability** read from CRDS `jwst_nirspec_msaoper_*.json` —
  failed-open shutters always disperse and contribute to the
  spec-overlap calculation.
- **Spectral overlap** — `|Δs| ≤ 1` cross-quadrant via NRS1 (Q1↔Q3)
  and NRS2 (Q2↔Q4) detector pairing; per-grating V2 half-extent
  (PRISM 35″, M-gratings 200″, H-gratings 500″).
- **Wavelength endpoints** per disperser+filter, clamped to the
  grating's intrinsic range (no spurious PRISM > 5.3 µm tooltips).
- **Source matching** uses APT's *Unconstrained* Source Centering
  rule (full shutter pitch including bars).
- **WCS Jacobian** uses `astropy.SkyCoord.spherical_offsets_to` —
  cos(Dec) factor handled correctly at non-equatorial fields.

### Example data shipped

- `example_a370/` (43 MB) — JWST NIRCam F182M+F200W+F210M FITS of
  Abell 370, target catalog, GTO-1208 APT MPT plan, shutter-mask CSV.
- `example_r0600/` (21 MB) — JWST NIRCam F090W+F200W+F444W JPG of
  RXCJ0600 + WCS sidecar + 28k-source target catalog. JPG re-encoded
  at quality 85 (was 251 MB) without changing WCS.

### Tests

- 63 tests, ~5 s. Run with `pytest tests/`.
- Coverage: session bundle round-trip, MPT plan parser (incl. .aptx
  archives), eMPT format byte-compatibility, MPT catalog writer
  format guard, wavelength model, image loaders, end-to-end export.

### Known limitations

- `plannerSpecification` block in `MPT_plan.json` carries sensible
  defaults (matching APT's reference schema) but its dither /
  search-grid parameters don't reflect any vMPT internal state —
  APT uses them only as starting values for re-planning.
- Bokeh single-session state: opening the same server in two
  browser tabs lets picks bleed across them. Use one tab per user.
- Older `pysiaf` PRD (PRDOPSSOC-068) lags the online version by
  ~0.05″ for some apertures; safe to ignore unless you need
  milli-arcsec geometry.

### Acknowledgements

Export-bundle format calibrated against [eMPT](https://github.com/esdc-esac-esa-int/eMPT_v1)
(Bonaventura et al. 2023, A&A 672, A40). Coordinate plumbing builds
on `pysiaf` (NIRSpec apertures) and `astropy.wcs`. Visibility
windows queried via [`jwst_gtvt`](https://github.com/spacetelescope/jwst_gtvt).
MPT catalog and plan JSON schemas follow the
[JDox MPT documentation](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template).

[1.0.0]: https://github.com/fengwusun/vMPT/releases/tag/v1.0.0
