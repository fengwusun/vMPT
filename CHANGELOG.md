# Changelog

All notable changes to vMPT are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### MPT-style spec-overlap colors

vMPT now uses the three-colour encoding that APT MPT shows on its
shutter map. The single orange "spec-overlap" layer from v1.0–v1.3.0
is split into three independent layers, each with its own colour
and its own slider in **Settings → Overlay appearance**:

- **Mask Stuck (pink)** — a shutter whose spectrum overlaps with
  a **stuck-open** shutter's dispersion. Computed once per
  (disperser, PA) — visible even before the user picks anything.
- **Masked (orange)** — overlap with a **user-open** shutter's
  dispersion. Updates on every shutter pick.
- **Mask Conflict (purple)** — overlap with **both** at once. A
  shutter that was pink turns purple when a user pick adds a
  second overlap source.

Per-polygon alpha now stacks with the conflict count:
`alpha = min(1, base × n_conflicts)`. Base is 0.20 (= the alpha
of a shutter overlapping a single dispersion source). A shutter
hit by 2 sources renders at 0.40; 5+ saturates at 1.0. The base
per category is adjustable in Settings and persists via the
preferences system (`~/.vmpt/preferences.json`); migration from
the old `Overlapping shutters` prefs key drops the saved value
onto all three new categories.

Implementation: `src_spec_overlap_stuck` / `_user` / `_both`
ColumnDataSources with a `fill_alpha` field column; three
`multi_polygons` glyphs read it field-referenced
(`fill_alpha="fill_alpha"`). The overlap calculation now tracks
per-shutter source counts in `stuck_counts` / `user_counts` dicts
rather than the previous `np.unique(idx)` flatten — so the
"alpha compositing stacks" comment that lived in the v1.0 code
finally has its intent realised.

The top stats-bar "Conflict shutters" cell counts the union of all
three categories. The aggregate `overlap_idx` is still computed
for the silver-edge layer's "exclude affected shutters" filter.

The pre-existing `src_spec_overlap` and `spec_overlap_glyph`
module symbols are kept as backwards-compat aliases pointing at
the orange (user-overlap) source/glyph, so external scripts and
tests that imported them keep working.

## [1.3.0] — 2026-06-04

Big feature release. Adds per-target spectral constraints (a whole
new constraint system), per-target centration overrides, persistent
user preferences, draggable modal dialogs, a confirm-overwrite gate
on file saves, and a long list of UX polish + bug fixes.

### Per-target spectral constraints (headline feature)

The catalog editor gains a per-row **Constraints…** button (gray
when unset, primary-blue when ≥1 field is set). Clicking it opens
a popover that captures spectral requirements for that source. At
every candidate pointing the optimizer fetches the source's
centre-shutter wavelength endpoints via `vmpt.wavelengths.cutoffs`
and drops the source if any constraint fails — same machinery as
the v1.2.0 collision drop, just with new reason codes.

Five per-target fields:

`required_lam`
: List of `(λ_lo, λ_hi)` ranges in μm that must land on the
  detector for this source (gap excluded — a range straddling the
  NRS1/NRS2 gap fails). Empty list = no constraint. Editor input
  format: `"1.0-1.3; 1.5-1.8"`.

`no_gap`
: Boolean. When True, the NRS gap must not fall inside the
  source's spectrum. PRISM at non-central shutters has no gap →
  passes. H gratings have a gap inside every shutter's spectrum →
  fails.

`extend_blue` / `extend_red`
: Booleans. Shutter must reach the disperser/filter's MSA-wide
  best blue / red wavelength (within 20 nm tolerance). Useful when
  you need the full blue/red end of the bandpass.

`protect`
: Per-target equivalent of the v1.2.0 catalog-wide collision
  protection. OR'd with the catalog-wide cutoff.

**Tolerant filter** (v1.3.0 final): a required-λ range entirely
outside the current disperser's wavelength bounds is silently
filtered out at evaluator-init rather than dropping the source.
Lets users pre-stage `1.0–1.2 μm` while G395H is selected without
losing every tagged source. Partial-overlap ranges are kept — the
standard `interval_covered` check handles the achievable portion.

### Per-target source-centering override

The Constraints… popover gains a `Source centering override`
dropdown with the five canonical levels (UNCONSTRAINED →
TIGHTLY_CONSTRAINED) plus `(use global)`. Whatever you pick wins
**unconditionally** over the optimizer modal's global
Source-centering Select for that one row — even when it's laxer
than the global. A small italic line under the global Select
shows how many rows carry overrides. Round-trips through the
catalog CSV as a new `centration` column.

### Customisable stats bar + catalog hover

Settings tab adds two pop-up dialogs, `Customise stats bar…` and
`Customise catalog hover…`. Both let the user reorder / hide
fields via MultiChoice chips — drop a chip with its × to hide,
click in the dropdown to re-add, drag to reorder. Affects the top
status bar above the figure and the tooltip shown when hovering a
catalog target on the canvas.

### Draggable modals + uniform header bar

Every modal dialog (optimizer config / results / advanced, catalog
editor, per-target Constraints, Customise pickers, Overwrite
confirmation) is now a draggable card. Each carries a distinct
**light-blue header bar** with the title on the left and an
explicit ✕ close button on the right — the **header is the only
drag handle**, so form controls and SlickGrid cells in the body
stay fully interactive.

### Confirm-overwrite dialog for file writes

Catalog editor's Save-as-CSV and the Save-session button now route
through a generic overwrite-confirmation modal when the target
path already exists. Red "Overwrite" danger button + default
"Cancel"; Cancel preserves the existing file. eMPT bundle export
already creates a timestamped subdirectory each call, so no
change there.

### Compact canvas X/Y inputs

Replaced the two full-width canvas-size sliders with a compact
`row(Spinner, Spinner)` pair ("Width (X)" / "Height (Y)", 88 px
each, side-by-side). Each spinner has up/down arrows + free
typing; commit on blur/Enter/arrow click so the previous
`value_throttled` distinction is moot. The user can set any
(X, Y) — `match_aspect=True` enforces the per-pixel-square
invariant by inflating the short data axis symmetrically around
the image, so image pixels stay square and the NIRSpec FoV stays
at correct geometric ratio at any canvas aspect.

### Canvas resize feedback

Slider release now triggers the same full-page loading overlay
(gold spinner on translucent backdrop) that file loads use,
labelled "Resizing canvas…". Stays up for 1.2 s after the
Python-side resize so the browser has time to repaint the image
before the spinner fades.

### Persistent user preferences

New file `~/.vmpt/preferences.json` stores the full Settings tab
across sessions: canvas X/Y, slitlet size, snap, layer visibility,
overlay alpha + stroke per layer (all 5), stats-bar order,
catalog-hover order, help-panel visibility. Auto-saved on every
widget change via `vmpt/preferences.py` (atomic temp+rename
writes; corrupt files don't break startup). A new "Reset display
to defaults" button in the Settings tab wipes the file and
restores hard-coded defaults. Override the file location via the
`VMPT_PREFS_PATH` env var (test hook).

### Help panel default-on

The right-side help panel now opens **expanded** at first launch
(persists via prefs). One click on "Hide help" collapses it; the
new value is saved.

### hMPT / eMPT attribution rewording

README, RTD docs, and the `vmpt/optimizer.py` module docstring all
clarified: vMPT's optimizer is a lightweight Python module
**inspired by** hMPT (which is itself inspired by ESA's eMPT) —
not a direct port of either. The MSA shutter geometry, V2/V3 ↔
(s, d) transforms, gnomonic projection, and constraint machinery
were written fresh; the search algorithm is a simpler variant
than hMPT's.

### Catalog editor — z as float-sortable

The `z` (redshift) column was string-stored, so the header click
sorted lexicographically (`"10.5"` < `"2.3"`). Now stored as
float (NaN for missing), rendered with a 3-decimal HTML formatter,
and coerced back to float on cell edit — same pattern as
priority + weight in v1.1.1.

### Bug fixes

- **Image aspect ratio stretched for non-square images.** With
  the v1.3.0 independent X/Y canvas sliders, `match_aspect=True`
  silently failed whenever the canvas aspect didn't match the
  image's W:H, because the data ranges were pinned on both axes.
  Fixed by pre-computing data ranges from the canvas aspect so
  match_aspect's check passes trivially.
- **`DataRange1d.start=None` crash on canvas slider change.** A
  prior fix tried `fig.x_range.update(start=None, end=None)` to
  reset to auto-bound; Bokeh 3.7.2 rejects None on explicit
  assignment to DataRange1d.start/end. Replaced with the
  aspect-matched explicit-Float fix above.
- **Optimizer crash on multi-source catalog with `centration`
  field.** `getattr(cat, "centration", None) or []` triggered
  `bool()` on a length-N numpy array, raising
  `ValueError: The truth value of an array with more than one
  element is ambiguous`. Fixed with explicit `is None` check;
  added a regression test that drives a 5-source Catalog with
  mixed centration overrides through PointingEvaluator end-to-end.

### Editor / loader

- The catalog loader recognises the new column aliases:
  `lam_req` / `lambda_required`, `no_gap` / `gapless`,
  `extend_blue` / `bluest`, `extend_red` / `reddest`,
  `protect` / `protected`, `centration` / `source_centering`.
- `save_catalog(cat, path, include_constraints="auto|always|never")`
  is now a public function for programmatic catalog writes; the
  catalog editor's Save-as-CSV uses the same path.
- The optimizer driver passes the constraint arrays to
  `PointingEvaluator` automatically — no extra wiring per-call.

### Optimizer

- `PointingEvaluator.__init__` accepts five new constraint kwargs
  + `centration_per_target`. A new method
  `evaluate_with_reasons(...)` returns the per-pointing
  drop-reason dict keyed by `DROP_REASONS` (`collision`,
  `required_lam`, `no_gap`, `extend_blue`, `extend_red`).
- `evaluate_with_stats(...)` keeps its v1.2 4-tuple shape — the
  int it returns is now the sum of the per-reason counts.
- The Optimizer-results modal Score-cell hover tooltip breaks
  the `−K` count down by reason:

  ```
  Top 10 placed sources at this pointing:
    …
  −6 dropped:
     3× spectral collision
     2× required λ-range missing
     1× detector gap inside spectrum
  ```

### RXCJ0600 example catalog slimmed

`example_r0600/v01_fsun.cat` reduced from 28,569 → 2,000 rows
(F200W or F444W magnitude in 24–27, random sample seed=42). Loads
in ~5 s on a laptop instead of ~70 s, with no science loss for a
demo file.

### Tests

`tests/test_optimizer_constraints.py` — 30 new tests for the
constraint system (parse-string, disperser_range, interval_covered,
Catalog dataclass defaults, optimizer keys + size-mismatch, per-
constraint behaviour, drop-reason invariants).

`tests/test_optimizer.py` — 6 new tests for per-target centration
override (unconditional rule, size-mismatch raise, blank/None/
unknown fall back to global), tolerant required_lam filter, and
the multi-source numpy-truthiness regression.

`tests/test_catalog.py` — 4 new tests for the constraint CSV
round-trip (auto / always / never emission policies) and 4 for the
centration column.

**183 passed, 5 skipped** total (up from 139/4 at v1.2.2). The
skips are synthetic-geometry cases where no source happens to
land in the MSA at the test pointing.

[1.3.0]: https://github.com/fengwusun/vMPT/releases/tag/v1.3.0

---

## [1.3.0-superseded] — earlier in-flight notes

The text below was the initial scoping note for v1.3.0 (the
per-target spectral constraints work alone). Preserved verbatim
for the project log; the consolidated v1.3.0 release notes above
supersede it.

Feature: **per-target spectral constraints**.

The catalog editor gains a per-row **Constraints…** button (gray
when unset, primary-blue when ≥1 field is set). Clicking it opens a
popover that captures four independent spectral requirements for
that source. At every candidate pointing the optimizer fetches the
source's centre-shutter wavelength endpoints via
`vmpt.wavelengths.cutoffs` and drops the source if any constraint
fails — same machinery as the v1.2.0 collision drop, just with new
reason codes.

### Per-target constraints

`required_lam`
: List of `(λ_lo, λ_hi)` ranges in μm that **must** land on the
  detector for this source (the gap is excluded, so a range that
  straddles the NRS1/NRS2 detector gap fails). Empty list = no
  constraint. Editor input format: `"1.0-1.3; 1.5-1.8"`.

`no_gap`
: Boolean. When True, the NRS detector gap must not fall inside
  the source's spectrum. PRISM at non-central shutters has no gap
  → passes. H gratings have a gap inside every shutter's spectrum
  → fails.

`extend_blue`
: Boolean. When True, the shutter must reach the disperser/filter's
  MSA-wide best blue wavelength (within a 20 nm tolerance to absorb
  per-shutter wavelength-solution variation). Useful for science
  that needs the full blue end of the bandpass.

`extend_red`
: Boolean. Same on the red side.

`protect`
: Boolean. Per-target equivalent of the v1.2.0 catalog-wide
  "collision protection" cutoff. Either source making a row
  protected enables the v1.2 collision rules for that row (logical
  OR).

### Optimizer

`PointingEvaluator.__init__` accepts the five new arrays. A new
method `evaluate_with_reasons(...)` returns the per-pointing
drop-reason dict keyed by the constants in
`vmpt.optimizer.DROP_REASONS` (`collision`, `required_lam`,
`no_gap`, `extend_blue`, `extend_red`). `evaluate_with_stats(...)`
keeps its v1.2 4-tuple shape — the int it returns is now the sum of
those per-reason counts.

### Catalog loader

The CSV / FITS / ASCII loader recognises five new column aliases:
`lam_req` (or `lambda_required`, `wavelength_required`, …),
`no_gap` / `gapless`, `extend_blue` / `bluest`, `extend_red` /
`reddest`, `protect` / `protected`. Boolean values accept any of
`{1, true, yes, ✓, ✔, on}` (case-insensitive). The wavelength-range
column stores as the same `"1.0-1.3; 1.5-1.8"` string the editor
uses, parsed at load.

### Results modal

The Score-cell hover tooltip now breaks `−K` down by reason:

```
Top 10 placed sources at this pointing:
  …
−6 dropped:
   3× spectral collision
   2× required λ-range missing
   1× detector gap inside spectrum
```

The hover top-10 keeps the 🛡 prefix for protected sources from
v1.2.0.

### Persistence

The Constraints… popover's edits stay scoped to the catalog
editor's working copy. Clicking **Apply changes & close** writes
the five new fields back to the in-memory `Catalog`. **Save as
CSV** emits the new columns *only* when at least one row has a
constraint set, so v1.2.x users who never touch the popover get
the same CSV format they had before.

Session save/load via `vMPT_workspace.json` is unchanged — the
workspace JSON references catalog **paths**, not contents. To
persist constraints across sessions, the user must **Save as CSV**
before **Save session**; the saved CSV becomes the canonical
catalog file for that session bundle.

### Tests

`tests/test_optimizer_constraints.py` — 30 new tests:

  - Parse-string helpers (8 parametrised cases for `lam_req` text
    format).
  - `disperser_range` / `disperser_min_lambda` / `disperser_max_lambda`.
  - 7 `interval_covered` semantics cases.
  - `Catalog` dataclass defaults.
  - Optimizer keys + size-mismatch validation (5 cases).
  - Required-λ in/out of range (2 cases).
  - `no_gap` under H gratings.
  - `total_dropped == sum(per_reason_counts)` invariant.
  - Per-target `protect` alone enables collision rules.

**169 passed, 5 skipped** total (up from 139/4 in v1.2.2). The
skips are synthetic-geometry cases where no source happens to land
in the MSA at the test pointing.

## [1.2.2] — 2026-06-03

Packaging release: vMPT is now **pip-installable** from PyPI as
[`jwst-vmpt`](https://pypi.org/project/jwst-vmpt/).

### `pip install`

```bash
pip install jwst-vmpt
vmpt                               # opens at http://localhost:5006/app
```

The console script `vmpt` accepts the same flags as `run.sh`:
`--port`, `--fits`, `--jpg`, `--wcs`, `--catalog` (repeatable). A
new `vmpt examples download [DIR]` subcommand pulls the two example
datasets (`example_a370`, `example_r0600` — together ~64 MB) from a
GitHub release asset on demand, so the pip wheel itself stays at
~20 MB (only the required MSA grid + per-shutter dispersion table
are bundled).

### Repo restructuring

- The Bokeh app directory renamed from `app/` to `vmpt/` so it
  doubles as the Python import package. All `from app.X` imports
  rewritten to `from vmpt.X` (17 files, ~56 references).
- `data/*.npz` moved to `vmpt/data/*.npz` so the wheel ships them
  alongside the modules. Path lookups in `vmpt/msa.py` and
  `vmpt/wavelengths.py` adjusted (one fewer `parent`).
- `run.sh` updated to `bokeh serve vmpt/` for source-tree users; no
  behavioural change.
- New top-level files: `pyproject.toml` (PEP 517/518 metadata + the
  `vmpt` console script), `MANIFEST.in` (sdist completeness),
  `vmpt/cli.py` (entry point).
- `.gitignore` gains `build/`, `dist/`, `.eggs/` for the pip build
  flow.

### No behavioural changes

Tests still pass at **139 / 4 skipped**. The Bokeh app behaves
identically to v1.2.1; only the install path / directory name
changed.

## [1.2.1] — 2026-06-03

Patch release. Two real bug fixes on top of v1.2.0's collision-
protection feature, plus a substantial UI cleanup that came out of
a hands-on review.

### Collision-protection fixes

- **Row tolerance is now slitlet-aware.** v1.2.0 hard-coded
  ``|Δs| ≤ 1`` between source centres, but a 3-shutter slitlet at
  row ``s_p`` already occupies rows ``{s_p−1, s_p, s_p+1}`` and the
  user-requested rule is "no other shutter at s_p±2 either." The
  evaluator now computes two tolerances at construction time:
    - Protected slitlet ↔ stuck-open (single shutter):
      ``|Δs| ≤ half + 1``
    - Protected ↔ another slitlet (same slit_length):
      ``|Δs| ≤ 2·half + 1``
  where ``half = slit_length // 2``. So the default ``slit_length=3``
  now correctly forbids stuck-open or other slitlets at rows
  ``s_p±2``. For ``slit_length=5`` the buffer scales up to
  ``s_p±3`` (stuck-open) or ``s_p±5`` (other slitlets).
  ``SHVAL_S_TOLERANCE = 1`` is preserved as the per-individual-
  shutter constant for the live-canvas orange overlap glyph (each
  opened shutter contributes its own ±1 zone, so the visualization
  already paints the correct envelope around a multi-shutter slitlet).
- **Advanced settings modal sits above the new config modal.**
  When the optimizer config dialog was added in this release the
  Advanced settings card stayed at the same z-index, so opening
  Advanced from inside Configure showed nothing — the config card
  drew on top. Bumped Advanced backdrop / card to ``z-index`` 1001 /
  1002.

### Pointing-tab UI moved into a dialog

The Pointing tab used to stack 10+ optimizer-config widgets, and the
``Run optimization`` button slid below the fold on any window under
~1200 px tall. The whole block now lives in a centered modal:

- The Pointing tab shows a single primary ``Open optimizer…`` button.
- The modal (``opt_config_modal_card``) contains every optimizer-
  config widget plus ``Run`` / ``Cancel`` and the live status line.
- The existing progress + results modal flow (``opt_modal_card``) is
  unchanged after ``Run`` is clicked — the config card just dismisses
  itself first.
- Both the optimizer config and the catalog editor modals gained a
  top-right ``×`` dismiss button.

### Help / status text — context-aware

- The help panel on the right side of the canvas is now **collapsed
  by default**. The toggle button stays in place; one click on
  ``Show help`` restores width to its v1.2.0 size with the Quick
  guide + rotating tip. The figure uses fixed
  ``frame_width``/``frame_height`` so the canvas pixel aspect doesn't
  change when the panel collapses / expands.
- The Method dropdown's three-line Democracy / Meritocracy /
  Hierarchy blurb is hidden by default; an ``ⓘ What do these mean?``
  toggle reveals it on demand. The dropdown's own option labels
  already carry the one-line summary.
- The status line under ``Run optimization`` was always reading
  *"Load a catalog with priorities, then click Run."* even when a
  catalog with priorities was loaded. ``_refresh_opt_status_div()``
  now updates it based on (catalog presence, method, priority /
  weight column availability):
    - no catalog → ``Load a catalog (Input tab) before running.``
    - catalog + Democracy → ``Ready · N sources.``
    - catalog + Meritocracy without ``weight`` →
      ``⚠ Meritocracy needs a weight column.``
    - catalog + Hierarchy without ``priority`` →
      ``⚠ Hierarchy needs a priority column.``

### Input / MPT tabs

- All path inputs across the Input + MPT tabs now use a unified
  ``_wrap_path_picker`` helper: the path ``TextInput`` is hidden
  behind an ``Edit path`` toggle when empty, and Browse buttons are
  promoted to primary blue. The path **auto-reveals** as soon as
  it's populated (by Browse, by autoload, or by typing), so users
  always see what's loaded — only the empty default is hidden.
- The MPT tab is grouped into four sections (Import / Save / Load /
  Export) separated by dashed and solid hr dividers, so the 10+
  widgets feel like coherent blocks instead of one long column.
- Renamed the ``Setting`` tab title to ``Settings`` (singular →
  plural).

### Tests

- ``tests/test_optimizer_protection.py`` gains 4 new tests:
  parametrize over ``N ∈ {1, 3, 5}`` to pin the cached tolerances,
  and a regression that an N=3 slitlet drops at least as many
  unprotected sources as N=1 under an H grating. **139 passed, 4
  skipped** in total (up from 135 / 4 in v1.2.0).

## [1.2.0] — 2026-06-03

Feature release: **shutter collision protection in the optimizer**.

Same-row sources on the same NIRSpec detector half (Q1/Q3 → NRS1,
Q2/Q4 → NRS2) disperse onto overlapping detector pixels when their
V2 separation is smaller than the spectrum's V2 half-extent
(`app.wavelengths.v2_overlap_distance` — 35″ for PRISM, ~500″ for
the H gratings). Until now the optimizer counted both members of
every such pair as observable; the live canvas already painted the
loser orange, but the score didn't reflect that downstream
penalty. v1.2.0 wires the same collision check into the optimizer's
per-pointing scoring so the user can mark high-priority targets as
**protected** and have the optimizer steer them into rows free of
collisions.

### Optimizer core (`app/optimizer.py`)

- `PointingEvaluator` accepts new keyword args: `protect_mask`,
  `priorities`, `weights`, `disperser`, `filt`, `reason`. When
  `protect_mask` is None (the default), behaviour is identical to
  v1.1.1 — the existing 16 optimizer tests are unchanged.
- A new method `evaluate_with_stats(...)` returns the existing
  3-tuple plus the count of sources dropped by the collision rules
  at this pointing. `evaluate(...)` still returns the 3-tuple but
  its `detected` mask is now the **kept** mask (post-drop) when
  protection is configured, so callers that score via `det.sum()`
  pick up collision filtering for free.
- Three rules, applied in order at every pointing:
  1. **Protected ↔ stuck-open** — a protected source landing on a
     row colliding with any shutter flagged as stuck-open (REASON
     == 2 in the CRDS `msaoper` file) is dropped. Stuck-opens
     always disperse light onto the detector regardless of which
     slitlets the user opens, so the protected target's spectrum is
     unavoidably contaminated.
  2. **Protected ↔ protected** — within each colliding cluster the
     lowest-priority-number source wins. Ties on priority break on
     higher weight; ties on weight break on lower index (stable).
     Losers are dropped; winners continue to provide collision
     pressure on the next rule.
  3. **Protected ↔ unprotected** — every unprotected source whose
     row collides with any still-kept protected source is dropped.
- Dropped protected sources do **not** propagate collision pressure
  to rules 2/3 — if a high-priority spectrum is already
  contaminated we won't compound the loss by also blocking
  unprotected sources in the same row.

### Pointing-tab UI (`app/main.py`)

- New **"Protect spectra from collision"** group in the optimizer
  sidebar (just below the existing Priority cutoff input):
  - Checkbox: **Enable collision protection**.
  - Radio: **By priority ≤** | **By weight ≥** (mutually exclusive).
  - **Threshold** text input.
  - Live status line: e.g. *"12 protected · 240 other (G140H /
    F100LP · V2 overlap ≈ 500″)"* — updates as you toggle the
    checkbox, switch the radio, type a threshold, or change the
    current Disperser/Filter.
- `_rebuild_merged_catalog` now propagates `weight` when multiple
  catalogs are stacked — previously single-catalog mode kept weight
  (it pointed at the original `Catalog` object) but merged mode
  dropped it, so the multi-catalog "By weight ≥" rule would have
  silently selected zero sources.

### Results modal

- When protection is enabled the Score cell gains a **`−K`**
  suffix where K = number of collision-dropped sources at this
  pointing. The Score column is widened by ~36 px so the suffix
  doesn't ellipsis-truncate.
- The hover top-10 prefixes protected sources with **🛡** so the
  user can verify which sources are providing collision pressure.
  A trailing line in the tooltip explains the −K count.
- Header summary line picks up a **"🛡 collision protection ON"**
  badge with a one-line explainer.

### Tests

- New `tests/test_optimizer_protection.py` — 11 tests covering
  backwards compatibility, input validation, all three drop rules,
  cross-detector-half non-collision, stuck-open handling, and the
  invariant that protection can only reduce a pointing's score
  (never increase it). 3 tests skip gracefully when synthetic
  sources don't happen to land in the geometry the test exercises.
- Existing test suite unchanged: **135 passed, 4 skipped** total
  (up from 124/1).

### Notes / known limitations

- For the H gratings the V2 overlap distance is ~500″ — comparable
  to the full MSA — so even one protected target rules out a large
  fraction of co-observable sources. The result is physically
  truthful, not a bug; the modal shows the lower kept count so
  expectations match reality.
- Unprotected sources whose rows collide with a stuck-open shutter
  are NOT dropped from scoring (only protected sources have the
  contamination penalty applied to them). This matches the
  user-requested semantic: protection is a high-priority-only
  feature, not a universal contamination filter.

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
