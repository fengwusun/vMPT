# Changelog

All notable changes to vMPT are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
project adheres to [Semantic Versioning](https://semver.org/).

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
