# Overnight build log — 2026-05-11 → 2026-05-12

User signed off at ~22:00. Target: deliver M0–M4 + M6 first draft by morning, M5 analytic model included.

## Status

- [x] M0 scaffolding — git init, dirs, eMPT clone, MSA grid copied, bokeh installed
- [ ] M0 `coords.py` + sanity reproduction
- [ ] M1 Bokeh app: FITS + JPG-sidecar loading + MSA overlay at any APA_V3
- [ ] M2 target catalog overlay
- [ ] M3 operability mask (CRDS download from `/Users/sunfengwu/crds_cache/`)
- [ ] M4 hand-picking + 3-shutter slitlets + conflict detection
- [ ] M5 analytic wavelength cutoffs (per-grating constants, no CRDS)
- [ ] M6 eMPT triple writer + byte-diff against `refs/eMPT_v1/trial_00_ref/...`
- [ ] Reviewer pass

## Decisions made tonight

- **M5 simplification (user, 22:30)**: analytic dispersion model from shutter V2 offset; skip `jwst.assign_wcs`. See PLAN.md M5.
- **eMPT shutter_mask.csv format (from reference file inspection)**: 1 padded header line + 730 × 342 grid of `{x,s,1,0}` cells. CSV "row" = dispersion `d` (1–365 = Q1/Q2, 366–730 = Q3/Q4). CSV "col" = shutter row `s` (1–171 = Q1/Q3, 172–342 = Q2/Q4). To be verified against Fortran writer.
- **Stack**: Bokeh server from M1 onward (no Streamlit step).

## Subagent assignments

A. Core data layer + analytic wavelengths (`coords.py`, `msa.py`, `wavelengths.py`, tests)
B. Image/catalog I/O (`image_io.py`, `catalog.py`, tests) — handles FITS, JPG+sidecar
C. eMPT export (`empt_io.py` + tests against `refs/eMPT_v1/trial_00_ref/...`)

Main/integration in `app/main.py` (Bokeh server) — done by lead.

Reviewer pass at the end.

## Open questions to surface to user when they're back

- Per-grating dispersion constants (λ_fid, dλ/dV2, detector boundaries) — sourcing from JDox; user should sanity-check values once.
- Slitlet height default — using **3** (s−1, s, s+1) with sidebar override.
