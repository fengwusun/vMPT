# Exporting & sharing

vMPT writes a self-contained **bundle directory** that captures
your full session. Same writer for `Save session` and
`Export eMPT bundle`; the only difference is where the files land.

## Files written

| File | Format | Consumed by |
|---|---|---|
| `MPT_plan.json` | JSON | APT MPT importer |
| `<catalog_stem>.cat` | Text | APT MPT target list |
| `vMPT_workspace.json` | JSON | vMPT-only (full state) |
| `eMPT_pointings.txt` | Text | eMPT pipeline step 0 |
| `eMPT_observed.cat` | Text | eMPT pipeline step 1 |
| `eMPT_shutter_mask.csv` | CSV | eMPT pipeline step 2 |

`vMPT_workspace.json` is the **canonical save**: every other file
can be regenerated from it. The recipient just needs to point
`Load session` at the workspace JSON (or `MPT_plan.json` — the
sibling is auto-discovered).

## Loading into APT

1. Open APT, create or open a NIRSpec MOS observation.
2. **File → Import MPT plan** → point at the `MPT_plan.json` from
   your bundle.
3. The pointing (RA, Dec, V3 PA), the disperser/filter, and the
   open-shutter pattern all land in APT.
4. Run `MPT-Verify` in APT to confirm vMPT's slitlets match APT's
   expectations.

## Loading into the eMPT pipeline

The three `eMPT_*` files mirror the format ESA's eMPT pipeline
expects:

```bash
cd /your/empt/checkout/
./step_0_ipa.sh                          # uses eMPT_pointings.txt
./step_1_observed.sh                     # uses eMPT_observed.cat
./step_2_shutters.sh                     # uses eMPT_shutter_mask.csv
```

## Loading back into vMPT

`Load session` accepts either:

- `vMPT_workspace.json` (full state — picks, slitlet size, snap
  setting, layer toggles, undo history, custom catalog columns)
- `MPT_plan.json` (less complete — picks + pointing only; the
  sibling `vMPT_workspace.json` is auto-discovered if present)

This is the **collaboration loop**: pick → save session → share the
bundle → collaborator opens it, edits, saves their own bundle →
back to you.

## Round-trip guarantee

A `Save session` followed by `Load session` reproduces the canvas
exactly. The integration tests in `tests/test_session_io.py` and
`tests/test_end_to_end.py` pin this round trip.

## What lives outside the bundle

These stay on your machine and are **not** written into the bundle:

- The image file (FITS or JPG + WCS sidecar) — too large; vMPT
  remembers the path.
- The original catalog file — vMPT writes its own copy via the
  `.cat` exporter; the source CSV/FITS isn't bundled.
- Stuck-open shutter mask — comes from CRDS, refreshed every run.
- Operability snapshot — same.

So when a collaborator opens your bundle, they need the **same**
image file at the **same** path you used. If they're on a
different machine, they can either symlink to their local copy or
load the image themselves first, then `Load session` on top of it.
