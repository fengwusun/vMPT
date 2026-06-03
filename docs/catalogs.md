# Catalogs

vMPT accepts target lists in three formats:

- **CSV** with a header row (`id,ra_deg,dec_deg,…`)
- **ASCII** with whitespace separation (any [astropy.io.ascii](https://docs.astropy.org/en/stable/io/ascii/)
  format works)
- **FITS** binary table (`hdul[1].data`)

## Required columns

Two columns are mandatory: **RA** and **Dec** in degrees. Everything
else is optional.

## Loose column matching

vMPT identifies the canonical columns by **case-insensitive
substring**, trying multiple aliases:

| Canonical | Accepted aliases (case-insensitive) |
|---|---|
| ID | `id`, `source_id`, `name`, `no_cat`, `cat_id` |
| RA (deg) | `ra`, `ra_deg`, `alpha`, `alpha_j2000` |
| Dec (deg) | `dec`, `dec_deg`, `delta`, `delta_j2000` |
| Priority | `priority`, `p`, `pri`, `pclass`, `priority_class` |
| Weight | `weight`, `w`, `wt`, `weights` |
| Magnitude | `mag`, `magnitude`, `f150w`, … |
| Redshift | `z`, `zspec`, `zphot`, `redshift` |
| Label | `label`, `name`, `title`, `comment` |

So `Right_Ascension_deg` matches RA, `Dec_J2000` matches Dec, and
so on.

## ID resolution + the mod-10⁷ rule

When an ID is **purely numeric** and fits in an int32, vMPT keeps
it verbatim. When it's larger than 10⁷ (which happens with
HST/SExtractor catalogs that carry the source's pixel position
inside the ID), vMPT internally tags shutters using `id % 10_000_000`
so labels don't overflow the eMPT exporter's int32 field.

String IDs (e.g. `"RJ0600-12345678-P0"`) are kept verbatim through
the entire pipeline.

## Multi-catalog stacking

Click `Add` to layer multiple catalogs. Each gets its own colour
(yellow / cyan / magenta / …) and a toggle row in the sidebar so
you can show/hide individually. The merged catalog is what the
optimizer sees.

```{note}
The `weight` column is propagated correctly across multiple
catalogs from **v1.2.1+**. If you stack older catalogs without a
`weight` column, the merged array fills with NaN — Meritocracy
mode then silently treats them as weight 0.
```

## Editing in-app

Click `Edit catalog…` in the Input tab. The pop-up table lets you:

- Sort by any column (header click). Priority + Weight sort
  numerically.
- Double-click any cell to edit.
- Click 🗑️ to delete a row.
- Add custom columns via the input + `Add column` button.
- Compute one of weight/priority from the other (`Compute w from p` /
  `Compute p from w`).
- Edit per-target spectral constraints via the **Constraints…**
  button at the end of each row (v1.3.0+, see below).
- Save as CSV (standalone copy) or `Apply changes & close` (commit
  to the in-memory catalog).

### Per-target spectral constraints (v1.3.0+)

Every row has a `Constraints` button at the right edge of the
editor table. It's **gray** when no constraints are set, **primary
blue** when at least one is. Clicking it opens a popover with four
independent toggles that govern how the optimizer treats that
specific source:

`Required λ ranges`
: Text input with the format `"1.0-1.3; 1.5-1.8"` (μm,
  semicolon-separated). At every candidate pointing the optimizer
  drops the source unless **every** listed interval lands fully on
  the detector (the NRS1/NRS2 detector gap is excluded — a range
  bisected by the gap fails). A yellow inline warning surfaces if
  any range falls outside the current Disperser/Filter — the save
  is still accepted in case you're pre-staging for a future
  disperser.

`Forbid detector gap inside spectrum`
: Drops the source if the centre shutter's spectrum has the
  NRS1/NRS2 gap inside `[λ_blue, λ_red]`. For PRISM at non-central
  shutters there's typically no gap → kept. For the H gratings the
  gap is always present → fails.

`Extend to bluest λ of disperser`
: Drops the source if its shutter's `λ_blue` is more than 20 nm
  above the disperser's MSA-wide best blue. Useful when you need
  the full blue end of the bandpass (it forces the optimizer to
  pick a V2 position where the spectrum isn't truncated on the
  blue side).

`Extend to reddest λ of disperser`
: Same logic, red side.

`Protect this source from spectral collision`
: Same as the v1.2.0 catalog-wide "collision protection" cutoff,
  but applied per-row. Either source (the per-target flag OR a
  matching row from the catalog-wide cutoff in the optimizer
  modal) makes a target collision-protected.

When you click **Apply** the popover's values get written into the
editor's working copy and a single undo entry is pushed. The Edit…
button colour flips to blue for that row.

**Persistence**: per-target constraints survive `Save as CSV`
(extra columns `lam_req`, `no_gap`, `extend_blue`, `extend_red`,
`protect` are emitted when at least one row has a constraint set —
otherwise the CSV stays unchanged for backward compatibility). To
preserve constraints across `Save session` / `Load session`,
**Save as CSV first** so the workspace JSON's catalog-path
reference points at a CSV that includes the constraints.

## The optimizer's view

Different methods need different columns:

| Method | Required column |
|---|---|
| Democracy | none |
| Meritocracy | `weight` |
| Hierarchy | `priority` |

The status line under `Run optimization` will tell you if a required
column is missing.

## Visualising

The canvas paints targets as **coloured circles** (per-catalog
colour). When a slitlet is opened on a target's shutter, the
circle turns **green** to mark it as observed.

For very dense catalogs (>10k targets) the circles can overwhelm
the underlying image. Drop the layer alpha in Settings → Overlay
appearance, or use the `Show priority class ≤` / `Show mag ≤`
filters in the Input tab.
