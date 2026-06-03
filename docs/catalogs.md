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
- Save as CSV (standalone copy) or `Apply changes & close` (commit
  to the in-memory catalog).

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
