# Multiple MPT configurations

vMPT can plan up to **five MPT configurations** — each an independent
pointing plus its own set of open shutters, exactly like APT/MPT, where
a single MSA program can carry several configs. (Two configs arrived in
**v1.4.0**; the cap rose to five in **v1.5.0**.) Each config is a
separate exposure: their spectra never collide with each other (only
within a config), and a source can be split across them or deliberately
repeated.

## The configs block (Pointing tab)

Beneath **Open optimizer…** you'll find:

- **Number of configs** — a spinner, 1 (default) to 5. Bumping it
  creates the extra, empty configs, each starting at the current
  pointing. Dropping the count back down hides the higher configs but
  keeps their picks in memory.
- **Working on** — which config you're editing. Manual shutter clicks,
  the spec-overlap overlay, and the live status bar all act on the
  **active** config only. Switching configs swaps the pointing widgets
  to that config's saved RA/Dec/PA and redraws its shutters.
- **View MPT catalog…** — a read-only table of every selected source
  (see [below](#viewing-the-selected-sources)).

### Which config am I editing?

In multi-config mode the active config is shown three ways, all sharing
one colour per config (**1 = blue, 2 = magenta, 3 = green, 4 = orange,
5 = violet**):

- a bright **✎ CONFIG k / N** chip pinned to the always-visible top bar.
  The chip is also a **one-click switcher** — click it to cycle to the
  next config without leaving the canvas (the *Working on* dropdown stays
  in sync);
- an *Editing Config k of N* banner in the Pointing tab;
- on the image, each config's **MSA quadrant boundary** is drawn in its
  own colour — the active config as a **solid** grid, every idle config
  as a **faint dashed outline** (boundary only, no fill) at *its own*
  pointing, so you can see where the idle exposure(s) sit.

Each config remembers its own pointing, so you can:

- put configs at the **same** pointing but pick disjoint shutters
  — the canonical way to observe sources whose spectra would
  collide if opened together (split them across exposures); or
- put them at **different** pointings to cover more of the field.

The disperser/filter is shared across configs in this release.

## Auto-producing all configs with the optimizer

With **Number of configs = N** (up to 5), one click of **Run
optimization** produces a complete N-config plan using a
*sequential-greedy* search:

1. Optimize **Config 1** over the full candidate pool.
2. Charge Config 1 #1's observed sources against their
   max-observation budget.
3. Optimize **Config 2** over whatever sources still have budget,
   then **Config 3** over the residual after that, … and so on
   through Config N.

The results modal then shows a **single combined table** — one row per
*plan* (Config 1 #k + … + Config N #k), rendered as one colour-coded
line per config that all share a **combined score** on the left, with
each config's own score and ΔRA/ΔDec/ΔPA on its line. Above it sits a
**✔ Apply recommended plan #1** button (fills every config at once), and
each row has its own **Apply** so you can pick a different set.

Hover any **Score** to see the top placed sources and, where a config
dropped some, the breakdown by reason. A later config's line often shows
a "**−K**" suffix — those are sources skipped because an earlier config
already observed them under the max-configs cap; a one-line legend
appears above the table whenever any plan drops sources.

> Running five passes is practical because the per-trial shutter
> assignment uses a polynomial inverse of the MSA map (v1.5.0),
> ~20–40× faster than the previous interpolation — see the
> [optimizer notes](optimizer.md).

All three methods behave exactly as designed — Democracy maximizes the
count, Meritocracy the summed weight, Hierarchy the priority tiers — the
budget just removes already-observed sources from the pool before
scoring.

## The max-configs cap (avoiding duplicate pointings)

If a source may be observed in unlimited configs, nothing stops Config 2
from re-picking Config 1's best targets — which can hand you **two
identical pointings**. To prevent that, cap how many configs a source
may appear in:

- **Globally** — *Max configs per source* lives in the optimizer's
  **Advanced settings** (next to the rarely-used *Priority cutoff*). It
  **defaults to 1**, so out of the box each source is observed in at most
  one config and Config 2 must find *new* targets. Blank = unlimited.
- **Per source** — three ways, each overriding the global default:
  - the **Constraints…** popover in the catalog editor (a *Max MPT
    configs* selector);
  - a `max_configs` column in the catalog CSV (aliases: `max_obs`,
    `n_configs`, …);
  - a **bulk rule** in the catalog editor — *Set max configs = N for
    sources where `<condition>`*, where the condition is a sandboxed
    boolean expression over any catalog column, e.g.
    `(mag_f444w > 27) & (z > 6)` (functions `abs`, `log10`, `sqrt`,
    `isin`; combine with `&` / `|` / `~`). The expression is
    syntax-checked — unknown columns and bad syntax are reported inline —
    *before* anything changes, and only matching rows are updated. Use it
    for "observe my faint, high-z sources in both configs" cuts.

## Viewing the selected sources

**View MPT catalog…** (Pointing tab) opens a read-only table of every
selected source across all live configs:

| Column | Meaning |
|---|---|
| Cfg | which config — colour-coded (1 blue, 2 magenta, 3 green, 4 orange, 5 violet) to match the chip |
| Source ID | catalog id |
| RA / Dec / Pri / Weight / Mag / z | from the loaded catalog |
| λ_blue / λ_red | the spectrum's on-detector wavelength coverage (μm) for this source's shutter |
| Gap (μm) | the NRS1/NRS2 detector-gap range inside the spectrum (`lo–hi`), or `none` when the spectrum is gap-free |
| Q / s / d | the shutter the source sits in |
| Label | original catalog id / label |

The λ coverage + gap are computed per shutter under the **current
Disperser / Filter** when you open the viewer (a brief loading spinner
covers the compute). The numeric columns (Pri, Weight, Mag, z, λ_blue,
λ_red, RA/Dec, Q/s/d) are stored as numbers, so they **sort
numerically** — `2` before `15`, not the string order `15` before `2`.
A **Show columns** chip picker hides/reorders columns, and rows are
drag-reorderable.

A source that lands in more than one config appears once per config.
**Save as CSV** writes the table out for sharing. The viewer is empty
until you pick a shutter or run the optimizer.

## Multi-source shutters

Two very close sources can share a single shutter. vMPT records **both**
— the slitlet's primary id (the anchor) plus a full `target_ids` list
per shutter. Every selected source, including a co-shutter second
source, reaches:

- the MPT catalog viewer (one row each),
- the eMPT `observed_targets.cat` (one row each),
- the `sourceIds` of the exported plan's config block.

## Export layout

- **`<catalog>_MPT_plan.json`** — one config block per config, each with
  its own pointing and `sourceIds`. `allowMultiSourceShutters` is true.
- **`<catalog>_APT_catalog.cat`** — the APT source catalog.
- **`README.md`** — the APT/MPT import steps, filled in for this bundle.
- **`vMPT_workspace.json`** — round-trips the whole multi-config +
  multi-source state, so **Load session** restores every config exactly.
- **eMPT bundle** — the top-level files describe the active config;
  with more than one config, each also gets a `config_N/` subfolder
  (its own `observed_targets.cat` + `shutter_mask.csv` +
  `pointing_summary.txt`), since eMPT artifacts are per-pointing.

See **[Exporting & sharing](exporting.md)** for the full file list and
the step-by-step APT/MPT import.

Single-config, single-source exports are byte-compatible with v1.3.x.
