# MSA pointing optimizer

Reachable from the **Pointing** tab → `Open optimizer…`. Searches
for an (RA, Dec, V3 PA) that maximises observations of catalog
sources in operable, well-centred MSA shutters.

Algorithm: a 3-axis grid search followed by
`scipy.optimize.differential_evolution` polish of the top-N
candidates. `vmpt/optimizer.py` is a lightweight, self-contained
Python module inspired by
[**hMPT**](https://github.com/zihaowu-astro/hMPT) (Z. Wu et al.,
CfA/Harvard), which is itself inspired by ESA's eMPT
(Bonaventura et al. 2023). vMPT is **not** a direct port of
either — the MSA geometry, coordinate transforms, and constraint
machinery are independently implemented, and the search algorithm
is a simpler version than hMPT's.

```{note}
**Speed (v1.5.0).** Every trial pointing has to map each source's
(ax, ay) position to a shutter. That map was previously a
`CloughTocher2D` interpolation evaluated afresh on all ~10⁴ trial
pointings — the dominant cost. It is now a degree-4 polynomial fit
once per quadrant from the MSA grid (in normalised coordinates),
which reproduces the shutter assignment essentially exactly
(max residual 2×10⁻⁴ shutters) while running **~20–40× faster**.
This is what makes the multi-pass [N-config](multiconfig.md) search
practical — a five-config plan drops from minutes to ~20 s. The
reference interpolation is kept as an automatic fallback.
```

## Method

Three modes, picked from the **Method** dropdown:

::::{grid} 1 1 1 3
:gutter: 2

:::{grid-item-card} Democracy
Maximises raw count. Every catalog source counts the same. No extra
columns needed.
:::

:::{grid-item-card} Meritocracy
Maximises **Σ weight** of placed sources. Requires a populated
`weight` column. One source with weight 10 outranks five sources
with weight 1.
:::

:::{grid-item-card} Hierarchy
Strict priority tiers (eMPT-style). Optimises top-tier count
first; ties broken by next tier; and so on. A higher-priority
source is **never** traded for any number of lower-priority ones.
Requires a populated `priority` column.
:::
::::

The catalog editor's `Compute w from p` and `Compute p from w`
buttons let you derive one column from the other so you can switch
modes without re-annotating.

## Inputs

ΔRA / ΔDec / ΔPA
: Search box around the current pointing, in arcseconds (for ΔRA /
  ΔDec) and degrees (for ΔPA). **ΔX = 0 freezes that axis** — e.g.
  ΔPA = 0 to keep the current roll and only search RA/Dec.

Refine top N
: How many of the best grid candidates get a DE polish. Default 10.

Source centring
: APT-style buffer (UNCONSTRAINED → TIGHTLY_CONSTRAINED). Tighter
  classes can only **reduce** the number of successful placements.

### Advanced settings

Behind `Advanced settings…`:

Grid n_RA / n_Dec / n_PA
: Resolution of the brute-force grid. Default 20³ = 8 000 pointings.

DE max iter
: Max iterations per `differential_evolution` polish.

Objective
: `number` (count) vs `flux`-weighted. Superseded by the Method
  dropdown for normal use; kept for back-compat.

Source σ (arcsec)
: Gaussian PSF σ used by the throughput integration.

APT θ (DVA, deg)
: APT-style differential velocity aberration. Default 90 = no shift.

#### Source selection

Priority cutoff ≤
: Restrict the optimizer to catalog rows with `priority ≤ X` (e.g.
  P0/P1 first; do fillers by hand later). Different from the
  collision-protect threshold below. Blank = all rows.

Max configs per source
: How many MPT configs may observe each source. **Defaults to 1**, so
  in [multi-config mode](multiconfig.md) Config 2 must find sources
  Config 1 didn't take (no duplicate pointings). Blank = unlimited; a
  per-source override (Constraints… popover, `max_configs` column, or
  the catalog-editor rule) wins over this default. No effect with a
  single config.

## Shutter-collision protection

**v1.2.0+.** Optional — opt in via the `Enable collision protection`
checkbox.

When enabled, you mark a subset of catalog sources as "protected"
(by priority ≤ X or weight ≥ Y, mutually exclusive). For every
candidate pointing the optimizer then drops:

1. **Protected ↔ stuck-open** — a protected source on a row colliding
   with any stuck-open (CRDS `msaoper` REASON==2) shutter. Its
   spectrum is unavoidably contaminated.
2. **Protected ↔ protected** — within each colliding cluster, the
   lowest-priority-number source wins. Ties on priority break on
   higher weight; ties on weight break on lower source index.
3. **Protected ↔ unprotected** — every unprotected source whose
   row collides with any still-kept protected one.

Two sources "collide on the same detector row" iff:

- They lie on the **same detector half** (Q1+Q3 → NRS1, Q2+Q4 →
  NRS2; cross-half pairs image onto different detectors).
- Their |Δs| (row separation) is within the **slitlet-aware
  tolerance** (v1.2.1+): `half + 1` for stuck-open, `2·half + 1`
  for two slitlets, where `half = slit_length // 2`.
- Their V2 separation is below
  {py:func}`vmpt.wavelengths.v2_overlap_distance` — 35″ for PRISM,
  ~500″ for the H gratings.

### Score display

When protection or any per-target constraint is active, the results
table's Score cell appends `−K` where K = total sources dropped at
that pointing. Hover the cell to see:

- The top 10 placed sources, with **🛡** prefixing rows in the
  protected set.
- A breakdown of the K drops by reason (v1.3.0+):

```text
−6 dropped:
   3× spectral collision
   2× required λ-range missing
   1× detector gap inside spectrum
```

The five reason codes match the constants in
{py:data}`vmpt.optimizer.DROP_REASONS`:
`collision`, `required_lam`, `no_gap`, `extend_blue`, `extend_red`.

### A caveat for H gratings

For G140H / G235H / G395H the V2 overlap distance is ~500″ — wider
than the MSA itself — so even one protected target eliminates a
large fraction of co-observable sources. That's the truthful
answer; the modal shows a lower kept count so expectations match
reality.

## Per-target spectral constraints (v1.3.0+)

In addition to the catalog-wide "Protect spectra from collision"
toggle, each catalog row can carry its own per-target spectral
constraints via the **Constraints…** button in the catalog editor.
The four constraint types are documented in [Catalogs ›
Per-target spectral constraints](catalogs.md#per-target-spectral-constraints-v130).
Briefly:

- **Required λ ranges** — list of `(λ_lo, λ_hi)` ranges in μm that
  must land on the detector (gap-excluded).
- **No detector gap** — the NRS1/NRS2 gap may not fall inside the
  spectrum.
- **Extend to bluest / reddest** — the centre-shutter spectrum must
  reach the disperser/filter's MSA-wide best blue / red wavelength.
- **Protect** — per-target equivalent of the v1.2.0 catalog-wide
  cutoff. Either makes a row collision-protected.

The optimizer evaluates all per-target constraints **after** the
v1.2.0 collision rules, so a source can only be dropped once. The
results modal's `−K` breakdown reports which constraint type
caused the drop.

## Performance

A typical run (~500 sources, 20³ = 8 000 grid pointings, top-10 DE
refinement) finishes in **5–15 s** on a modern laptop. Grid
dominates; per-pointing eval is ~1 ms after the inverse MSA
mapping's Delaunay triangulation is cached on first call.

## Programmatic use

The optimizer is exposed as a regular Python module — you don't
need the Bokeh UI to call it:

```python
from vmpt.optimizer import PointingEvaluator, grid_search, refine_top
import numpy as np

ra  = np.array([53.16, 53.17, 53.18])      # deg
dec = np.array([-27.78, -27.79, -27.80])
pri = np.array([1, 1, 2], dtype=float)

ev = PointingEvaluator(
    ra, dec,
    centration="UNCONSTRAINED",
    slit_length=3,
    protect_mask=(pri <= 1),
    priorities=pri,
    weights=np.ones(3),
    disperser="PRISM", filt="CLEAR",
)
grid = grid_search(ev, 53.17, -27.79, 0.0, dra_arcsec=30, ddec_arcsec=30, dpa_deg=15)
top  = refine_top(ev, grid, n_top=5)
print(top["score"][:5], top["ra"][:5], top["dec"][:5], top["pa"][:5])
```

See [API reference](api.md) for the full signature.
