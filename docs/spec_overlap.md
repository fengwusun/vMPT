# Spec-overlap colours

When you pick a shutter, vMPT shades nearby shutters whose spectra
would land on the same detector pixels as yours. The shading
matches APT MPT's three-colour convention:

| Colour | Name | Meaning |
|---|---|---|
| 🟪 **purple** | Mask Conflict | Two open slitlets crowd each other with no operable buffer row between them. Only the rows where they crowd — the `±2`-row window between the two slitlets — turn purple; the rest of each band stays orange. |
| 🟧 **orange** | Masked | An operable shutter whose spectrum would overlap one of YOUR open picks. Opening it would put two spectra on the same pixels. |
| 🟪 **pink** | Mask Stuck | An operable shutter whose spectrum would overlap a STUCK-OPEN shutter's dispersion. Stuck-opens disperse light unconditionally; opening the pink shutter mixes its spectrum with the stuck-open's. |
| ⬜ silver outline | Operable | Empty operable shutter. Pick away. |
| 🟥 red fill | Your pick | You opened this shutter. |
| 🟥 dark red | Stuck-open | Cannot be closed; disperses light unconditionally. |

The three contamination colours each have their own alpha slider
(**Settings → Overlay appearance**); per-polygon alpha stacks with
the number of sources contaminating that shutter, so a shutter hit
by 3 dispersion sources looks ~3× darker than one hit by 1.

## When each colour appears

### Pink (Mask Stuck)

Shown the moment you pick a disperser. Independent of your picks —
every stuck-open shutter in the operability table disperses light;
the operable shutters whose spectra would collide with that
dispersion are tinted pink so you know not to pick them.

### Orange (Masked)

Appears on every shutter you'd masking-clobber by opening one of
its row-mates within the disperser's V2 window. Re-rendered on
every shutter pick. The orange stripe is always exactly **N+2 rows**
wide for an N-shutter slitlet (your N picked rows plus a one-row
buffer above and below, capturing the ±1-row tolerance from the
NIRSpec spectral trace).

### Purple (Mask Conflict)

Two open slitlets genuinely collide. There are two cases, and both are
**bounded** — purple is only ever a re-classification of an
already-contaminated shutter, so it never paints a clean (silver) shutter and
never exceeds the orange/pink that would exist without the conflict.

**Same-quadrant** — two slitlets TOUCHING (row ranges adjacent or overlapping
with no operable row between them). Purple is bounded to the rows where they
crowd:

* Let the lower slitlet span rows `[a_lo … a_hi]` and the upper one
  `[b_lo … b_hi]`. Purple covers `[b_lo − 2 … a_hi + 2]` — the gap
  between their near edges, extended 2 rows each way — across the
  full dispersion width (every operable shutter in those rows).
* Beyond that window, the rest of each slitlet's band reverts to
  orange (Masked). So two adjacent N=3 slitlets give a 2-orange /
  4-purple / 2-orange stack, matching APT MPT.

**Cross-quadrant** — two slitlets in *different* quadrants whose spectra fold
onto the same detector region (Q1↔Q3 on NRS1, Q2↔Q4 on NRS2). These share no
MSA-row frame, so purple is bounded to the **genuine overlap**: only the
shutters contaminated by BOTH slitlets (where the two spectra actually land on
the same pixels) turn purple; shutters hit by only one stay orange. The
conflicting slitlets' own **open shutters are also purple** — each is one of
the two colliding spectra, even though the partner's light falls on
neighbouring shutters rather than the pick itself. (Closing either slitlet
clears the conflict.)

A single isolated slitlet on a clean (silver-outlined) row never
produces purple — only ever orange / pink alpha-stacking. Note that a slitlet
is a CONTIGUOUS column of open shutters: two separate clusters in the same
column (common in a mask CSV, where opens carry no target id) are treated as
two slitlets, so the empty gap between them is never masked.

**Grating diagonal-step relaxation.** A *no-buffer adjacency* (the two
slitlets exactly 1 row apart, with no real row overlap) is only a true
Mask Conflict for PRISM and for **same-column** stacking. On the grating
side, stepping the second slitlet to a *different column* (Δd ≥ 1) is a
deliberate *diagonal* step — e.g. tracking an elongated galaxy — and is
demoted from purple to **orange**.

Why any nonzero column offset suffices (not a large margin): two slitlets
at adjacent rows stay a **fixed ~1-row apart in cross-dispersion no matter
how far apart they are in columns** — both spectra share the same trace
tilt, so a column step only slides them along dispersion (detector *x*),
never across it (detector *y*). And the spectra are far longer (M ~510
col, H ~1260 col) than the MSA is wide (365 col), so the *x*-overlap never
closes within reach either. The column-offset *magnitude* therefore can't
change the (marginal, 1-row) overlap; the only physically meaningful split
is **same-column (Δd=0, spectra fully stacked → purple)** vs
**different-column (Δd ≥ 1, a deliberate step → orange)**. The threshold
constants `GRATING_ADJ_MIN_COLSEP_H` / `_M` in `vmpt/wavelengths.py` are
therefore both **1** (raise them if you want to require a wider deliberate
step). Real row overlap, same-column adjacency, and PRISM always stay
purple. This makes diagonal slit-stepping practical; the optimizer's
collision protection honours the same rule so it won't block those steps.

## How the collision check works

vMPT does a **detector-pixel** intersection check for every
candidate shutter against every open / stuck shutter, on every
state change. The check has three pieces:

1. **MSA-row check (within-quadrant)**. Two shutters at the same s
   row in the same quadrant share a detector y row — if their
   V2 separation is within `v2_overlap_distance(disperser, filter)`
   they collide. ±1-row safety buffer (`SHVAL_S_TOLERANCE = 1`).

2. **Subtractive x-range filter**. When the precomputed
   `data/dispersion_cutoffs.npz` is available (all 9 supported
   combos in the v1.3.1 wheel), the MSA-row hits are filtered
   AGAIN by checking the actual detector x-range overlap on a
   shared detector. This drops V2-distance false positives — e.g.
   G140M/F070LP at ΔV2 ≈ 84″ where the spectra are physically
   ~60″ wide on detector and don't reach each other even though
   `v2_overlap_distance = 98″` would say they could.

3. **Additive cross-quadrant detector-y check**. Catches
   cross-quadrant pairs whose spectrum y-stripes coincide at the
   x-overlap region even when the MSA s coordinate differs. Fires
   only for `q_candidate ≠ q_open`, requires x-overlap ≥ 10 px AND
   `|Δy_local| ≤ 5 px` (slit thickness). Catches e.g.
   G140M/F100LP Q4 s=34 ↔ Q2 s=33 where MSA-row alone would miss
   the collision.

   The per-shutter detector-x/y used here is bilinearly interpolated
   from a coarse 10×10 quadrant sample (rows ≈16–155, cols ≈34–331).
   Off-grid **edge** rows/cols are **linearly extrapolated**, not
   clamped — a clamp would pin every edge row to one detector-y and
   make this check fire across a whole block of edge rows that don't
   actually share a detector row (it once flooded lower Q2 with
   spurious Masked shutters whenever Q4 had opens).

## Worked example: G395M Q1 s=112 d=312 vs Q3 s=111 d=214

User-supplied test case (one of the regression set):

* Open Q1 d=312 s=112 has spectrum on NRS1
  `x=[1735, 2014]` (the blue end, ~2.87–3.41 μm), and on NRS2
  `x=[10, 883]` (the red end, ~3.71–5.27 μm) — gap-spanning.
* Q3 d=214 s=111 has spectrum almost entirely on NRS1
  `x=[824, 2014]` (~2.87–5.03 μm).

The two NRS1 footprints overlap at `x=[1735, 2014]` (about 280 px
wide). At that x range, Q1 is delivering its blue end (~2.87–3.34
μm) and Q3 is delivering its red end (~4.51–5.03 μm) — different
wavelengths from the two shutters land on the same detector pixels.
That's contamination, and vMPT flags Q3 d=214 s=111 as a
buffer-hit (orange tint) when Q1 d=312 s=112 is opened.

## Worked example: G140M/F070LP — no overlap despite ΔV2 < v2_overlap

A V2-distance-only check would over-flag this:

* Open Q3 d=208 s=108 G140M/F070LP — narrower 0.70–1.27 μm range
  clips the spectrum to NRS1 only, `x=[320, 1216]`.
* Candidate Q1 d=240 s=109 — most of its spectrum is on NRS2
  (gap-spanning); its tiny blue tip on NRS1 is at `x=[1402, 2030]`.

ΔV2 between the two shutters is +84″, well inside the 98″
`v2_overlap_distance` for G140M/F070LP — so the V2 check passes.
But the spectra physically don't overlap (gap from x=1216 to
x=1402 on NRS1). The **subtractive x-range filter** drops this
candidate and the orange / purple band correctly stops at the
detector gap.

The same shutter pair under G140M/F100LP (wider 0.97–1.89 μm
range) DOES overlap — the longer spectrum closes the gap. vMPT
flags it as DIRECT.

## Behavior summary by disperser / filter

| Combo | v2_overlap | Notes |
|---|---|---|
| PRISM / CLEAR | 32″ | Smallest extent. Same-row pairs within ~30″ flag. |
| G140M / F070LP | 98″ | Narrow filter clips spectrum; many same-row pairs at large ΔV2 correctly don't flag. |
| G140M / F100LP | 109″ | Full M-grating extent; same-row pairs across most of a quadrant flag. |
| G235M / F170LP | 110″ | |
| G395M / F290LP | 103″ | |
| G140H / F070LP | 185″ | H-gratings span multiple detector regions. |
| G140H / F100LP | 307″ | Largest extent. Many cross-quadrant pairs flag. |
| G235H / F170LP | 300″ | |
| G395H / F290LP | 281″ | |

Per-combo numbers measured from direct
`jwst.assign_wcs.slit_frame → detector` traces. See
`tests/test_wavelengths.py` for the regression tests.

## Tilt: why the band stays N+2 rows wide

NIRSpec spectral traces actually drift slightly in the
cross-dispersion direction as you move along dispersion. If
naively applied to the row check, the orange band's union over a
full v2_overlap window can reach N+3 or N+4 rows wide at the
extremes — visually distracting (users reported this on the
v1.3.0 release).

Since v1.3.1, `row_offset` is clamped to 0: the band stays at
exactly s_open ± (N//2 + 1) rows across the entire dispersion
range. Trade-off: at very-far-d candidates where the actual
spectrum drift exceeds 0.5 rows, the rendered band doesn't track
the spectrum's physical row — but the dropped contamination at
the wings is small (the spectrum has rolled off by then).

The tilt-slope grid is still shipped in
`data/dispersion_cutoffs.npz` for the additive cross-quadrant
detector-y check (which DOES use the slope to compute local y at
the x-overlap midpoint) and for any future re-introduction of
tilt rendering.
