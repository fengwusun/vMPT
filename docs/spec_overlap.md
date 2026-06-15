# Spec-overlap colours

When you pick a shutter, vMPT shades nearby shutters whose spectra
would land on the same detector pixels as yours. The shading
matches APT MPT's three-colour convention:

| Colour | Name | Meaning |
|---|---|---|
| 🟪 **purple** | Mask Conflict | An open shutter (yours OR stuck-open) sitting in a touching-collision dispersion band. Chain-propagates downstream. |
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

Reserved for cases where two open slitlets are TOUCHING — their
row ranges are adjacent with no operable row between them. In that
event:

* The touching shutters on each slitlet render purple.
* The full dispersion band of each touching slitlet renders purple
  ("chain propagation"). Opening a third shutter downstream of a
  touching pair will also light up purple.

A single isolated slitlet on a clean (silver-outlined) row never
produces purple — only ever orange / pink alpha-stacking.

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
