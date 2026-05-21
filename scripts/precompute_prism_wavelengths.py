"""Pre-compute per-shutter PRISM/CLEAR detector wavelength bounds.

Runs msaviz's per-shutter dispersion integration over every MSA
shutter that has an entry in the pipeline reference table, then
extracts the four wavelength values vMPT shows in shutter tooltips:

    blue_edge   = min illuminated wavelength on either detector
    gap_lo      = max wavelength on NRS1 (just BEFORE the gap)
    gap_hi      = min wavelength on NRS2 (just AFTER the gap)
    red_edge    = max illuminated wavelength on either detector

Saved as `data/prism_cutoffs.npz` with four (4, 171, 365) float32
arrays — one slice per quadrant, with NaN where the shutter's
spectrum doesn't reach the corresponding detector. The table is
indexed by vMPT's 1-based (q, s, d) as `arr[q-1, s-1, d-1]`.

Usage:
    PYTHONPATH=/tmp/msaviz python scripts/precompute_prism_wavelengths.py

`/tmp/msaviz` is a clone of https://github.com/spacetelescope/msaviz
— used only to run the integration; vMPT does NOT depend on msaviz
at runtime, it loads the pre-computed npz.

Re-run this script if msaviz updates its reference files. The
output table is committed alongside the code.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np

# msaviz uses scipy.integrate.odeint and the pipeline reference
# FITS tables shipped in its `data/` directory. Set PYTHONPATH to a
# msaviz checkout before running.
from msaviz.msa import MSA


HERE = Path(__file__).resolve().parent.parent
OUTPUT = HERE / "data" / "prism_cutoffs.npz"


def main() -> None:
    print(f"Output: {OUTPUT}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    # Initialise msaviz for PRISM/CLEAR. Loads the dispersion polynomial
    # and the per-quadrant lookup tables of integration ICs.
    msa = MSA("clear", "prism")
    sci_lo, sci_hi = msa.sci_range
    print(f"sci_range = ({sci_lo:.3f}, {sci_hi:.3f})")

    # Output arrays — one slice per quadrant, indexed [q-1, s-1, d-1].
    NQ, NS, ND = 4, 171, 365
    blue = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
    gap_lo = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
    gap_hi = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
    red = np.full((NQ, NS, ND), np.nan, dtype=np.float32)

    # Walk each quadrant: pull the (I, J) shutters present in the
    # lookup table and run the integration in batches of 500. Batching
    # amortises Python-level overhead; ~3 ms per shutter end-to-end
    # via msaviz.__call__.
    t_start = time.time()
    BATCH = 500
    for q in range(1, 5):
        tab = msa._quadrants[q - 1]  # 0-based attribute on MSA
        cols = np.asarray(tab["I"], dtype=int)   # vMPT d (1-365)
        rows = np.asarray(tab["J"], dtype=int)   # vMPT s (1-171)
        n_q = len(cols)
        print(f"Q{q}: {n_q} populated shutters")

        for batch_start in range(0, n_q, BATCH):
            batch_end = min(batch_start + BATCH, n_q)
            cs = cols[batch_start:batch_end]
            rs = rows[batch_start:batch_end]
            # msaviz expects 0-based (quadrant, column, row) per
            # _prism_integrate's docstring (line 295), but the
            # __call__ wrapper takes them as-is and passes through.
            # We need 0-based here:
            qs_0 = np.full(len(cs), q - 1, dtype=int)
            cs_0 = cs - 1
            rs_0 = rs - 1
            coords = np.vstack([qs_0, cs_0, rs_0])

            waves = msa(coords)  # (2, N, 2048)
            # zeros means "not illuminated"; treat as NaN for clean min/max
            waves = np.where(waves == 0.0, np.nan, waves)

            with np.errstate(all="ignore"):
                nrs1_min = np.nanmin(waves[0], axis=1)
                nrs1_max = np.nanmax(waves[0], axis=1)
                nrs2_min = np.nanmin(waves[1], axis=1)
                nrs2_max = np.nanmax(waves[1], axis=1)

            # Clip to the disperser's intrinsic range.
            for arr in (nrs1_min, nrs1_max, nrs2_min, nrs2_max):
                np.clip(arr, sci_lo, sci_hi, out=arr, where=np.isfinite(arr))

            for k, (c1, r1) in enumerate(zip(cs, rs)):
                qi, si, di = q - 1, r1 - 1, c1 - 1
                v_b1, v_b2 = nrs1_min[k], nrs2_min[k]
                v_r1, v_r2 = nrs1_max[k], nrs2_max[k]
                blue[qi, si, di] = np.nanmin([v_b1, v_b2])
                red[qi, si, di] = np.nanmax([v_r1, v_r2])
                # Gap is meaningful only when BOTH detectors are lit.
                if np.isfinite(v_r1) and np.isfinite(v_b2):
                    gap_lo[qi, si, di] = v_r1
                    gap_hi[qi, si, di] = v_b2

            elapsed = time.time() - t_start
            done = batch_end + sum(
                len(msa._quadrants[qq - 1]) for qq in range(1, q)
            )
            total = sum(len(msa._quadrants[qq - 1]) for qq in range(1, 5))
            print(
                f"  Q{q} {batch_end}/{n_q} ({100*done/total:.1f}% overall, "
                f"{elapsed:.1f}s elapsed, ETA "
                f"{elapsed*(total-done)/max(done,1):.0f}s)"
            )

    print(f"Total integration time: {time.time()-t_start:.1f}s")
    n_with_gap = int(np.isfinite(gap_lo).sum())
    n_with_any = int(np.isfinite(blue).sum())
    print(f"Shutters with any spectrum: {n_with_any}")
    print(f"Shutters spanning the gap: {n_with_gap}")

    np.savez_compressed(
        OUTPUT,
        blue_edge=blue,
        gap_lo=gap_lo,
        gap_hi=gap_hi,
        red_edge=red,
    )
    print(f"Saved → {OUTPUT} ({OUTPUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
