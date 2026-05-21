"""Pre-compute per-shutter detector wavelength bounds for ALL
NIRSpec disperser/filter combinations.

Runs spacetelescope/msaviz's per-shutter dispersion integration over
every MSA shutter for each (disperser, filter) pair vMPT supports,
then extracts the four wavelength values vMPT shows in shutter
tooltips:

    blue_edge   = min illuminated wavelength on either detector
    gap_lo      = max wavelength on NRS1 (just BEFORE the gap)
    gap_hi      = min wavelength on NRS2 (just AFTER the gap)
    red_edge    = max illuminated wavelength on either detector

Saved as `data/dispersion_cutoffs.npz` with one (4, 171, 365)
float32 array per combo per quantity:

    {DISPERSER}_{FILTER}_blue_edge
    {DISPERSER}_{FILTER}_gap_lo
    {DISPERSER}_{FILTER}_gap_hi
    {DISPERSER}_{FILTER}_red_edge

Values are NaN where the shutter's spectrum doesn't reach the
corresponding detector (e.g. Q3/Q4 PRISM shutters have no gap; some
H-grating edge shutters don't reach a particular detector).

Usage:
    PYTHONPATH=/tmp/msaviz python scripts/precompute_prism_wavelengths.py

`/tmp/msaviz` is a clone of https://github.com/spacetelescope/msaviz
— used only to run the integration; vMPT does NOT depend on msaviz
at runtime, it loads the pre-computed npz.

PRISM integrates per shutter (slow; ~12 minutes for 194k shutters);
gratings integrate in batches per quadrant (very fast; under a
minute for all 8 grating combos combined).

Re-run this script if msaviz updates its reference files. The
output table is committed alongside the code.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

# msaviz uses scipy.integrate.odeint and the pipeline reference
# FITS tables shipped in its `data/` directory. Set PYTHONPATH to a
# msaviz checkout before running.
from msaviz.msa import MSA


HERE = Path(__file__).resolve().parent.parent
OUTPUT = HERE / "data" / "dispersion_cutoffs.npz"

# (vMPT disperser, vMPT filter, msaviz disperser, msaviz filter).
# msaviz uses lowercase names everywhere; vMPT uses uppercase.
COMBOS = [
    ("PRISM", "CLEAR",  "prism", "clear"),
    ("G140M", "F070LP", "g140m", "f070lp"),
    ("G140M", "F100LP", "g140m", "f100lp"),
    ("G235M", "F170LP", "g235m", "f170lp"),
    ("G395M", "F290LP", "g395m", "f290lp"),
    ("G140H", "F070LP", "g140h", "f070lp"),
    ("G140H", "F100LP", "g140h", "f100lp"),
    ("G235H", "F170LP", "g235h", "f170lp"),
    ("G395H", "F290LP", "g395h", "f290lp"),
]

NQ, NS, ND = 4, 171, 365


def _shutter_bounds_per_quadrant(msa: MSA, q1based: int) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray,
]:
    """Return four (NS, ND) arrays of (blue, gap_lo, gap_hi, red)
    for every shutter in the given quadrant of the supplied MSA.
    NaN cells = the spectrum doesn't reach that detector for that
    shutter.
    """
    sci_lo, sci_hi = msa.sci_range

    # Build the (q, col, row) coordinate stack for every shutter.
    # msaviz expects 0-based indices.
    cs, rs = np.meshgrid(np.arange(1, ND + 1), np.arange(1, NS + 1), indexing="xy")
    cs = cs.ravel()
    rs = rs.ravel()
    N = cs.size
    coords = np.vstack([
        np.full(N, q1based - 1, dtype=int),
        cs - 1,
        rs - 1,
    ])

    blue = np.full((NS, ND), np.nan, dtype=np.float32)
    gap_lo = np.full((NS, ND), np.nan, dtype=np.float32)
    gap_hi = np.full((NS, ND), np.nan, dtype=np.float32)
    red = np.full((NS, ND), np.nan, dtype=np.float32)

    if msa.disperser == "prism":
        # PRISM only integrates for the shutters in its lookup table;
        # missing shutters never produce a spectrum.
        tab = msa._quadrants[q1based - 1]
        pop_d = np.asarray(tab["I"], dtype=int)
        pop_s = np.asarray(tab["J"], dtype=int)
        pop_idx = (pop_d - 1) + (pop_s - 1) * ND
        mask = np.zeros(N, dtype=bool)
        mask[pop_idx] = True
        coords = coords[:, mask]
        if coords.shape[1] == 0:
            return blue, gap_lo, gap_hi, red

    # Batch the integration. PRISM is per-shutter inside msaviz so
    # batching it large is fine; gratings batch the whole quadrant
    # at once in a single odeint.
    BATCH = 5000 if msa.disperser == "prism" else N
    keep_idx_cs = coords[1] + 1  # 1-based d
    keep_idx_rs = coords[2] + 1  # 1-based s
    for start in range(0, coords.shape[1], BATCH):
        end = min(start + BATCH, coords.shape[1])
        chunk = coords[:, start:end]
        waves = msa(chunk)  # (2, n, 2048)
        # For gratings: out-of-range wavelengths are still "integrated"
        # values but they don't correspond to real detector pixels.
        # Clip by sci_range, then derive min/max.
        # For PRISM: zeros mean "not illuminated".
        waves = np.where(waves == 0.0, np.nan, waves)
        in_range = (waves >= sci_lo) & (waves <= sci_hi)
        waves = np.where(in_range, waves, np.nan)

        with np.errstate(all="ignore"):
            nrs1_min = np.nanmin(waves[0], axis=1)
            nrs1_max = np.nanmax(waves[0], axis=1)
            nrs2_min = np.nanmin(waves[1], axis=1)
            nrs2_max = np.nanmax(waves[1], axis=1)

        ds_chunk = keep_idx_cs[start:end]
        ss_chunk = keep_idx_rs[start:end]
        for k in range(end - start):
            si = ss_chunk[k] - 1
            di = ds_chunk[k] - 1
            v_b1, v_r1 = nrs1_min[k], nrs1_max[k]
            v_b2, v_r2 = nrs2_min[k], nrs2_max[k]
            # blue = the bluer of the two detector mins; red = redder of maxes.
            b_candidates = [x for x in (v_b1, v_b2) if np.isfinite(x)]
            r_candidates = [x for x in (v_r1, v_r2) if np.isfinite(x)]
            if b_candidates:
                blue[si, di] = min(b_candidates)
            if r_candidates:
                red[si, di] = max(r_candidates)
            # Gap is meaningful only when BOTH detectors carry a real
            # piece of the spectrum.
            if np.isfinite(v_r1) and np.isfinite(v_b2):
                gap_lo[si, di] = v_r1
                gap_hi[si, di] = v_b2

    return blue, gap_lo, gap_hi, red


def main() -> None:
    print(f"Output: {OUTPUT}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    out: dict[str, np.ndarray] = {}
    total_t0 = time.time()
    for vmpt_disp, vmpt_filt, m_disp, m_filt in COMBOS:
        t0 = time.time()
        msa = MSA(m_filt, m_disp)
        sci_lo, sci_hi = msa.sci_range
        print(f"\n{vmpt_disp}/{vmpt_filt}  sci=({sci_lo:.3f}, {sci_hi:.3f})")
        blue = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
        gap_lo = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
        gap_hi = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
        red = np.full((NQ, NS, ND), np.nan, dtype=np.float32)
        for q in range(1, 5):
            b, glo, ghi, r = _shutter_bounds_per_quadrant(msa, q)
            blue[q - 1] = b
            gap_lo[q - 1] = glo
            gap_hi[q - 1] = ghi
            red[q - 1] = r
        n_any = int(np.isfinite(blue).sum())
        n_gap = int(np.isfinite(gap_lo).sum())
        elapsed = time.time() - t0
        print(f"  populated: {n_any}, gap-spanning: {n_gap}, "
              f"time: {elapsed:.1f}s")

        key = f"{vmpt_disp}_{vmpt_filt}"
        out[f"{key}_blue_edge"] = blue
        out[f"{key}_gap_lo"] = gap_lo
        out[f"{key}_gap_hi"] = gap_hi
        out[f"{key}_red_edge"] = red

    print(f"\nTotal: {time.time()-total_t0:.1f}s")
    np.savez_compressed(OUTPUT, **out)
    print(f"Saved → {OUTPUT} ({OUTPUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
