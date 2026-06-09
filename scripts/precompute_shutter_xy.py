"""Faster x/y-range precompute for the spec-overlap detector-pixel
check.

Same 10×10 grid per quadrant as precompute_trace_tilt.py but WITHOUT
reference shutters — the tilt computation needs refs, the x/y range
doesn't. Dropping refs cuts the per-MSA shutter count from ~300 to
100, which (because AssignWcsStep.validate_open_slits is super-linear
in slit count) roughly cubes the per-call speedup. Empirically ~5×
faster overall.

The output is appended to ``data/dispersion_cutoffs.npz`` alongside
the existing keys; pre-existing tilt arrays are preserved.

Usage:
    CRDS_CONTEXT=jwst_1464.pmap python scripts/precompute_shutter_xy.py
    VMPT_XY_ONLY=G140M_F100LP,G395M_F290LP python scripts/precompute_shutter_xy.py
"""

from __future__ import annotations
import logging
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits

# Re-use the existing make_msa_file / patch_rate_file plumbing.
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "scripts"))
import precompute_trace_tilt as ptt  # noqa: E402


OUTPUT = HERE / "vmpt" / "data" / "dispersion_cutoffs.npz"
NQ, NS, ND = 4, 171, 365

# Same 10×10 grid as the tilt precompute so the interpolation keys
# match across the two products.
GRID_ROWS = ptt.GRID_ROWS
GRID_COLS = ptt.GRID_COLS

# Same combo list as the tilt precompute.
COMBOS = ptt.COMBOS


def compute_xy_for_combo(template: Path, disperser: str, filt: str,
                          tmpdir: Path) -> dict[str, np.ndarray]:
    """Per (disperser, filter), open the 10×10 grid in each quadrant,
    run AssignWcsStep on both detectors, record the on-detector x/y
    range of each grid shutter."""
    Ngr = len(GRID_ROWS)
    Ngc = len(GRID_COLS)
    x_lo_nrs1 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_hi_nrs1 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    y_nrs1   = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_lo_nrs2 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_hi_nrs2 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    y_nrs2   = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)

    from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
    from jwst.assign_wcs.nirspec import nrs_wcs_set_input

    n_landed = {"NRS1": 0, "NRS2": 0}
    t_total = {"NRS1": 0.0, "NRS2": 0.0}

    for q in (1, 2, 3, 4):
        # ONLY grid shutters — no reference columns. This is the speed
        # win over precompute_trace_tilt.py.
        shutters = [(q, int(r), int(c)) for r in GRID_ROWS for c in GRID_COLS]
        msa_path = tmpdir / f"msa_xy_{disperser}_{filt}_Q{q}.fits"
        slit_id_to_qrc = ptt.make_msa_file(msa_path, shutters)
        qrc_to_slit_id = {v: k for k, v in slit_id_to_qrc.items()}

        for detector in ("NRS1", "NRS2"):
            rate_path = tmpdir / f"rate_xy_{disperser}_{filt}_{detector}_Q{q}.fits"
            ptt.patch_rate_file(template, rate_path, msa_path,
                                 grating=disperser, filt=filt, detector=detector)
            t0 = time.time()
            try:
                result = AssignWcsStep.call(str(rate_path), save_results=False)
            except Exception as exc:
                if "NoDataOnDetector" in type(exc).__name__:
                    pass
                else:
                    print(f"      ! AssignWcsStep failed on Q{q}/{detector}: "
                          f"{type(exc).__name__}: {exc}")
                continue

            try:
                open_slits = list(
                    result.meta.wcs.get_transform("gwa", "slit_frame").slits
                )
            except Exception:
                open_slits = []
            wstart = float(result.meta.wcsinfo.waverange_start)
            wend = float(result.meta.wcsinfo.waverange_end)
            waves = np.linspace(wstart, wend, 21)
            target_xlo = x_lo_nrs1 if detector == "NRS1" else x_lo_nrs2
            target_xhi = x_hi_nrs1 if detector == "NRS1" else x_hi_nrs2
            target_y   = y_nrs1   if detector == "NRS1" else y_nrs2

            n_set = 0
            for slit in open_slits:
                try:
                    slit_wcs = nrs_wcs_set_input(result, slit.name)
                    s2d = slit_wcs.get_transform("slit_frame", "detector")
                except Exception:
                    continue
                sx0 = float(getattr(slit, "source_xpos", 0.0))
                sy0 = float(getattr(slit, "source_ypos", 0.0))
                try:
                    xd, yd = s2d(np.full_like(waves, sx0),
                                  np.full_like(waves, sy0), waves)
                except Exception:
                    continue
                ok = (np.isfinite(xd) & np.isfinite(yd)
                      & (xd >= 0) & (xd <= 2047)
                      & (yd >= 0) & (yd <= 2047))
                if ok.sum() < 1:
                    continue
                slit_id = int(slit.name)
                if slit_id not in slit_id_to_qrc:
                    continue
                q_g, vmpt_r, vmpt_c = slit_id_to_qrc[slit_id]
                if q_g != q:
                    continue
                try:
                    ri = int(np.where(GRID_ROWS == vmpt_r)[0][0])
                    ci = int(np.where(GRID_COLS == vmpt_c)[0][0])
                except IndexError:
                    continue
                target_xlo[q - 1, ri, ci] = float(xd[ok].min())
                target_xhi[q - 1, ri, ci] = float(xd[ok].max())
                target_y  [q - 1, ri, ci] = float(np.median(yd[ok]))
                n_set += 1

            n_landed[detector] += n_set
            t_total[detector] += time.time() - t0
            result.close()
            print(f"      Q{q}/{detector}: n_set={n_set}/{Ngr*Ngc}, "
                  f"wall={time.time()-t0:.1f}s")

    for d in ("NRS1", "NRS2"):
        print(f"    {d}: {n_landed[d]:3d}/400 grid shutters landed "
              f"({t_total[d]:.1f}s wall)")

    return {"x_lo_nrs1": x_lo_nrs1, "x_hi_nrs1": x_hi_nrs1, "y_nrs1": y_nrs1,
            "x_lo_nrs2": x_lo_nrs2, "x_hi_nrs2": x_hi_nrs2, "y_nrs2": y_nrs2}


def main() -> None:
    os.environ.setdefault("CRDS_CONTEXT", "jwst_1464.pmap")
    template_env = os.environ.get("VMPT_TRACE_TEMPLATE")
    if template_env:
        template = Path(template_env)
    else:
        candidates = list(
            Path("/Users/sunfengwu").glob(
                "jwst_cycle4/emerald_cy4/working/rate*/JWST/"
                "*_nrs1/*_nrs1_rate.fits"
            )
        )
        template = None
        for c in candidates:
            try:
                h = fits.getheader(c, 0)
                if h.get("EXP_TYPE", "") == "NRS_MSASPEC":
                    template = c
                    break
            except Exception:
                pass
        if template is None:
            raise SystemExit(
                "Could not find an NRS_MSASPEC rate template. "
                "Set VMPT_TRACE_TEMPLATE=/path/to/_nrs1_rate.fits."
            )
    print(f"Using template: {template}")
    print(f"Grid: {len(GRID_ROWS)} rows × {len(GRID_COLS)} cols per quadrant")

    for nm in ("jwst", "stpipe", "CRDS", "stpipe.AssignWcsStep"):
        logging.getLogger(nm).setLevel(logging.ERROR)
    warnings.filterwarnings("ignore")

    tmpdir = Path(tempfile.mkdtemp(prefix="vmpt_xy_"))
    print(f"Working dir: {tmpdir}")

    # Output path: by default APPEND to dispersion_cutoffs.npz (single
    # sequential run, current behaviour). If VMPT_XY_OUTPUT is set,
    # write ONLY the new keys to that path — for use in parallel
    # multi-process runs where a separate merger script combines per-
    # process npz files at the end without race conditions.
    out_path_env = os.environ.get("VMPT_XY_OUTPUT", "").strip()
    if out_path_env:
        out_path = Path(out_path_env)
        existing: dict = {}
        print(f"VMPT_XY_OUTPUT={out_path} (writing ONLY new keys here)")
    else:
        if not OUTPUT.exists():
            raise SystemExit(f"Expected {OUTPUT} to exist already.")
        out_path = OUTPUT
        existing = dict(np.load(OUTPUT, allow_pickle=False))
        print(f"Loaded {len(existing)} keys from existing {OUTPUT.name}")

    only = os.environ.get("VMPT_XY_ONLY", "").strip()
    if only:
        wanted = {c.strip().upper() for c in only.split(",")}
        combos_to_run = [(d, f) for (d, f) in COMBOS if f"{d}_{f}" in wanted]
        print(f"VMPT_XY_ONLY={only}  → {len(combos_to_run)} combo(s)")
    else:
        combos_to_run = COMBOS

    for disperser, filt in combos_to_run:
        print(f"\n=== {disperser} {filt} ===")
        t0 = time.time()
        out = compute_xy_for_combo(template, disperser, filt, tmpdir)
        existing[f"{disperser}_{filt}_tilt_grid_rows"] = GRID_ROWS.astype(np.int16)
        existing[f"{disperser}_{filt}_tilt_grid_cols"] = GRID_COLS.astype(np.int16)
        for key in ("x_lo_nrs1", "x_hi_nrs1", "y_nrs1",
                    "x_lo_nrs2", "x_hi_nrs2", "y_nrs2"):
            existing[f"{disperser}_{filt}_{key}"] = out[key]
        np.savez(out_path, **existing)
        print(f"  → saved {len(existing)} arrays "
              f"(combo took {time.time()-t0:.1f}s)")

    print(f"\nFINAL: {len(existing)} arrays in {out_path}")


if __name__ == "__main__":
    main()
