"""Pre-compute the cross-dispersion tilt of NIRSpec MSA spectral
traces, on a 10×10 grid per MSA quadrant, for every disperser/filter
combo vMPT supports.

The runtime overlap calc in `vmpt.main` flags candidate shutters that
share an MSA row with a dispersing source. With the simple model
`s_expected = s_o`, the row-of-the-spectrum is independent of how far
along the spectrum we look — which is only an approximation. In
reality the spectral trace tilts slightly: at the edges of the
spectrum the spectrum's row drifts by up to ±1 MSA row from the
shutter row. The exact tilt depends on the disperser/filter and the
shutter's position in the field.

This script computes, for each (disperser, filter) combo, a tilt
slope `dΔs / dΔv2` (MSA rows per arcsec V2) at 10×10 grid points per
quadrant. The runtime overlap calc bilinearly interpolates this
slope from the grid to the open shutter's actual (row, col), and
augments its row check from `|s_c - s_o| ≤ 1` to
`|s_c - s_o - slope·(v2_c - v2_o)| ≤ 1`.

How the trace is computed:

  1. We synthesise a minimal MSA metadata FITS marking the 400
     grid-point shutters open (10×10 × 4 quadrants).
  2. We clone a real NRS_MSASPEC rate file as a header template,
     patch GRATING/FILTER/MSAMETFL/MSAMETID/DETECTOR, and zero the
     data arrays — assign_wcs only reads the header.
  3. We run `jwst.assign_wcs.AssignWcsStep` once per (combo, detector),
     get each open slit's `slit_frame → detector` transform, sample
     at several wavelengths spanning the disperser's range, and
     translate detector (Δx, Δy) into MSA (Δs, Δv2) using the
     `slit_frame → v2v3` transform from the same WCS pipeline.

This costs ~18 AssignWcsStep calls (9 combos × 2 detectors), each
~10–30s — so ~5–10 minutes wall clock total.

The script appends to the existing `data/dispersion_cutoffs.npz`,
adding these keys per (DISPERSER, FILTER):

  {DISP}_{FILT}_tilt_grid_rows   shape (Ngrid_r,)   vMPT row centres
  {DISP}_{FILT}_tilt_grid_cols   shape (Ngrid_c,)   vMPT col centres
  {DISP}_{FILT}_tilt_slope       shape (4, Ngrid_r, Ngrid_c)  rows/arcsec
  {DISP}_{FILT}_tilt_det         shape (4, Ngrid_r, Ngrid_c)  int detector flag

Usage:
    conda activate stenv
    VMPT_TRACE_TEMPLATE=/path/to/any/_nrs1_rate.fits \
        python scripts/precompute_trace_tilt.py
    # outputs land in data/dispersion_cutoffs.npz next to the existing keys
"""
from __future__ import annotations

import logging
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from astropy.io import fits


# CRDS config — must be set BEFORE importing jwst. Pin to a specific
# context so we don't depend on STScI network reachability at run time.
os.environ.setdefault("CRDS_PATH", str(Path.home() / "crds_cache"))
os.environ.setdefault("CRDS_SERVER_URL", "https://jwst-crds.stsci.edu")
os.environ.setdefault("CRDS_CONTEXT", "jwst_1464.pmap")


HERE = Path(__file__).resolve().parent.parent
OUTPUT = HERE / "vmpt" / "data" / "dispersion_cutoffs.npz"
NQ, NS, ND = 4, 171, 365
MSA_METADATA_ID = 1
DITHER_POINT_INDEX = 1

# Same (disperser, filter) list as the per-shutter cutoffs precompute.
COMBOS = [
    ("PRISM", "CLEAR"),
    ("G140M", "F070LP"),
    ("G140M", "F100LP"),
    ("G235M", "F170LP"),
    ("G395M", "F290LP"),
    ("G140H", "F070LP"),
    ("G140H", "F100LP"),
    ("G235H", "F170LP"),
    ("G395H", "F290LP"),
]

# Detector each MSA quadrant primarily disperses onto (NIRSpec layout):
# Q1 + Q2 → NRS2, Q3 + Q4 → NRS1.  Quadrant codes are 1-based.
QUAD_TO_DETECTOR = {1: "NRS2", 2: "NRS2", 3: "NRS1", 4: "NRS1"}

# 10×10 grid centres per quadrant — uniform 1-based indices in vMPT
# convention (row ∈ 1..171, col ∈ 1..365).
def _grid(n_total: int, n_grid: int) -> np.ndarray:
    """Pick `n_grid` evenly-spaced 1-based indices spanning 1..n_total
    with the first/last samples comfortably away from the edges."""
    return np.linspace(1, n_total, n_grid + 2, dtype=int)[1:-1]


GRID_ROWS = _grid(NS, 10)   # length 10, e.g. [15, 31, ..., 155]
GRID_COLS = _grid(ND, 10)   # length 10, e.g. [33, 66, ..., 332]


def make_msa_file(path: Path, shutters: list[tuple[int, int, int]]) -> dict[int, tuple[int, int, int]]:
    """Build a minimal MSA metadata FITS with the given shutters marked
    open. Returns a slit_id → (quadrant, row_vMPT, col_vMPT) map for
    re-keying results coming back from the WCS pipeline.

    vMPT row ∈ 1..171 (cross-dispersion) → MSA file shutter_column
    vMPT col ∈ 1..365 (dispersion)       → MSA file shutter_row
    """
    n = len(shutters)
    slit_ids = np.arange(1, n + 1, dtype=np.int32)
    quads = np.array([s[0] for s in shutters], dtype=np.int16)
    rows_jwst = np.array([s[2] for s in shutters], dtype=np.int16)  # 1..365
    cols_jwst = np.array([s[1] for s in shutters], dtype=np.int16)  # 1..171
    src_ids = slit_ids.copy()

    shutter_info = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="slitlet_id", format="J", array=slit_ids),
            fits.Column(name="msa_metadata_id", format="I",
                        array=np.full(n, MSA_METADATA_ID, dtype=np.int16)),
            fits.Column(name="shutter_quadrant", format="I", array=quads),
            fits.Column(name="shutter_row", format="I", array=rows_jwst),
            fits.Column(name="shutter_column", format="I", array=cols_jwst),
            fits.Column(name="source_id", format="J", array=src_ids),
            fits.Column(name="background", format="A1",
                        array=np.array(["N"] * n)),
            fits.Column(name="shutter_state", format="A6",
                        array=np.array(["OPEN"] * n)),
            fits.Column(name="estimated_source_in_shutter_x", format="E",
                        array=np.full(n, 0.5, dtype=np.float32)),
            fits.Column(name="estimated_source_in_shutter_y", format="E",
                        array=np.full(n, 0.5, dtype=np.float32)),
            fits.Column(name="dither_point_index", format="I",
                        array=np.full(n, DITHER_POINT_INDEX, dtype=np.int16)),
            fits.Column(name="primary_source", format="A1",
                        array=np.array(["Y"] * n)),
        ],
        name="SHUTTER_INFO",
    )
    source_info = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="source_id", format="J", array=src_ids),
            fits.Column(name="source_name", format="A32",
                        array=np.array([f"SRC{i}" for i in src_ids])),
            fits.Column(name="alias", format="A32",
                        array=np.array([f"S{i}" for i in src_ids])),
            fits.Column(name="ra", format="D",
                        array=np.zeros(n, dtype=np.float64)),
            fits.Column(name="dec", format="D",
                        array=np.zeros(n, dtype=np.float64)),
            fits.Column(name="stellarity", format="E",
                        array=np.zeros(n, dtype=np.float32)),
        ],
        name="SOURCE_INFO",
    )
    img = fits.ImageHDU(data=np.zeros((730, 342), dtype=np.int16),
                        name="SHUTTER_IMAGE")
    fits.HDUList(
        [fits.PrimaryHDU(), shutter_info, source_info, img]
    ).writeto(path, overwrite=True)

    return {int(sid): (int(q), int(c_jwst), int(r_jwst))
            for sid, q, r_jwst, c_jwst in zip(
                slit_ids, quads, rows_jwst, cols_jwst)}


def patch_rate_file(template: Path, out: Path, msa_path: Path, *,
                    grating: str, filt: str, detector: str) -> None:
    """Clone a real NRS_MSASPEC rate file, patch headers, zero data."""
    with fits.open(template, memmap=False) as hdul:
        h = hdul[0].header
        h["MSAMETFL"] = str(msa_path)
        h["MSAMETID"] = MSA_METADATA_ID
        h["PATT_NUM"] = DITHER_POINT_INDEX
        h["NUMDTHPT"] = 1
        h["GRATING"] = grating
        h["FILTER"] = filt
        h["DETECTOR"] = detector
        for ext in hdul[1:]:
            if ext.data is not None and hasattr(ext.data, "shape"):
                if ext.data.dtype.kind == "f":
                    ext.data[:] = 0.0
                elif ext.data.dtype.kind in ("i", "u"):
                    ext.data[:] = 0
        hdul.writeto(out, overwrite=True)


def compute_slopes_for_combo(template: Path, disperser: str, filt: str,
                              tmpdir: Path) -> dict[str, np.ndarray]:
    """For one (disperser, filter), run AssignWcsStep on both detectors
    and return arrays sized (4, len(GRID_ROWS), len(GRID_COLS)).

    The slope captures the COMBINED tilt that APT MPT shows in its
    spec-overlap stripe: the answer to "at column d_B (i.e. V2 offset
    Δv2 from the open shutter at d_A), what MSA s-row does the open's
    spectrum land on?" This is the field-tilt-plus-spectrum-curvature
    that drives the discrete row jumps in MPT's overlap rendering.

    To measure it we trace the open shutter's spectrum across many
    wavelengths AND also trace a "reference" shutter at the same s
    but a different d, then compute

        Δy_combined = y_open_spectrum_at_x(d_ref)  −  y_ref_at_central_λ
        Δs_combined = −Δy_combined / PX_PER_MSA_ROW
        slope       =  Δs_combined / Δv2(d_open, d_ref)

    Decomposed, Δy_combined is the sum of (a) the single-shutter
    spectrum drift between λ_central and the wavelength that puts the
    spectrum at d_ref's x, and (b) the small field-tilt between two
    same-row shutters at different d. Earlier versions of this script
    only used (a), missing ~30–40 % of the magnitude (the user verified
    against APT MPT: q1 d247 s93 → d153 should snap to −1 row, but the
    spectrum-only slope rounded to 0).
    """
    Ngr = len(GRID_ROWS)
    Ngc = len(GRID_COLS)
    slope = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    # Diagnostics: detector-pixel Δy at the reference column, and the
    # V2 distance used. Useful for sanity-checking the slope against
    # MPT screenshots after the fact.
    dy_combined = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    dv2_used = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    det = np.zeros((NQ, Ngr, Ngc), dtype=np.int8)

    # Per-shutter on-detector x/y range, per detector. NaN means the
    # spectrum doesn't reach that detector for the sampled grid
    # shutter. At runtime, bilinear interpolation across the 10×10
    # quadrant grid gives a per-shutter (x_lo, x_hi, y) estimate that
    # the spec-overlap check uses to find actual detector-pixel
    # collisions between any two shutters — not just same-quadrant
    # ones.
    x_lo_nrs1 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_hi_nrs1 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    y_nrs1   = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_lo_nrs2 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    x_hi_nrs2 = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)
    y_nrs2   = np.full((NQ, Ngr, Ngc), np.nan, dtype=np.float32)

    # NIRSpec MSA cross-dispersion pitch ≈ 5.0 detector pixels per MSA
    # s-row — measured empirically from adjacent same-d shutters
    # (e.g. q1 d247 rows 88–98 had y_det pitch of 4.99 px/row).
    PX_PER_MSA_ROW = 5.0

    # Column offsets for the reference shutters. PRISM's full
    # spectrum spans only ~160 cols on the detector (asymmetric,
    # roughly −50 to +80 from the shutter), so ±100 puts the
    # reference outside the spectrum's x-range for most grid points.
    # ±60 (≈ 12″ V2) bracketing works for almost every grid+ref pair:
    # both refs land on the same detector with their central x inside
    # the spectrum trace, so the slope averaging captures both blue-
    # and red-side drift. The slope is local; non-linearity at larger
    # Δv2 is < 10 % for PRISM and smaller for the gratings (where
    # v2_overlap is much larger than the spectrum's on-detector
    # extent), so a single local slope reproduces MPT's snap.
    REF_D_OFFSETS = (-60, +60)

    # Wavelength range comes from CRDS wavelengthrange model.
    from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
    from jwst.assign_wcs.nirspec import nrs_wcs_set_input

    # Cache MSA V2 from the shipped npz — used to compute Δv2 between
    # the open grid shutter and its reference.
    msa_v2_path = HERE / "vmpt" / "data" / "nirspec_msa_v2v3.npz"
    msa = np.load(msa_v2_path)
    V2_MSA = msa["v2_msa"]  # (4, 171, 365) absolute V2 arcsec

    n_landed_total = {"NRS1": 0, "NRS2": 0}
    t_total = {"NRS1": 0.0, "NRS2": 0.0}
    for q in (1, 2, 3, 4):
        # Each grid shutter ALSO gets up to two references at the same
        # s but at d ± REF_D_OFFSETS — bracketing the open shutter
        # so the slope captures the average drift across a baseline
        # comparable to v2_overlap. Edge grid columns may have only
        # one valid reference (the other would fall off the MSA); we
        # accept that and use whatever is available.
        shutters: list[tuple[int, int, int]] = []
        refs_for_grid: dict[
            tuple[int, int, int], list[tuple[int, int, int]]
        ] = {}
        for r in GRID_ROWS:
            for c in GRID_COLS:
                grid_key = (q, int(r), int(c))
                shutters.append(grid_key)
                refs: list[tuple[int, int, int]] = []
                for off in REF_D_OFFSETS:
                    d_ref = int(c) + off
                    if 1 <= d_ref <= ND:
                        refs.append((q, int(r), d_ref))
                if refs:
                    refs_for_grid[grid_key] = refs
                    shutters.extend(refs)
        # Deduplicate while preserving order.
        seen = set()
        unique_shutters: list[tuple[int, int, int]] = []
        for sh in shutters:
            if sh not in seen:
                seen.add(sh)
                unique_shutters.append(sh)
        shutters = unique_shutters

        msa_path = tmpdir / f"msa_{disperser}_{filt}_Q{q}.fits"
        slit_id_to_qrc = make_msa_file(msa_path, shutters)
        qrc_to_slit_id = {v: k for k, v in slit_id_to_qrc.items()}

        for detector in ("NRS1", "NRS2"):
            det_flag = 1 if detector == "NRS1" else 2
            rate_path = tmpdir / f"rate_{disperser}_{filt}_{detector}_Q{q}.fits"
            patch_rate_file(template, rate_path, msa_path,
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
            # Sample more wavelengths so interpolation at the reference
            # column's x is accurate.
            waves = np.linspace(wstart, wend, 15)
            wave_central = 0.5 * (wstart + wend)

            # Cache: (slit_name) → (x_arr, y_arr, x_at_central, y_at_central)
            slit_data: dict[int, tuple[np.ndarray, np.ndarray, float, float]] = {}
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
                    x_c, y_c = s2d(np.array([sx0]), np.array([sy0]),
                                    np.array([wave_central]))
                except Exception:
                    continue
                if not np.any(np.isfinite(xd)):
                    continue
                slit_data[int(slit.name)] = (
                    np.asarray(xd), np.asarray(yd),
                    float(x_c[0]), float(y_c[0]),
                )

            # Record the spectrum's on-detector x range and median y
            # for every grid shutter that landed on THIS detector.
            # Only samples within the detector's [0, 2047] bounds
            # count — points outside are clipped (e.g. the spectrum
            # extending beyond x=2047 wraps to the other detector).
            xy_target = (x_lo_nrs1, x_hi_nrs1, y_nrs1) if detector == "NRS1" \
                else (x_lo_nrs2, x_hi_nrs2, y_nrs2)
            for grid_key in {(q, int(r), int(c))
                              for r in GRID_ROWS for c in GRID_COLS}:
                _q_g, vmpt_r, vmpt_c = grid_key
                grid_id = qrc_to_slit_id.get(grid_key)
                if grid_id is None or grid_id not in slit_data:
                    continue
                xd_g, yd_g, _xc_g, _yc_g = slit_data[grid_id]
                ok = (np.isfinite(xd_g) & np.isfinite(yd_g)
                      & (xd_g >= 0) & (xd_g <= 2047)
                      & (yd_g >= 0) & (yd_g <= 2047))
                if ok.sum() < 1:
                    continue
                try:
                    ri = int(np.where(GRID_ROWS == vmpt_r)[0][0])
                    ci = int(np.where(GRID_COLS == vmpt_c)[0][0])
                except IndexError:
                    continue
                xy_target[0][q - 1, ri, ci] = float(xd_g[ok].min())
                xy_target[1][q - 1, ri, ci] = float(xd_g[ok].max())
                xy_target[2][q - 1, ri, ci] = float(np.median(yd_g[ok]))

            # Iterate over grid shutters and combine with their
            # reference shutters. We average across all available
            # references (one per side) so the slope captures the
            # average drift over a baseline comparable to v2_overlap.
            n_set_this = 0
            for grid_key, ref_keys in refs_for_grid.items():
                _q_g, vmpt_r, vmpt_c = grid_key
                grid_id = qrc_to_slit_id.get(grid_key)
                if grid_id is None or grid_id not in slit_data:
                    continue
                xd_g, yd_g, _xc_g, _yc_g = slit_data[grid_id]
                ok = np.isfinite(xd_g) & np.isfinite(yd_g)
                if ok.sum() < 4:
                    continue
                x_g = xd_g[ok]
                y_g = yd_g[ok]
                order = np.argsort(x_g)
                x_sorted = x_g[order]
                y_sorted = y_g[order]
                v2_open = float(V2_MSA[q - 1, vmpt_r - 1, vmpt_c - 1])

                slopes_this: list[float] = []
                dy_samples: list[float] = []
                dv2_samples: list[float] = []
                for ref_key in ref_keys:
                    _q_r, _r_r, vmpt_c_ref = ref_key
                    ref_id = qrc_to_slit_id.get(ref_key)
                    if ref_id is None or ref_id not in slit_data:
                        continue
                    _xd_r, _yd_r, xc_r, yc_r = slit_data[ref_id]
                    if not (np.isfinite(xc_r) and np.isfinite(yc_r)):
                        continue
                    if not (x_sorted[0] <= xc_r <= x_sorted[-1]):
                        # Reference column's x is outside the open
                        # shutter's spectrum trace on this detector.
                        continue
                    y_at_xc_r = float(np.interp(xc_r, x_sorted, y_sorted))
                    dy_comb = y_at_xc_r - yc_r           # det px
                    ds_comb = -dy_comb / PX_PER_MSA_ROW  # MSA rows
                    v2_ref = float(V2_MSA[q - 1, vmpt_r - 1, vmpt_c_ref - 1])
                    dv2 = v2_ref - v2_open
                    if abs(dv2) < 0.1:
                        continue
                    slopes_this.append(ds_comb / dv2)
                    dy_samples.append(dy_comb)
                    dv2_samples.append(dv2)

                if not slopes_this:
                    continue
                k_slope = float(np.mean(slopes_this))
                dy_avg = float(np.mean(dy_samples))
                dv2_avg = float(np.mean(np.abs(dv2_samples)))

                try:
                    ri = int(np.where(GRID_ROWS == vmpt_r)[0][0])
                    ci = int(np.where(GRID_COLS == vmpt_c)[0][0])
                except IndexError:
                    continue

                if det[q - 1, ri, ci] == 0:
                    slope[q - 1, ri, ci] = k_slope
                    dy_combined[q - 1, ri, ci] = dy_avg
                    dv2_used[q - 1, ri, ci] = dv2_avg
                    det[q - 1, ri, ci] = det_flag
                    n_landed_total[detector] += 1
                    n_set_this += 1
            print(f"      [{detector}] Q{q}: slit_data={len(slit_data)}, "
                  f"n_set={n_set_this}, "
                  f"first non-NaN slope on this pass: "
                  f"{np.nanmedian(slope[q-1])}")

            elapsed = time.time() - t0
            t_total[detector] += elapsed
            result.close()

    for d in ("NRS1", "NRS2"):
        print(f"    {d}: {n_landed_total[d]:3d}/400 slits landed "
              f"({t_total[d]:.1f}s wall)")

    return {"tilt_slope": slope, "tilt_dy_combined": dy_combined,
            "tilt_dv2": dv2_used, "tilt_det": det.astype(np.int8),
            "x_lo_nrs1": x_lo_nrs1, "x_hi_nrs1": x_hi_nrs1, "y_nrs1": y_nrs1,
            "x_lo_nrs2": x_lo_nrs2, "x_hi_nrs2": x_hi_nrs2, "y_nrs2": y_nrs2}


def main() -> None:
    template_env = os.environ.get("VMPT_TRACE_TEMPLATE")
    if template_env:
        template = Path(template_env)
    else:
        # Look for any local NRS_MSASPEC rate file
        candidates = list(
            Path("/Users/sunfengwu").glob(
                "jwst_cycle4/emerald_cy4/working/rate*/JWST/"
                "*_nrs1/*_nrs1_rate.fits"
            )
        )
        # Pick the first that's NRS_MSASPEC.
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
    print(f"Grid rows (vMPT 1..171): {list(GRID_ROWS)}")
    print(f"Grid cols (vMPT 1..365): {list(GRID_COLS)}")

    # Quiet the JWST pipeline noise.
    for nm in ("jwst", "stpipe", "CRDS", "stpipe.AssignWcsStep"):
        logging.getLogger(nm).setLevel(logging.ERROR)

    tmpdir = Path(tempfile.mkdtemp(prefix="vmpt_trace_tilt_"))
    print(f"Working dir: {tmpdir}")

    # Load existing npz to preserve all existing keys.
    if not OUTPUT.exists():
        raise SystemExit(f"Expected {OUTPUT} to exist already.")
    existing = dict(np.load(OUTPUT, allow_pickle=False))
    print(f"Loaded {len(existing)} keys from existing {OUTPUT.name}")

    # Optional env-var filter for partial reruns: e.g.
    #   VMPT_TRACE_ONLY=G395M_F290LP python scripts/precompute_trace_tilt.py
    # processes ONLY the matching combo (comma-separated for multiple).
    only = os.environ.get("VMPT_TRACE_ONLY", "").strip()
    if only:
        wanted = {c.strip().upper() for c in only.split(",")}
        combos_to_run = [(d, f) for (d, f) in COMBOS if f"{d}_{f}" in wanted]
        print(f"VMPT_TRACE_ONLY={only}  → {len(combos_to_run)} combo(s)")
    else:
        combos_to_run = COMBOS

    for disperser, filt in combos_to_run:
        print(f"\n=== {disperser} {filt} ===")
        t0 = time.time()
        out = compute_slopes_for_combo(template, disperser, filt, tmpdir)
        existing[f"{disperser}_{filt}_tilt_grid_rows"] = GRID_ROWS.astype(np.int16)
        existing[f"{disperser}_{filt}_tilt_grid_cols"] = GRID_COLS.astype(np.int16)
        existing[f"{disperser}_{filt}_tilt_slope"] = out["tilt_slope"]
        existing[f"{disperser}_{filt}_tilt_dy_combined"] = out["tilt_dy_combined"]
        existing[f"{disperser}_{filt}_tilt_dv2"] = out["tilt_dv2"]
        existing[f"{disperser}_{filt}_tilt_det"] = out["tilt_det"]
        # Per-shutter on-detector x/y range (one (10×10) grid per
        # quadrant, per detector). The live overlay bilinearly
        # interpolates these to any (s, d) at runtime and checks
        # x-range intersection on each detector for the spec-overlap
        # rule.
        existing[f"{disperser}_{filt}_x_lo_nrs1"] = out["x_lo_nrs1"]
        existing[f"{disperser}_{filt}_x_hi_nrs1"] = out["x_hi_nrs1"]
        existing[f"{disperser}_{filt}_y_nrs1"]    = out["y_nrs1"]
        existing[f"{disperser}_{filt}_x_lo_nrs2"] = out["x_lo_nrs2"]
        existing[f"{disperser}_{filt}_x_hi_nrs2"] = out["x_hi_nrs2"]
        existing[f"{disperser}_{filt}_y_nrs2"]    = out["y_nrs2"]
        # Strip old keys from the previous (single-shutter) formula so
        # they don't get carried forward and confuse downstream code.
        for old_key in (f"{disperser}_{filt}_tilt_dy_total",
                         f"{disperser}_{filt}_tilt_dx_total"):
            existing.pop(old_key, None)
        # Incremental save: write the npz after EVERY combo so a
        # killed run still leaves usable data for the combos that
        # already finished. The wall-clock cost of np.savez is small
        # (~1 s) compared to the per-combo AssignWcsStep cost
        # (10–40 min), so we do it eagerly.
        np.savez(OUTPUT, **existing)
        print(f"  → saved {len(existing)} arrays to "
              f"{OUTPUT.name} (cumulative)")
        # Quick stats
        vals = out["tilt_slope"][np.isfinite(out["tilt_slope"])]
        if vals.size:
            print(f"  slope (rows/arcsec): "
                  f"median={np.median(vals):+.4f} "
                  f"range=[{vals.min():+.4f}, {vals.max():+.4f}] "
                  f"abs_p90={np.percentile(np.abs(vals), 90):.4f}  "
                  f"n={vals.size}")
        else:
            print("  no slopes computed for this combo")
        print(f"  combo took {time.time()-t0:.1f}s")

    np.savez(OUTPUT, **existing)
    print(f"\nSaved {len(existing)} arrays to {OUTPUT}")


if __name__ == "__main__":
    main()
