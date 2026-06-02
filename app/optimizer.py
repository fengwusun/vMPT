"""MSA pointing optimizer.

Searches for (RA, Dec, V3 PA) maximising the number — or weighted flux —
of catalog sources that fall in operable, well-centred MSA shutters.

The algorithm is derived from **hMPT** by Daniel Eisenstein, Samuel
McCarty, and Zihao Wu (Harvard / CfA), itself a Python re-implementation
of ESA's eMPT (Bonaventura et al. 2023, A&A 672, A40):
<https://github.com/zihaowu-astro/hMPT>.

vMPT does NOT depend on hMPT — this module re-implements the core
algorithm in vMPT's style with attribution, so it composes cleanly
with our existing MSA grid (`data/nirspec_msa_v2v3.npz`), CRDS
operability loader, and Bokeh UI.

Algorithm summary
-----------------
1. **`radec_to_axy`** — vectorised gnomonic projection of source
   (RA, Dec) onto the MSA aperture plane (ax, ay), with optional
   differential-velocity-aberration scaling and the PA rotation.
2. **`axy_to_shutter`** — per-quadrant CloughTocher2D interpolation
   maps (ax, ay) → fractional shutter indices (quad, s_row, d_col).
   Built lazily from the shutter centres vMPT already loads.
3. **`PointingEvaluator.evaluate`** — combines the above with the
   operability mask (incl. a 3-shutter vertical slit constraint),
   a configurable APT-style centration buffer, and a Gaussian-PSF
   throughput fraction.
4. **`grid_search`** — brute-force ranking over a (ΔRA, ΔDec, ΔPA)
   cube.
5. **`refine_top`** — `scipy.optimize.differential_evolution` polish
   of the top-N grid candidates inside a small box.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.optimize import differential_evolution
from scipy.special import erfc

from .coords import MSA_V2_REF, MSA_V3_REF, V3_IDL_Y_ANGLE
from .msa import load_msa_grid, load_operability


# Physical shutter dimensions on the focal plane (arcsec).
# `SHUTTER_X` is the dispersion direction (columns, 0..364); `SHUTTER_Y`
# is the spatial direction (rows, 0..170). Values from hMPT, which
# matches APT's MSA model.
SHUTTER_X_ARCSEC = 0.2679
SHUTTER_Y_ARCSEC = 0.5294

# Centration buffer classes (inset from shutter edge, arcsec).
# Mirrors APT's source-centering modes; values from hMPT.
CENTRATION_BUFFERS = {
    "UNCONSTRAINED":       0.000,
    "ENTIRE_OPEN":         0.035,
    "MIDPOINT":            0.059,
    "CONSTRAINED":         0.072,
    "TIGHTLY_CONSTRAINED": 0.091,
}

# The MSA frame rotates from V2/V3 into the aperture (ax, ay) frame by
# angle Φ. pysiaf reports V3IdlYAngle for NRS_FULL_MSA ≈ 138.575°;
# hMPT writes this as PHI = 41.42543 with the convention PHI = 180 − V3IdlYAngle.
# We use V3_IDL_Y_ANGLE from coords.py as the source of truth.
_ROT_AXY_DEG: float = 180.0 - V3_IDL_Y_ANGLE


# Lazy caches.
_inverse_cache: list[dict] | None = None


# ---------------------------------------------------------------------
# Coordinate maths
# ---------------------------------------------------------------------


def _rotation_matrix(theta_rad: float) -> np.ndarray:
    """2×2 rotation by `theta_rad`. Conventions match hMPT — applied
    via right-multiplication: ``axy = v23 @ R``."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, s], [-s, c]])


def _build_inverse_interpolators() -> list[dict]:
    """Build per-quadrant Axy→(s, d) interpolators from the MSA grid.

    Cached after the first call. Construction takes ~1–3 s on a laptop
    (Delaunay triangulation over ~62k points per quadrant), so this is
    deferred to first lookup rather than import-time.
    """
    global _inverse_cache
    if _inverse_cache is not None:
        return _inverse_cache

    v2, v3 = load_msa_grid()                       # (4, 171, 365) each
    rot = _rotation_matrix(np.deg2rad(_ROT_AXY_DEG))

    interpolators: list[dict] = []
    for q in range(4):
        dv2 = v2[q] - MSA_V2_REF
        dv3 = v3[q] - MSA_V3_REF
        v23 = np.stack([dv2, dv3], axis=-1)        # (171, 365, 2)
        axy = v23 @ rot                            # (171, 365, 2)

        n_s, n_d, _ = axy.shape
        ss, dd = np.meshgrid(np.arange(n_s), np.arange(n_d), indexing="ij")
        points = axy.reshape(-1, 2)
        s_vals = ss.ravel().astype(float)
        d_vals = dd.ravel().astype(float)

        interp_s = CloughTocher2DInterpolator(points, s_vals)
        interp_d = CloughTocher2DInterpolator(points, d_vals)
        ax_lo, ax_hi = float(points[:, 0].min()), float(points[:, 0].max())
        ay_lo, ay_hi = float(points[:, 1].min()), float(points[:, 1].max())
        interpolators.append({
            "interp_s": interp_s, "interp_d": interp_d,
            "ax_bounds": (ax_lo, ax_hi),
            "ay_bounds": (ay_lo, ay_hi),
        })

    _inverse_cache = interpolators
    return _inverse_cache


def radec_to_axy(
    ra: np.ndarray,
    dec: np.ndarray,
    ra_p: float,
    dec_p: float,
    pa_v3_deg: float,
    theta_deg: float = 90.0,
) -> np.ndarray:
    """Project (RA, Dec) onto MSA aperture coords (ax, ay) in arcsec.

    ``theta_deg`` is the APT differential-velocity-aberration parameter
    (date-dependent — exported from APT's XML). The default 90° is the
    no-correction case used by hMPT during planning, which agrees with
    APT to ≲ 1 mas at typical pointings.
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    dra = np.deg2rad(ra - ra_p)
    dec_r = np.deg2rad(dec)
    dec_pr = np.deg2rad(dec_p)
    denom = (np.sin(dec_r) * np.sin(dec_pr)
             + np.cos(dec_r) * np.cos(dec_pr) * np.cos(dra))
    denom_arcsec = denom * np.pi / 3600.0 / 180.0
    # Small-angle gnomonic projection (west→east, south→north in arcsec).
    x = np.cos(dec_r) * np.sin(dra) / denom_arcsec
    y = ((np.sin(dec_r) * np.cos(dec_pr)
          - np.cos(dec_r) * np.sin(dec_pr) * np.cos(dra))
         / denom_arcsec)
    # Differential velocity aberration (small magnification).
    m_dva = 1.0 / (1.0 - 30.0 / 3e5 * np.cos(np.deg2rad(theta_deg - pa_v3_deg)))
    x *= m_dva
    y *= m_dva
    # PA rotation into V2/V3.
    th = np.deg2rad(pa_v3_deg)
    v2 = np.cos(th) * x - np.sin(th) * y
    v3 = np.sin(th) * x + np.cos(th) * y
    # V2/V3 → aperture (ax, ay).
    rot = _rotation_matrix(np.deg2rad(_ROT_AXY_DEG))
    v23 = np.stack([v2, v3], axis=-1)
    return v23 @ rot


def axy_to_shutter(
    axy: np.ndarray,
    interpolators: Optional[list[dict]] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(quad, s_frac, d_frac)`` per source.

    ``quad`` is 1–4 for sources inside a quadrant, 0 for outside.
    ``s_frac`` and ``d_frac`` are fractional shutter indices
    (``s ∈ [0, 170]``, ``d ∈ [0, 364]``). NaN where ``quad == 0``.

    Vignetting cutoffs at each quadrant's inner corner mirror hMPT's
    `find_shutter_from_Axy` (lines 437–445 of msa_planner.py).
    """
    if interpolators is None:
        interpolators = _build_inverse_interpolators()
    axy = np.atleast_2d(axy)
    n = axy.shape[0]
    quad = np.zeros(n, dtype=int)
    s_frac = np.full(n, np.nan)
    d_frac = np.full(n, np.nan)

    for q in range(4):
        m = interpolators[q]
        ax_lo, ax_hi = m["ax_bounds"]
        ay_lo, ay_hi = m["ay_bounds"]
        in_box = ((axy[:, 0] >= ax_lo) & (axy[:, 0] <= ax_hi) &
                  (axy[:, 1] >= ay_lo) & (axy[:, 1] <= ay_hi))
        if not in_box.any():
            continue
        pr = m["interp_s"](axy[in_box])
        pc = m["interp_d"](axy[in_box])
        valid = ((pr >= -0.5) & (pr <= 170.5) &
                 (pc >= -0.5) & (pc <= 364.5) &
                 ~np.isnan(pr) & ~np.isnan(pc))
        # Inner-corner vignetting per quadrant (hMPT values).
        if q == 0:
            valid &= (pr >= 11.5) & (pc >= 8.5)
        elif q == 1:
            valid &= (pr <= 158.5) & (pc >= 8.5)
        elif q == 2:
            valid &= (pr >= 11.5) & (pc <= 356.5)
        elif q == 3:
            valid &= (pr <= 157.5) & (pc <= 359.5)
        # Stash results at the original positions.
        good_idx = np.where(in_box)[0][valid]
        s_frac[good_idx] = pr[valid]
        d_frac[good_idx] = pc[valid]
        quad[good_idx] = q + 1
    return quad, s_frac, d_frac


# ---------------------------------------------------------------------
# PSF / centration helpers
# ---------------------------------------------------------------------


def _integrate_gaussian(mean: np.ndarray, sigma: float,
                        lo: float, hi: float) -> np.ndarray:
    """∫_lo^hi 𝒩(mean, σ²) dx."""
    lo_z = (lo - mean) / sigma
    hi_z = (hi - mean) / sigma
    return 0.5 * (erfc(lo_z * np.sqrt(0.5)) - erfc(hi_z * np.sqrt(0.5)))


def _gaussian_through_shutter(s_frac: np.ndarray, d_frac: np.ndarray,
                              sigma_arcsec: float) -> np.ndarray:
    """Fraction of a circular Gaussian PSF (σ arcsec) transmitted by
    a single shutter at the given fractional (s, d) position within
    the shutter."""
    off_y = (s_frac - np.rint(s_frac)) * SHUTTER_Y_ARCSEC
    off_x = (d_frac - np.rint(d_frac)) * SHUTTER_X_ARCSEC
    # Shutter clear aperture: ~±0.23″ vertically, ~±0.10″ horizontally.
    half_y, half_x = 0.23, 0.10
    return (_integrate_gaussian(off_y, sigma_arcsec, -half_y, half_y)
            * _integrate_gaussian(off_x, sigma_arcsec, -half_x, half_x))


# ---------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------


class PointingEvaluator:
    """One catalog × one MSA = a re-usable per-pointing scorer.

    Caches the interpolators and operability mask so repeated
    ``evaluate(ra, dec, pa)`` calls are fast — the grid search runs
    this hundreds-of-thousands of times.

    Parameters
    ----------
    ra_sources, dec_sources : array-like
        Source positions in degrees.
    flux_sources : array-like, optional
        Source fluxes (linear units). Used for the ``"flux"`` objective.
    sigma_arcsec : float
        Gaussian PSF σ for the throughput integration.
    centration : str
        One of the keys in ``CENTRATION_BUFFERS``.
    slit_length : int
        Vertical extent of the slitlet (1, 2, 3 or 5 shutters); every
        shutter in the slitlet must be operable for the source to count.
    operable : ndarray, optional
        Pre-loaded (4, 171, 365) operability mask. Loaded lazily if None.
    """

    def __init__(
        self,
        ra_sources,
        dec_sources,
        flux_sources=None,
        sigma_arcsec: float = 0.06,
        centration: str = "UNCONSTRAINED",
        slit_length: int = 3,
        operable: Optional[np.ndarray] = None,
    ):
        self.ra = np.asarray(ra_sources, dtype=float)
        self.dec = np.asarray(dec_sources, dtype=float)
        self.flux = (np.ones_like(self.ra) if flux_sources is None
                     else np.asarray(flux_sources, dtype=float))
        self.sigma = float(sigma_arcsec)
        self.buffer = CENTRATION_BUFFERS.get(
            centration.upper(), CENTRATION_BUFFERS["UNCONSTRAINED"])
        self.slit_length = int(slit_length)
        if operable is None:
            operable, _ = load_operability()
        self.operable = np.asarray(operable, dtype=bool)
        self.interpolators = _build_inverse_interpolators()

    def evaluate(
        self,
        ra_p: float, dec_p: float, pa_v3: float,
        theta_deg: float = 90.0,
    ) -> tuple[np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Return ``(detected_bool, throughput, (quad, s, d))`` per source."""
        axy = radec_to_axy(self.ra, self.dec, ra_p, dec_p, pa_v3, theta_deg)
        quad, s_frac, d_frac = axy_to_shutter(axy, self.interpolators)
        operable_mask = self._check_operable(quad, s_frac, d_frac)
        centered = self._check_centration(s_frac, d_frac)
        with np.errstate(invalid="ignore"):
            tp = _gaussian_through_shutter(s_frac, d_frac, self.sigma)
        tp = np.where(operable_mask & centered, tp, 0.0)
        detected = tp > 0
        return detected, tp, (quad, s_frac, d_frac)

    # -- internals ---------------------------------------------------

    def _check_operable(
        self, quad: np.ndarray, s_frac: np.ndarray, d_frac: np.ndarray,
    ) -> np.ndarray:
        """All `slit_length` consecutive rows centred on the source must
        be operable. Off-grid sources (``quad == 0``) fail by default."""
        out = np.zeros(len(quad), dtype=bool)
        in_grid = quad > 0
        if not in_grid.any():
            return out
        idx = np.where(in_grid)[0]
        q0 = quad[idx] - 1
        s0 = np.rint(s_frac[idx]).astype(int)
        d0 = np.rint(d_frac[idx]).astype(int)

        valid = np.ones(len(idx), dtype=bool)
        half = self.slit_length // 2
        for ds in range(-half, half + 1):
            s_off = s0 + ds
            in_range = ((s_off >= 0) & (s_off < 171)
                        & (d0 >= 0) & (d0 < 365))
            this_ok = np.zeros(len(idx), dtype=bool)
            if in_range.any():
                ir = np.where(in_range)[0]
                this_ok[ir] = self.operable[q0[ir], s_off[ir], d0[ir]]
            valid &= this_ok
        out[idx] = valid
        return out

    def _check_centration(
        self, s_frac: np.ndarray, d_frac: np.ndarray,
    ) -> np.ndarray:
        """Source falls within the centration buffer (in BOTH axes)."""
        with np.errstate(invalid="ignore"):
            row_limit = 0.5 - (self.buffer / SHUTTER_Y_ARCSEC)
            col_limit = 0.5 - (self.buffer / SHUTTER_X_ARCSEC)
            off_r = np.abs(s_frac - np.rint(s_frac))
            off_c = np.abs(d_frac - np.rint(d_frac))
            return ((off_r < row_limit) & (off_c < col_limit)
                    & np.isfinite(s_frac) & np.isfinite(d_frac))


# ---------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------


def grid_search(
    evaluator: PointingEvaluator,
    ra0: float, dec0: float, pa0: float,
    *,
    dra_arcsec: float = 30.0,
    ddec_arcsec: float = 30.0,
    dpa_deg: float = 30.0,
    n_ra: int = 20,
    n_dec: int = 20,
    n_pa: int = 20,
    weights: Optional[np.ndarray] = None,
    objective: str = "number",
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> dict:
    """Brute-force ranking over a (ΔRA, ΔDec, ΔPA) cube.

    The ΔRA / ΔDec arguments are in *arcseconds*; the ΔRA span is
    automatically scaled by 1/cos(Dec) so the box is roughly square on
    the sky. ``progress_cb(done, total)`` is invoked at ~2 % increments
    so the UI can report progress.

    If any of ``dra_arcsec``, ``ddec_arcsec``, ``dpa_deg`` is ≤ 0, that
    axis is FROZEN at the central value (``n`` is forced to 1, no
    sweep). This is the convention the UI uses to mean "keep this
    coordinate at its current value."
    """
    # Freeze axes whose delta is zero or negative — corresponds to
    # "do not search this dimension".
    if dra_arcsec <= 0:
        n_ra = 1
    if ddec_arcsec <= 0:
        n_dec = 1
    if dpa_deg <= 0:
        n_pa = 1

    cos_dec = max(np.cos(np.deg2rad(dec0)), 1e-3)
    dra_deg = dra_arcsec / 3600.0 / cos_dec
    ddec_deg = ddec_arcsec / 3600.0
    # `linspace(0, 0, 1)` returns [0.0] — exactly the centre, which is
    # what "frozen" should produce.
    ras = ra0 + (np.array([0.0]) if n_ra == 1
                 else np.linspace(-dra_deg, dra_deg, n_ra))
    decs = dec0 + (np.array([0.0]) if n_dec == 1
                   else np.linspace(-ddec_deg, ddec_deg, n_dec))
    pas = pa0 + (np.array([0.0]) if n_pa == 1
                 else np.linspace(-dpa_deg, dpa_deg, n_pa))
    if weights is None:
        weights = np.ones_like(evaluator.ra)
    weights = np.asarray(weights, dtype=float)

    n_total = n_ra * n_dec * n_pa
    scores = np.empty(n_total, dtype=float)
    ras_out = np.empty(n_total, dtype=float)
    decs_out = np.empty(n_total, dtype=float)
    pas_out = np.empty(n_total, dtype=float)

    report_every = max(1, n_total // 50)
    use_flux = (objective == "flux")
    idx = 0
    for ra in ras:
        for dec in decs:
            for pa in pas:
                det, tp, _ = evaluator.evaluate(ra, dec, pa)
                if use_flux:
                    s = float(np.sum(tp * evaluator.flux * weights))
                else:
                    s = float(np.sum(det * weights))
                scores[idx] = s
                ras_out[idx] = ra
                decs_out[idx] = dec
                pas_out[idx] = pa
                idx += 1
                if progress_cb is not None and (idx % report_every == 0):
                    progress_cb(idx, n_total)
    if progress_cb is not None:
        progress_cb(n_total, n_total)

    order = np.argsort(-scores)
    return {
        "score": scores[order],
        "ra": ras_out[order],
        "dec": decs_out[order],
        "pa": pas_out[order],
    }


def refine_top(
    evaluator: PointingEvaluator,
    grid_results: dict,
    *,
    n_top: int = 10,
    dra_arcsec: float = 2.0,
    ddec_arcsec: float = 2.0,
    dpa_deg: float = 2.0,
    maxiter: int = 200,
    weights: Optional[np.ndarray] = None,
    objective: str = "number",
    progress_cb: Optional[Callable[[int, int], None]] = None,
    dedup_tol: tuple[float, float, float] = (0.3, 0.3, 0.05),
) -> dict:
    """Differential-evolution polish of the top-N grid candidates.

    Each candidate is refined inside a small (dra, ddec, dpa) box.
    Returns a fresh ranked dict in the same schema as `grid_search`.

    ``dedup_tol`` is ``(arcsec_ra, arcsec_dec, deg_pa)``: refined
    solutions within these tolerances of an earlier (higher-scoring)
    solution are dropped. Without this the user often sees N
    near-identical rows when the score landscape has a wide plateau.

    Any of ``dra_arcsec``, ``ddec_arcsec``, ``dpa_deg`` that is ≤ 0
    freezes the corresponding axis: scipy's
    ``differential_evolution`` doesn't accept zero-width bounds, so we
    drop the frozen variable from the optimisation and patch it back
    in afterwards.
    """
    if weights is None:
        weights = np.ones_like(evaluator.ra)
    weights = np.asarray(weights, dtype=float)
    cos_dec_med = max(np.cos(np.deg2rad(np.median(evaluator.dec))), 1e-3)
    dra_deg = max(dra_arcsec, 0.0) / 3600.0 / cos_dec_med
    ddec_deg = max(ddec_arcsec, 0.0) / 3600.0
    dpa = max(dpa_deg, 0.0)
    use_flux = (objective == "flux")

    # Which axes are searched vs frozen at the candidate value.
    free = [dra_arcsec > 0, ddec_arcsec > 0, dpa_deg > 0]
    widths = [dra_deg, ddec_deg, dpa]

    n_top = int(min(n_top, len(grid_results["score"])))
    refined_scores: list[float] = []
    refined_params: list[np.ndarray] = []

    for i in range(n_top):
        ra0 = float(grid_results["ra"][i])
        dec0 = float(grid_results["dec"][i])
        pa0 = float(grid_results["pa"][i])

        if not any(free):
            # All axes frozen — nothing to optimise; keep the grid value.
            det, tp, _ = evaluator.evaluate(ra0, dec0, pa0)
            s = (float(np.sum(tp * evaluator.flux * weights)) if use_flux
                 else float(np.sum(det * weights)))
            refined_scores.append(s)
            refined_params.append(np.array([ra0, dec0, pa0]))
            if progress_cb is not None:
                progress_cb(i + 1, n_top)
            continue

        # DE only over the free axes; frozen axes are passed in via
        # closure and reconstructed before each `evaluate` call.
        free_idx = [k for k, f in enumerate(free) if f]
        bounds = []
        for k in free_idx:
            base = (ra0, dec0, pa0)[k]
            bounds.append((base - widths[k], base + widths[k]))

        def neg_score(free_params, _free_idx=free_idx,
                      _ra0=ra0, _dec0=dec0, _pa0=pa0):
            ra, dec, pa = _ra0, _dec0, _pa0
            for j, k in enumerate(_free_idx):
                v = float(free_params[j])
                if k == 0:
                    ra = v
                elif k == 1:
                    dec = v
                else:
                    pa = v
            try:
                det, tp, _ = evaluator.evaluate(ra, dec, pa)
            except Exception:
                return 1e6
            if use_flux:
                return -float(np.sum(tp * evaluator.flux * weights))
            return -float(np.sum(det * weights))

        # `seed` is fixed for repeatable optimisation runs in tests.
        result = differential_evolution(
            neg_score, bounds=bounds,
            maxiter=int(maxiter), popsize=10, seed=42, tol=1e-4,
            polish=True,
        )
        # Reconstruct the full (ra, dec, pa) from the DE result + frozen base.
        ra_p, dec_p, pa_p = ra0, dec0, pa0
        for j, k in enumerate(free_idx):
            v = float(result.x[j])
            if k == 0:
                ra_p = v
            elif k == 1:
                dec_p = v
            else:
                pa_p = v
        refined_scores.append(-float(result.fun))
        refined_params.append(np.array([ra_p, dec_p, pa_p]))
        if progress_cb is not None:
            progress_cb(i + 1, n_top)

    # Sort then dedup.
    refined_scores_arr = np.asarray(refined_scores, dtype=float)
    refined_params_arr = np.asarray(refined_params, dtype=float)
    order = np.argsort(-refined_scores_arr)
    refined_scores_arr = refined_scores_arr[order]
    refined_params_arr = refined_params_arr[order]

    ra_tol_deg = dedup_tol[0] / 3600.0 / cos_dec_med
    dec_tol_deg = dedup_tol[1] / 3600.0
    pa_tol_deg = dedup_tol[2]
    keep: list[int] = []
    for i in range(len(refined_scores_arr)):
        ra_i, dec_i, pa_i = refined_params_arr[i]
        is_dup = False
        for j in keep:
            ra_j, dec_j, pa_j = refined_params_arr[j]
            if (abs(ra_i - ra_j) <= ra_tol_deg
                    and abs(dec_i - dec_j) <= dec_tol_deg
                    and abs(((pa_i - pa_j + 180.0) % 360.0) - 180.0) <= pa_tol_deg):
                is_dup = True
                break
        if not is_dup:
            keep.append(i)

    return {
        "score": refined_scores_arr[keep],
        "ra": refined_params_arr[keep, 0],
        "dec": refined_params_arr[keep, 1],
        "pa": refined_params_arr[keep, 2],
    }
