"""MSA pointing optimizer.

Searches for (RA, Dec, V3 PA) maximising the number — or weighted flux —
of catalog sources that fall in operable, well-centred MSA shutters.

The algorithm is derived from **hMPT** by Daniel Eisenstein, Samuel
McCarty, and Zihao Wu (Harvard / CfA), itself a Python re-implementation
of ESA's eMPT (Bonaventura et al. 2023, A&A 672, A40):
<https://github.com/zihaowu-astro/hMPT>.

vMPT does NOT depend on hMPT — this module re-implements the core
algorithm in vMPT's style with attribution, so it composes cleanly
with our existing MSA grid (`vmpt/data/nirspec_msa_v2v3.npz`), CRDS
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
from .wavelengths import (
    cutoffs as _wavelength_cutoffs,
    disperser_range,
    interval_covered,
    v2_overlap_distance,
)


# Drop-reason codes used by the per-pointing reason tally returned by
# :meth:`PointingEvaluator.evaluate_with_reasons`. Stable identifiers
# so downstream callers (results modal, tests) can key off them.
DROP_COLLISION    = "collision"
DROP_REQUIRED_LAM = "required_lam"
DROP_NO_GAP       = "no_gap"
DROP_EXTEND_BLUE  = "extend_blue"
DROP_EXTEND_RED   = "extend_red"
DROP_REASONS = (
    DROP_COLLISION,
    DROP_REQUIRED_LAM,
    DROP_NO_GAP,
    DROP_EXTEND_BLUE,
    DROP_EXTEND_RED,
)


# Physical shutter dimensions on the focal plane (arcsec).
# `SHUTTER_X` is the dispersion direction (columns, 0..364); `SHUTTER_Y`
# is the spatial direction (rows, 0..170). Values from hMPT, which
# matches APT's MSA model.
SHUTTER_X_ARCSEC = 0.2679
SHUTTER_Y_ARCSEC = 0.5294

# Maximum |Δs| between two individual shutters for them to be
# considered on the same detector y-row (and thus possible spectral
# collisions). eMPT uses shval ≈ s exactly; we allow ±1 to be on the
# safe side. Matches `SHVAL_S_TOLERANCE` in `vmpt/main.py`
# (live-canvas orange overlap).
#
# NOTE: this is the *per-shutter* tolerance. The optimizer's
# collision protection compares two source positions (each opens an
# N-shutter slitlet centred on the source row), so it uses a wider
# slitlet-aware tolerance computed once in
# `PointingEvaluator._init_protection`: `half + 1` for protected ↔
# stuck-open and `2·half + 1` for protected ↔ slitlet-source, with
# `half = slit_length // 2`.
SHVAL_S_TOLERANCE: int = 1

# Detector-half assignment: Q1+Q3 → NRS1, Q2+Q4 → NRS2. Cross-half
# pairs image onto different detectors and therefore never overlap.
NRS1_QUADS = frozenset({1, 3})
NRS2_QUADS = frozenset({2, 4})

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
    protect_mask : ndarray of bool, optional
        Parallel to ``ra_sources``; True marks a source whose spectrum
        must be protected from same-row collisions under the current
        (disperser, filter). Requires ``disperser`` + ``filt`` to be
        meaningful. When None or all-False, no protection is applied
        and ``evaluate`` behaves exactly as before.
    priorities, weights : ndarray, optional
        Per-source priority / weight, used only to break ties when two
        protected sources collide (lower priority number wins; on tie,
        higher weight wins). NaN-tolerant. Falls back to source index
        order if neither is provided.
    disperser, filt : str, optional
        e.g. ``"PRISM"`` and ``"CLEAR"``. Only consulted when
        ``protect_mask`` flags any source. The V2 half-extent of the
        spectrum is looked up from :func:`v2_overlap_distance`.
    reason : ndarray, optional
        (4, 171, 365) operability-reason array from
        :func:`vmpt.msa.load_operability`. Cells equal to 2 are
        stuck-open shutters, which act as always-on dispersion sources
        even when no slitlet is opened there. When provided AND
        protection is enabled, a protected source landing on a row
        colliding with any stuck-open shutter is dropped (its spectrum
        is unavoidably contaminated).
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
        *,
        protect_mask: Optional[np.ndarray] = None,
        priorities: Optional[np.ndarray] = None,
        weights: Optional[np.ndarray] = None,
        disperser: Optional[str] = None,
        filt: Optional[str] = None,
        reason: Optional[np.ndarray] = None,
        # Per-target spectral constraints (v1.3.0+). Each is a length-N
        # array parallel to ra_sources/dec_sources. Defaults preserve
        # v1.2.x behaviour. See :meth:`_apply_constraint_drops` for the
        # exact rules.
        required_lam: Optional[np.ndarray] = None,
        no_gap: Optional[np.ndarray] = None,
        extend_blue: Optional[np.ndarray] = None,
        extend_red: Optional[np.ndarray] = None,
        protect: Optional[np.ndarray] = None,
        # v1.3.1+: per-target centration override. Length-N array of
        # strings (or None). Each non-empty cell wins **unconditionally**
        # over the global ``centration`` argument for that one source —
        # even when it's a laxer level. Empty / None cells use the global.
        centration_per_target: Optional[np.ndarray] = None,
    ):
        self.ra = np.asarray(ra_sources, dtype=float)
        self.dec = np.asarray(dec_sources, dtype=float)
        self.flux = (np.ones_like(self.ra) if flux_sources is None
                     else np.asarray(flux_sources, dtype=float))
        self.sigma = float(sigma_arcsec)
        self.buffer = CENTRATION_BUFFERS.get(
            centration.upper(), CENTRATION_BUFFERS["UNCONSTRAINED"])
        # Per-target centration buffer. Defaults to a length-N array
        # filled with the global ``self.buffer``; any source with a
        # non-empty entry in ``centration_per_target`` overrides its
        # cell to that level's buffer. ``_check_centration`` consumes
        # ``self.buffer_per_source`` (vector) instead of ``self.buffer``
        # (scalar); the latter is retained only for tests / introspection.
        n_src = len(self.ra)
        self.buffer_per_source = np.full(n_src, self.buffer, dtype=float)
        if centration_per_target is not None:
            arr = np.asarray(centration_per_target, dtype=object)
            if arr.size != n_src and arr.size != 0:
                raise ValueError(
                    f"centration_per_target size {arr.size} != "
                    f"ra_sources size {n_src}"
                )
            for i in range(min(n_src, arr.size)):
                v = arr[i]
                if v is None:
                    continue
                s = str(v).strip()
                if not s:
                    continue
                buf = CENTRATION_BUFFERS.get(s.upper())
                if buf is not None:
                    self.buffer_per_source[i] = buf
                # Unrecognised label → silently leaves the global
                # buffer in place (defensive — the loader/UI already
                # normalises, this is just belt-and-braces).
        self.slit_length = int(slit_length)
        if operable is None:
            operable_loaded, reason_loaded = load_operability()
            operable = operable_loaded
            if reason is None:
                reason = reason_loaded
        self.operable = np.asarray(operable, dtype=bool)
        self.interpolators = _build_inverse_interpolators()
        # Effective protect mask = (catalog-wide v1.2 cutoff)
        #                       ∪ (per-target editor flag).
        # Either source making a row protected enables the v1.2.x
        # collision-protection rules for that row.
        effective_protect = self._merge_protect_masks(protect_mask, protect)
        self._init_protection(
            protect_mask=effective_protect,
            priorities=priorities,
            weights=weights,
            disperser=disperser,
            filt=filt,
            reason=reason,
        )
        self._init_spectral_constraints(
            required_lam=required_lam,
            no_gap=no_gap,
            extend_blue=extend_blue,
            extend_red=extend_red,
            disperser=disperser,
            filt=filt,
        )

    # -- per-target protect helpers ---------------------------------

    def _merge_protect_masks(
        self,
        catalog_wide: Optional[np.ndarray],
        per_target: Optional[np.ndarray],
    ) -> Optional[np.ndarray]:
        """Build the effective protect mask = OR of the two inputs.

        Returns None when both inputs are None or all-False — the
        common case, which skips :meth:`_init_protection`'s heavier
        setup entirely.
        """
        n = len(self.ra)
        a = (np.zeros(n, dtype=bool) if catalog_wide is None
             else np.asarray(catalog_wide, dtype=bool))
        b = (np.zeros(n, dtype=bool) if per_target is None
             else np.asarray(per_target, dtype=bool))
        if a.size != n and a.size != 0:
            raise ValueError(
                f"protect_mask size {a.size} != ra_sources size {n}"
            )
        if b.size != n and b.size != 0:
            raise ValueError(
                f"per-target protect size {b.size} != ra_sources size {n}"
            )
        if a.size == 0:
            a = np.zeros(n, dtype=bool)
        if b.size == 0:
            b = np.zeros(n, dtype=bool)
        merged = a | b
        return merged if merged.any() else None

    # -- collision-protection setup ---------------------------------

    def _init_protection(
        self,
        protect_mask: Optional[np.ndarray],
        priorities: Optional[np.ndarray],
        weights: Optional[np.ndarray],
        disperser: Optional[str],
        filt: Optional[str],
        reason: Optional[np.ndarray],
    ) -> None:
        """Cache everything the per-pointing collision check needs.

        Lazily skipped when ``protect_mask`` is None or all-False — the
        common case (no protection) pays no construction cost beyond
        a couple of cheap None assignments.
        """
        self._protect_enabled = False
        self._protect_mask = None
        self._collision_rank = None
        self._v2_lut = None
        self._v2_overlap = 0.0
        self._stuck_open_half = np.empty(0, dtype=np.int8)
        self._stuck_open_s = np.empty(0, dtype=np.int32)
        self._stuck_open_v2 = np.empty(0, dtype=float)
        # Slitlet-aware spatial tolerances. Computed once from
        # `self.slit_length`. A protected target at row s_p occupies
        # 2*half+1 shutters; the user wants ALSO the rows s_p±half±1
        # to be empty so no neighbouring shutter (own- or stuck-open)
        # disperses onto the same detector row.
        #
        #   For protected ↔ stuck-open (single shutter at row s_o):
        #       collide iff |s_o − s_p| ≤ half + 1
        #
        #   For protected ↔ another-source-with-slitlet at row s_q
        #   (both with the same slit_length, so half_q = half_p):
        #       collide iff |s_q − s_p| ≤ 2·half + 1
        #
        # (Slitlets touching one row outside each other = one row of
        # mutual buffer; the +1 in each formula encodes that buffer.)
        self._slit_half = self.slit_length // 2
        self._sd_tol_ps = self._slit_half + 1
        self._sd_tol_pp = 2 * self._slit_half + 1

        if protect_mask is None:
            return
        pm = np.asarray(protect_mask, dtype=bool)
        if pm.size != len(self.ra):
            raise ValueError(
                "protect_mask size must match ra_sources "
                f"({pm.size} vs {len(self.ra)})"
            )
        if not pm.any():
            return
        if disperser is None or filt is None:
            raise ValueError(
                "protect_mask is set; disperser and filt are required "
                "to look up the spectral-collision V2 half-extent."
            )

        self._protect_enabled = True
        self._protect_mask = pm

        # Build a tie-break rank: smaller = wins. Primary key is
        # priority (NaN sinks to the back), secondary is -weight (so
        # higher weight wins), tertiary is index (stable).
        n = len(self.ra)
        pri = (np.asarray(priorities, dtype=float) if priorities is not None
               else np.full(n, np.nan))
        wgt = (np.asarray(weights, dtype=float) if weights is not None
               else np.full(n, np.nan))
        pri = np.where(np.isnan(pri), np.inf, pri)
        wgt = np.where(np.isnan(wgt), -np.inf, wgt)
        order = np.lexsort((np.arange(n), -wgt, pri))
        self._collision_rank = np.empty(n, dtype=np.int64)
        self._collision_rank[order] = np.arange(n, dtype=np.int64)

        # V2 lookup per (q, s, d). (4, 171, 365) array of floats.
        v2_grid, _ = load_msa_grid()
        self._v2_lut = np.asarray(v2_grid, dtype=float)

        # Spectral overlap half-extent for the requested (disperser,
        # filter). Falls back to a conservative default inside
        # `v2_overlap_distance` for unknown configs.
        self._v2_overlap = float(v2_overlap_distance(disperser, filt))

        # Stuck-open shutters (REASON == 2). Cache their (det_half, s,
        # V2). When `reason` isn't provided we leave the arrays empty;
        # collision protection then ignores stuck-opens.
        if reason is not None:
            r_arr = np.asarray(reason, dtype=np.int8)
            stuck = np.argwhere(r_arr == 2)
            if stuck.size > 0:
                q_idx = stuck[:, 0]  # 0..3
                s_idx = stuck[:, 1]  # 0..170
                d_idx = stuck[:, 2]  # 0..364
                # NRS1 = quads 1,3 → q_idx 0,2; NRS2 = quads 2,4 → 1,3.
                half = np.where(np.isin(q_idx, [0, 2]), 1, 2).astype(np.int8)
                v2_vals = self._v2_lut[q_idx, s_idx, d_idx].astype(float)
                self._stuck_open_half = half
                self._stuck_open_s = s_idx.astype(np.int32)
                self._stuck_open_v2 = v2_vals

    # -- spectral-constraint setup ----------------------------------

    def _init_spectral_constraints(
        self,
        required_lam: Optional[np.ndarray],
        no_gap: Optional[np.ndarray],
        extend_blue: Optional[np.ndarray],
        extend_red: Optional[np.ndarray],
        disperser: Optional[str],
        filt: Optional[str],
    ) -> None:
        """Cache per-target spectral constraints.

        Lazily skipped when no target has a constraint set — the
        common case (and v1.2.x behaviour) pays no construction cost.
        """
        n = len(self.ra)
        self._constraint_enabled = False
        self._required_lam = None
        self._no_gap = np.zeros(n, dtype=bool)
        self._extend_blue = np.zeros(n, dtype=bool)
        self._extend_red = np.zeros(n, dtype=bool)
        self._constraint_disperser = None
        self._constraint_filt = None
        self._disperser_lam_lo = None
        self._disperser_lam_hi = None

        # Per-target arrays. NaN-tolerant: a missing array stays as
        # the all-False default; a mismatched-size array raises.
        def _bool_array(arr, name):
            if arr is None:
                return np.zeros(n, dtype=bool)
            a = np.asarray(arr, dtype=bool)
            if a.size != n:
                raise ValueError(
                    f"{name} size {a.size} != ra_sources size {n}"
                )
            return a

        if no_gap is not None:
            self._no_gap = _bool_array(no_gap, "no_gap")
        if extend_blue is not None:
            self._extend_blue = _bool_array(extend_blue, "extend_blue")
        if extend_red is not None:
            self._extend_red = _bool_array(extend_red, "extend_red")

        # required_lam is ragged (list of (lo, hi) tuples per source).
        if required_lam is not None:
            rl = np.asarray(required_lam, dtype=object)
            if rl.size != n:
                raise ValueError(
                    f"required_lam size {rl.size} != ra_sources size {n}"
                )
            # Coerce each entry to a list of (float, float) tuples for
            # uniform downstream access; sanitise garbage entries.
            self._required_lam = np.empty(n, dtype=object)
            for i in range(n):
                entry = rl[i]
                if entry is None or (isinstance(entry, float)
                                     and np.isnan(entry)):
                    self._required_lam[i] = []
                    continue
                try:
                    cleaned = [
                        (float(lo), float(hi)) for lo, hi in entry
                        if np.isfinite(float(lo)) and np.isfinite(float(hi))
                    ]
                except (TypeError, ValueError):
                    cleaned = []
                self._required_lam[i] = cleaned
        else:
            # `np.array([[] for ...], dtype=object)` produces a 2D
            # array of shape (n, 0) — element access then returns an
            # empty 1D numpy array, and `bool()` of an empty numpy
            # array warns. Build a true 1D object array explicitly.
            self._required_lam = np.empty(n, dtype=object)
            for i in range(n):
                self._required_lam[i] = []

        # Any constraint flagged on at least one target?
        has_required = any(len(r) > 0 for r in self._required_lam)
        has_any = (
            has_required
            or self._no_gap.any()
            or self._extend_blue.any()
            or self._extend_red.any()
        )
        if not has_any:
            return
        if disperser is None or filt is None:
            raise ValueError(
                "per-target spectral constraints are set; disperser "
                "and filt are required to evaluate them."
            )

        self._constraint_enabled = True
        self._constraint_disperser = str(disperser).upper()
        self._constraint_filt = str(filt).upper()
        # Cache the disperser/filter nominal range for the extend_blue
        # and extend_red checks. None when the combo isn't recognised —
        # extend_* constraints then always fail (which is what the user
        # would want: "extend to the bluest of an unsupported combo"
        # has no satisfiable answer).
        rng = disperser_range(self._constraint_disperser,
                              self._constraint_filt)
        if rng is not None:
            self._disperser_lam_lo, self._disperser_lam_hi = rng

        # v1.3.2+ tolerant filter on `required_lam` ranges:
        # silently drop any user-supplied (lo, hi) range that lies
        # entirely outside the current disperser/filter's wavelength
        # bounds. Rationale — a user editing constraints often
        # pre-stages for a different disperser (say PRISM at 1.0–1.2
        # μm) and would otherwise see every such source dropped under
        # G395H. We treat impossible-under-this-grating ranges as
        # "no constraint" rather than "always fails". When a range
        # only PARTIALLY overlaps the disperser, we keep it — the
        # standard `interval_covered` check then enforces the
        # achievable portion.
        self._required_lam_dropped = 0
        if rng is not None and self._required_lam is not None:
            lo_d, hi_d = rng
            for i in range(n):
                orig = self._required_lam[i]
                if not orig:
                    continue
                kept = [
                    (lo, hi) for (lo, hi) in orig
                    # Overlap test: range [lo, hi] overlaps the
                    # disperser [lo_d, hi_d] iff lo < hi_d AND hi > lo_d.
                    if (lo < hi_d and hi > lo_d)
                ]
                if len(kept) != len(orig):
                    self._required_lam_dropped += (
                        len(orig) - len(kept)
                    )
                    self._required_lam[i] = kept
        # Scan the precomputed per-shutter dispersion table to find
        # the actual "best" lam_blue and lam_red ACHIEVABLE for this
        # (disp, filt) across the MSA. extend_blue passes iff the
        # source's shutter lam_blue is within EDGE_TOL of this
        # MSA-best value (so per-shutter variation doesn't spuriously
        # fail every source). Falls back to the nominal range if the
        # table isn't loadable.
        self._table_best_lam_blue = self._disperser_lam_lo
        self._table_best_lam_red = self._disperser_lam_hi
        if (self._extend_blue.any() or self._extend_red.any()):
            from .wavelengths import _load_dispersion_table
            tbl = _load_dispersion_table()
            if tbl is not None:
                key_blue = f"{self._constraint_disperser}_{self._constraint_filt}_blue_edge"
                key_red  = f"{self._constraint_disperser}_{self._constraint_filt}_red_edge"
                if key_blue in tbl and key_red in tbl:
                    arr_b = np.asarray(tbl[key_blue], dtype=float)
                    arr_r = np.asarray(tbl[key_red], dtype=float)
                    # NaN-aware min/max — np.nanmin handles all-NaN
                    # by raising; guard with finite check.
                    if np.isfinite(arr_b).any():
                        self._table_best_lam_blue = float(np.nanmin(arr_b))
                    if np.isfinite(arr_r).any():
                        self._table_best_lam_red = float(np.nanmax(arr_r))

    # -- evaluation --------------------------------------------------

    def evaluate(
        self,
        ra_p: float, dec_p: float, pa_v3: float,
        theta_deg: float = 90.0,
    ) -> tuple[np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Return ``(detected_bool, throughput, (quad, s, d))`` per source.

        When collision protection was configured at construction time,
        ``detected`` is the **kept** mask — sources dropped by the
        protection rules are zeroed in both ``detected`` and
        ``throughput``. The raw pre-drop mask is not returned; use
        :meth:`evaluate_with_stats` if you also need the drop count.
        """
        kept, tp, idx, _ = self.evaluate_with_stats(
            ra_p, dec_p, pa_v3, theta_deg=theta_deg,
        )
        return kept, tp, idx

    def evaluate_with_stats(
        self,
        ra_p: float, dec_p: float, pa_v3: float,
        theta_deg: float = 90.0,
    ) -> tuple[np.ndarray, np.ndarray,
               tuple[np.ndarray, np.ndarray, np.ndarray], int]:
        """Like :meth:`evaluate` plus an ``n_dropped`` count of sources
        that landed in operable, centred shutters but were excluded by
        either the collision-protection rules (v1.2.0+) or the per-
        target spectral constraints (v1.3.0+) at this pointing.

        ``n_dropped == 0`` when neither protection nor constraints are
        configured. See :meth:`evaluate_with_reasons` for a per-reason
        breakdown.
        """
        kept, tp, idx, reasons = self.evaluate_with_reasons(
            ra_p, dec_p, pa_v3, theta_deg=theta_deg,
        )
        return kept, tp, idx, int(sum(reasons.values()))

    def evaluate_with_reasons(
        self,
        ra_p: float, dec_p: float, pa_v3: float,
        theta_deg: float = 90.0,
    ) -> tuple[np.ndarray, np.ndarray,
               tuple[np.ndarray, np.ndarray, np.ndarray], dict]:
        """Like :meth:`evaluate_with_stats` but returns a dict of
        per-reason drop counts instead of just a scalar.

        Keys are the constants in :data:`DROP_REASONS`
        (``collision``, ``required_lam``, ``no_gap``, ``extend_blue``,
        ``extend_red``). Values sum to the same scalar
        ``evaluate_with_stats`` returns.

        Empty dict when neither protection nor constraints are
        configured.
        """
        axy = radec_to_axy(self.ra, self.dec, ra_p, dec_p, pa_v3, theta_deg)
        quad, s_frac, d_frac = axy_to_shutter(axy, self.interpolators)
        operable_mask = self._check_operable(quad, s_frac, d_frac)
        centered = self._check_centration(s_frac, d_frac)
        with np.errstate(invalid="ignore"):
            tp = _gaussian_through_shutter(s_frac, d_frac, self.sigma)
        tp = np.where(operable_mask & centered, tp, 0.0)
        detected = tp > 0
        reasons: dict = {r: 0 for r in DROP_REASONS}

        if not (self._protect_enabled or self._constraint_enabled):
            return detected, tp, (quad, s_frac, d_frac), reasons

        kept = detected.copy()
        # Collision rules first (their losers don't get re-checked by
        # the spectral rules — they're already dropped).
        if self._protect_enabled:
            after_collision = self._apply_collision_drops(
                kept, quad, s_frac, d_frac,
            )
            reasons[DROP_COLLISION] = int(kept.sum() - after_collision.sum())
            kept = after_collision
        # Spectral constraints (per target).
        if self._constraint_enabled:
            kept, spec_reasons = self._apply_constraint_drops(
                kept, quad, s_frac, d_frac,
            )
            for k, v in spec_reasons.items():
                reasons[k] = reasons.get(k, 0) + int(v)
        tp_kept = np.where(kept, tp, 0.0)
        return kept, tp_kept, (quad, s_frac, d_frac), reasons

    # -- collision rules ---------------------------------------------

    def _apply_collision_drops(
        self,
        detected: np.ndarray,
        quad: np.ndarray,
        s_frac: np.ndarray,
        d_frac: np.ndarray,
    ) -> np.ndarray:
        """Apply the three collision-protection rules at one pointing.

        Every check is **explicitly per-shutter**: each protected
        source's slitlet is expanded into its N constituent shutters
        and every one of them is checked against the other party. Two
        individual shutters collide on the detector when all three of
        the following hold:

        - Same detector half (Q1+Q3 → NRS1; Q2+Q4 → NRS2; cross-half
          pairs image onto different detectors).
        - Same row to within ``SHVAL_S_TOLERANCE = 1`` (the
          per-individual-shutter eMPT convention).
        - V2 separation < :func:`v2_overlap_distance` for the current
          (disperser, filter) — i.e. their dispersed spectra share
          some detector x-pixel range.

        For a protected target with slit_length=N, the slitlet covers
        rows ``{s_p − half, …, s_p + half}`` (``half = N // 2``). The
        slitlet has a collision iff **any one** of its N shutters
        collides with **any one** of the other side's shutters
        (a single shutter for stuck-open, a slitlet for another
        source). This drops the historical "widened tolerance"
        shortcut from v1.2.1 — same math for column-aligned slitlets
        (intra-slitlet V2 variation ≲ 0.4″ ≪ V2-overlap, so the
        shortcut was numerically equivalent), but the explicit form
        uses the actual per-shutter V2 from the MSA grid and reads
        directly as "every slitlet member is checked".

        Rules:

        1. **Protected ↔ stuck-open**: a protected source whose
           slitlet has any shutter colliding with any stuck-open
           shutter is dropped (its spectrum is unavoidably
           contaminated). Stuck-open is a single-shutter dispersion
           source independent of pointing.
        2. **Protected ↔ protected**: when two protected sources'
           slitlets collide, the lower-rank one (higher priority
           number → lower weight → higher index) is dropped. The
           winner stays and continues to provide collision pressure
           on rule 3.
        3. **Protected ↔ unprotected**: any detected unprotected
           source whose slitlet collides with a still-kept protected
           source's slitlet is dropped.

        Returns the kept (post-drop) boolean mask. Dropped protected
        sources do **not** provide collision pressure for steps 2/3 —
        if a high-priority spectrum is already contaminated we won't
        compound the loss by also blocking the unprotected sources.
        """
        kept = detected.copy()
        if not kept.any():
            return kept

        det_half = np.zeros(len(kept), dtype=np.int8)
        is_q1q3 = (quad == 1) | (quad == 3)
        is_q2q4 = (quad == 2) | (quad == 4)
        det_half[is_q1q3] = 1
        det_half[is_q2q4] = 2

        # ── Build per-source slitlet arrays (rows + V2) ─────────────
        # Each detected source occupies N=2*half+1 shutters at rows
        # s_int + k for k in [-half, +half]. We fetch the V2 of each
        # of those shutters individually — they differ by ~0.4″ across
        # the N rows (MSA is rotated relative to V2/V3). Out-of-MSA
        # slitlet members (s_offset < 0 or ≥ 171) carry NaN, which
        # makes any subsequent V2 comparison return False (NumPy NaN
        # semantics) — i.e. they're treated as "no shutter exists
        # there to collide with", which is correct.
        half = self._slit_half
        n_slit = 2 * half + 1
        n_sources = len(kept)
        slit_rows = np.full((n_sources, n_slit), -10_000, dtype=np.int32)
        slit_v2 = np.full((n_sources, n_slit), np.nan, dtype=float)

        in_grid_full = (quad > 0) & detected
        active_idx = np.where(in_grid_full)[0]
        if active_idx.size > 0:
            qi = quad[active_idx] - 1  # 0-based, (n_active,)
            si = np.rint(s_frac[active_idx]).astype(int)
            di = np.rint(d_frac[active_idx]).astype(int)
            # Reject sources whose centre column is off-grid; their
            # slitlet rows still need bounds-checked below.
            di_ok = (di >= 0) & (di < 365)
            for k, ds in enumerate(range(-half, half + 1)):
                s_off = si + ds
                row_ok = di_ok & (s_off >= 0) & (s_off < 171)
                # Where the slitlet shutter exists on the MSA, record
                # its row and V2; else leave the (−10_000, NaN)
                # sentinels in place.
                good = active_idx[row_ok]
                if good.size > 0:
                    slit_rows[good, k] = s_off[row_ok]
                    slit_v2[good, k] = self._v2_lut[qi[row_ok],
                                                    s_off[row_ok],
                                                    di[row_ok]]

        # -------- Rule 1: protected slitlet ↔ stuck-open --------
        if self._stuck_open_v2.size > 0:
            prot_idx = np.where(kept & self._protect_mask & (det_half > 0))[0]
            if prot_idx.size > 0:
                # Broadcast shape: (n_prot, n_slit, n_stuck)
                #   ph        : (n_prot, 1, 1)
                #   slit_rows : (n_prot, n_slit, 1)
                #   slit_v2   : (n_prot, n_slit, 1)
                #   so_h      : (1, 1, n_stuck)
                #   so_s      : (1, 1, n_stuck)
                #   so_v2     : (1, 1, n_stuck)
                ph = det_half[prot_idx][:, None, None]
                rows_p = slit_rows[prot_idx][:, :, None]
                v2_p = slit_v2[prot_idx][:, :, None]
                so_h = self._stuck_open_half[None, None, :]
                so_s = self._stuck_open_s[None, None, :]
                so_v2 = self._stuck_open_v2[None, None, :]
                collide = (
                    (ph == so_h)
                    & (np.abs(rows_p - so_s) <= SHVAL_S_TOLERANCE)
                    & (np.abs(v2_p - so_v2) < self._v2_overlap)
                )
                # Drop if ANY (slit shutter, stuck shutter) pair collides.
                kept[prot_idx[collide.any(axis=(1, 2))]] = False

        # -------- Rule 2: protected slitlet ↔ protected slitlet --------
        kept_prot_idx = np.where(kept & self._protect_mask & (det_half > 0))[0]
        if kept_prot_idx.size >= 2:
            rank = self._collision_rank[kept_prot_idx]
            order = np.argsort(rank)
            ordered = kept_prot_idx[order]
            alive = np.ones(ordered.size, dtype=bool)
            for k in range(ordered.size):
                if not alive[k]:
                    continue
                wi = ordered[k]
                tail = ordered[k + 1:]
                tail_alive = alive[k + 1:]
                if tail.size == 0 or not tail_alive.any():
                    continue
                live_tail_idx = np.where(tail_alive)[0]
                lt = tail[live_tail_idx]
                # Broadcast: (n_lt, n_slit_l, n_slit_w)
                w_rows = slit_rows[wi]            # (n_slit,)
                w_v2 = slit_v2[wi]                # (n_slit,)
                w_half = det_half[wi]             # scalar
                l_rows = slit_rows[lt][:, :, None]  # (n_lt, n_slit_l, 1)
                l_v2 = slit_v2[lt][:, :, None]    # (n_lt, n_slit_l, 1)
                l_half = det_half[lt][:, None, None]  # (n_lt, 1, 1)
                collide = (
                    (l_half == w_half)
                    & (np.abs(l_rows - w_rows[None, None, :])
                       <= SHVAL_S_TOLERANCE)
                    & (np.abs(l_v2 - w_v2[None, None, :])
                       < self._v2_overlap)
                )
                drops = collide.any(axis=(1, 2))
                if drops.any():
                    losers = lt[drops]
                    kept[losers] = False
                    alive[k + 1 + live_tail_idx[drops]] = False

        # -------- Rule 3: protected slitlet ↔ unprotected slitlet --------
        kept_prot_idx = np.where(kept & self._protect_mask & (det_half > 0))[0]
        if kept_prot_idx.size > 0:
            unprot_idx = np.where(
                kept & ~self._protect_mask & (det_half > 0)
            )[0]
            if unprot_idx.size > 0:
                # Broadcast: (n_u, n_slit_u, n_p, n_slit_p)
                u_rows = slit_rows[unprot_idx][:, :, None, None]
                u_v2 = slit_v2[unprot_idx][:, :, None, None]
                u_half = det_half[unprot_idx][:, None, None, None]
                p_rows = slit_rows[kept_prot_idx][None, None, :, :]
                p_v2 = slit_v2[kept_prot_idx][None, None, :, :]
                p_half = det_half[kept_prot_idx][None, None, :, None]
                collide = (
                    (u_half == p_half)
                    & (np.abs(u_rows - p_rows) <= SHVAL_S_TOLERANCE)
                    & (np.abs(u_v2 - p_v2) < self._v2_overlap)
                )
                # Drop unprot if any (its shutter, prot's shutter)
                # pair collides for any prot source.
                kept[unprot_idx[collide.any(axis=(1, 2, 3))]] = False

        return kept

    # -- spectral-constraint rules ----------------------------------

    def _apply_constraint_drops(
        self,
        detected: np.ndarray,
        quad: np.ndarray,
        s_frac: np.ndarray,
        d_frac: np.ndarray,
    ) -> tuple[np.ndarray, dict]:
        """Apply the per-target spectral constraints (v1.3.0+).

        For every detected source with at least one of
        ``required_lam``, ``no_gap``, ``extend_blue``, ``extend_red``
        set, look up the centre shutter's wavelength endpoints via
        :func:`vmpt.wavelengths.cutoffs` and drop the source if any
        flagged constraint fails.

        Returns ``(kept, reasons)`` where ``reasons`` is a dict
        keyed by the relevant ``DROP_*`` constants. Each source can
        only be dropped once — the first constraint that fails wins
        the bookkeeping. (We still check all constraints for one
        source so a future "explain why this source dropped" UI has
        access to the full list; that needs the per-source detail
        which we don't currently expose.)
        """
        reasons = {
            DROP_REQUIRED_LAM: 0,
            DROP_NO_GAP: 0,
            DROP_EXTEND_BLUE: 0,
            DROP_EXTEND_RED: 0,
        }
        if not self._constraint_enabled:
            return detected, reasons

        kept = detected.copy()
        disp = self._constraint_disperser
        filt = self._constraint_filt
        lam_lo = self._disperser_lam_lo
        lam_hi = self._disperser_lam_hi
        # Tolerance (μm) for the extend_blue / extend_red comparison.
        # 20 nm absorbs per-shutter wavelength-solution variation —
        # for PRISM typical centre shutters land at 0.604 while the
        # table-wide minimum is 0.600 (a few resolution elements'
        # spread). Edge truncation pushes lam_blue much further red
        # than this so the constraint still fires in the case it's
        # designed to catch.
        EDGE_TOL = 0.020

        in_grid = (quad > 0) & detected
        if not in_grid.any():
            return kept, reasons
        idx = np.where(in_grid)[0]

        for i in idx:
            req = (self._required_lam[i]
                   if self._required_lam is not None else [])
            # `req` is a python list of (lo, hi) tuples; bool(list)
            # is safe. The bool fields are numpy scalars — cast
            # explicitly so future NumPy "ambiguous truthiness"
            # rules don't bite us.
            need_no_gap = bool(self._no_gap[i])
            need_blue = bool(self._extend_blue[i])
            need_red = bool(self._extend_red[i])
            # `len()` is safe on both Python lists and NumPy arrays;
            # avoids `bool(arr)`'s ambiguous-truthiness deprecation.
            has_req = (req is not None
                       and hasattr(req, "__len__") and len(req) > 0)
            if not (has_req or need_no_gap or need_blue or need_red):
                continue

            q = int(quad[i])
            s_int = int(round(float(s_frac[i])))
            d_int = int(round(float(d_frac[i])))
            # cutoffs() takes 1-based shutter indices.
            try:
                bounds = _wavelength_cutoffs(
                    0.0, 0.0, disp, filt,
                    q=q, s=s_int + 1, d=d_int + 1,
                )
            except Exception:  # noqa: BLE001
                bounds = None
            if bounds is None:
                # Disperser/filter combination not in the per-shutter
                # table — every spectral constraint fails by default.
                kept[i] = False
                # Attribute the drop to the first flagged reason in a
                # stable order.
                if has_req:
                    reasons[DROP_REQUIRED_LAM] += 1
                elif need_no_gap:
                    reasons[DROP_NO_GAP] += 1
                elif need_blue:
                    reasons[DROP_EXTEND_BLUE] += 1
                elif need_red:
                    reasons[DROP_EXTEND_RED] += 1
                continue
            # cutoffs() returns None for values that are NaN OR below
            # the filter blue cutoff. We treat None as NaN throughout
            # the constraint checks — interval_covered etc. handle
            # NaN bounds correctly.
            def _to_f(v):
                return float("nan") if v is None else float(v)
            blue = _to_f(bounds["lam_blue"])
            gap_lo = _to_f(bounds["lam_gap_lo"])
            gap_hi = _to_f(bounds["lam_gap_hi"])
            red = _to_f(bounds["lam_red"])

            dropped_reason = None
            # Required wavelength ranges — every interval must be
            # covered (with the gap excluded). First failure wins.
            if has_req:
                for (lo, hi) in req:
                    if not interval_covered(lo, hi, blue, gap_lo, gap_hi, red):
                        dropped_reason = DROP_REQUIRED_LAM
                        break
            if dropped_reason is None and need_no_gap:
                if np.isfinite(gap_lo) and np.isfinite(gap_hi):
                    dropped_reason = DROP_NO_GAP
            if dropped_reason is None and need_blue:
                # The "extend to bluest" comparison is against the
                # MSA-WIDE best lam_blue cached at init time (not the
                # disperser's nominal range), so per-shutter
                # dispersion variation doesn't spuriously fail every
                # source. None / NaN → can't satisfy → drop.
                ref_blue = self._table_best_lam_blue
                if (ref_blue is None or not np.isfinite(blue)
                        or blue > ref_blue + EDGE_TOL):
                    dropped_reason = DROP_EXTEND_BLUE
            if dropped_reason is None and need_red:
                ref_red = self._table_best_lam_red
                if (ref_red is None or not np.isfinite(red)
                        or red < ref_red - EDGE_TOL):
                    dropped_reason = DROP_EXTEND_RED

            if dropped_reason is not None:
                kept[i] = False
                reasons[dropped_reason] += 1

        return kept, reasons

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
        """Source falls within the centration buffer (in BOTH axes).

        Uses ``self.buffer_per_source`` so per-target overrides
        (v1.3.1+) take effect element-wise; the limits are computed
        per row instead of as scalars. The global ``self.buffer`` is
        unused here — it survives only as the default fill of
        ``self.buffer_per_source``.
        """
        with np.errstate(invalid="ignore"):
            row_limit = 0.5 - (self.buffer_per_source / SHUTTER_Y_ARCSEC)
            col_limit = 0.5 - (self.buffer_per_source / SHUTTER_X_ARCSEC)
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
