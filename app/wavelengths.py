"""Analytic per-grating dispersion model for NIRSpec MSA shutters."""

from .coords import MSA_V2_REF

# Per-grating fiducial wavelength ranges (microns), sourced from JDox
# (jwst-docs.stsci.edu -> NIRSpec -> Dispersers and Filters).
GRATING_RANGES: dict[str, dict[str, tuple[float, float]]] = {
    "PRISM": {"CLEAR": (0.60, 5.30)},
    "G140M": {"F070LP": (0.70, 1.27), "F100LP": (0.97, 1.84)},
    "G235M": {"F170LP": (1.66, 3.07)},
    "G395M": {"F290LP": (2.87, 5.14)},
    "G140H": {"F070LP": (0.81, 1.27), "F100LP": (0.97, 1.84)},
    "G235H": {"F170LP": (1.66, 3.07)},
    "G395H": {"F290LP": (2.87, 5.14)},
}

# Filter blue cutoffs (microns) — hard lower bounds.
FILTER_BLUE_CUTOFF: dict[str, float] = {
    "CLEAR": 0.60,
    "F070LP": 0.70,
    "F100LP": 0.97,
    "F170LP": 1.66,
    "F290LP": 2.87,
}

# Projected MSA-to-detector dispersion span along V2 (arcsec). Placeholder per PLAN.md.
V2_DISP_EXTENT: float = 180.0

LAMBDA_PER_ARCSEC: float = (5.30 - 0.60) / V2_DISP_EXTENT  # ~0.026 μm/arcsec


# Half-width of the on-detector spectrum, in V2 arcsec, per disperser.
#
# Source: eMPT (Bonaventura et al. 2023, A&A 672, A40),
# `reference_files/prism_sep.dat`, which tabulates per-shutter +/- column
# separations for the PRISM spectrum. At a central shutter the values are
# sep_p ≈ +176, sep_m ≈ -177 columns. Each MSA column is ~0.20" along V2,
# so the projected V2 half-extent of PRISM's spectrum is ~0.20 × 176 ≈ 35".
# Within ~10 % across the MSA.
#
# For the grating modes, eMPT applies no column cutoff at all — any pair
# of shutters at the same detector-y row collide (the spectrum spans the
# entire detector). For an H-grating the spectrum reaches beyond the MSA
# in V2, so cross-quadrant pairs collide as well. We approximate both
# behaviours with a single (large) V2 half-extent below.
#
# IMPORTANT — what the visualisation actually shows. We compute the
# orange spectral-conflict shutters as "same q AND same s in the MSA grid,
# within `v2_overlap_distance` of the open shutter in V2." So:
#
#   - PRISM at a center shutter (d ~ 200): ~94 % of the row in the same
#     quadrant. (Yes, "most of the row" — physical: PRISM spectrum on
#     detector is ~70" wide in V2, the row is 73" wide.)
#   - PRISM at an edge shutter (d ~ 10): only ~51 % of the row. (The
#     spectrum runs off the end.)
#   - M / H gratings: ~100 % of the row in same quadrant.
#
# Not yet modelled: cross-quadrant collisions for grating modes. Two
# shutters in different quadrants but at the same detector y can collide;
# we'd need the eMPT shval tables to do this exactly. Reasonable
# enhancement for a future release.
SPECTRUM_V2_HALFEXTENT: dict[str, float] = {
    "PRISM": 35.0,
    "G140M": 200.0,
    "G235M": 200.0,
    "G395M": 200.0,
    "G140H": 500.0,
    "G235H": 500.0,
    "G395H": 500.0,
}


def v2_overlap_distance(disperser: str, filt: str) -> float:
    """V2 half-extent (arcsec) of the spectrum on the detector. Two same-row
    shutters whose V2 separation is less than this value collide on the
    detector. See module-level comment for the eMPT reference."""
    return SPECTRUM_V2_HALFEXTENT.get(disperser.upper(), 180.0)

# Per-disperser NRS1/NRS2 detector-gap wavelengths (microns), at the
# MSA fiducial (central) shutter. The gap is the projection of the
# NRS1/NRS2 detector boundary back through the disperser into the
# wavelength axis.
#
# PRISM/CLEAR — sourced from spacetelescope/msaviz (commit on main
# circa 2026), which numerically integrates the pipeline PRISM
# dispersion polynomial per MSA shutter. At the central Q1 shutter
# (I=311, J=86) the gap edges are (1.87, 3.93) μm. PRISM dispersion
# is highly non-linear so the per-shutter spread is large (gap_lo
# 5-95 % ≈ [0.65, 3.59], gap_hi 5-95 % ≈ [3.03, 5.02] μm); we render
# the fiducial values as an approximation and note the spread in
# the UI tooltip rather than mis-shifting via the linear model.
#
# Gratings — approximate values. The previous code used a uniform
# "gap at 50 % of span, 10 % wide" fudge; we keep the linear shift
# model for gratings (dispersion IS roughly linear there) but at a
# narrower width that better matches the on-detector pixel gap.
# Verifying these against msaviz/calibration tables is a follow-up.
DETECTOR_GAP_FIDUCIAL: dict[str, tuple[float, float]] = {
    "PRISM": (1.87, 3.93),
}

# Fallback fractional gap parameters used when a disperser isn't in
# DETECTOR_GAP_FIDUCIAL. Width tightened from the previous 10 % so
# gratings show a more believably narrow detector gap.
GAP_CENTER_REL: float = 0.50
GAP_WIDTH_REL: float = 0.04


def cutoffs(v2_arcsec: float, v3_arcsec: float, disperser: str, filt: str) -> dict:
    """Wavelength endpoints of the dispersed spectrum on the detector for
    a shutter at (V2, V3). All values are CLAMPED to the disperser's
    intrinsic [lam_min, lam_max] range — we never report wavelengths
    beyond what the grating can usefully observe (so PRISM is capped at
    5.3 μm even for shutters far from V2_REF).

    For PRISM the gap location is held at the fiducial msaviz value
    rather than shifted linearly with V2: PRISM dispersion is too
    non-linear for the linear model to capture, and the fiducial value
    is a much better approximation than the linearly-shifted result.
    For the gratings the linear shift model is retained.
    """
    disperser = disperser.upper()
    filt = filt.upper()
    if disperser not in GRATING_RANGES or filt not in GRATING_RANGES[disperser]:
        raise ValueError(f"Unsupported (disperser, filter) = ({disperser}, {filt})")

    lam_min, lam_max = GRATING_RANGES[disperser][filt]
    span = lam_max - lam_min
    dlam_dv2 = span / V2_DISP_EXTENT

    # Wavelength shift induced by shutter V2 offset from the fiducial.
    # PRISM's non-linear dispersion makes this shift a poor model
    # (msaviz shows the actual per-shutter shift is several × larger
    # in the red than the blue); we use it only for the gratings.
    shift = (v2_arcsec - MSA_V2_REF) * dlam_dv2
    if disperser == "PRISM":
        lam_blue_raw = lam_min
        lam_red_raw = lam_max
    else:
        lam_blue_raw = lam_min + shift
        lam_red_raw = lam_max + shift

    # Clamp the shifted spectrum to the grating's intrinsic range — the
    # detector can't pick up wavelengths the grating doesn't usefully
    # disperse.
    lam_blue = max(lam_min, min(lam_max, lam_blue_raw))
    lam_red = max(lam_min, min(lam_max, lam_red_raw))

    fixed_gap = DETECTOR_GAP_FIDUCIAL.get(disperser)
    if fixed_gap is not None:
        # Gap held at msaviz fiducial values; no V2 shift applied.
        lam_gap_lo_raw, lam_gap_hi_raw = fixed_gap
    else:
        gap_half = 0.5 * GAP_WIDTH_REL * span
        lam_gap_center_raw = lam_min + GAP_CENTER_REL * span + shift
        lam_gap_lo_raw = lam_gap_center_raw - gap_half
        lam_gap_hi_raw = lam_gap_center_raw + gap_half

    # The gap is physical — if the shifted gap falls outside the
    # observable range, there's effectively no detector gap on this
    # shutter's spectrum (it's seen contiguously, or not at all).
    if lam_gap_hi_raw < lam_blue or lam_gap_lo_raw > lam_red:
        lam_gap_lo = None
        lam_gap_hi = None
    else:
        lam_gap_lo = max(lam_blue, min(lam_red, lam_gap_lo_raw))
        lam_gap_hi = max(lam_blue, min(lam_red, lam_gap_hi_raw))

    blue_cut = FILTER_BLUE_CUTOFF.get(filt, 0.0)

    def _maybe(lam, is_blue_edge: bool = False):
        if lam is None:
            return None
        if lam < blue_cut:
            if is_blue_edge and lam_red > blue_cut:
                return float(blue_cut)
            return None
        return float(lam)

    return {
        "lam_blue": _maybe(lam_blue, is_blue_edge=True),
        "lam_gap_lo": _maybe(lam_gap_lo),
        "lam_gap_hi": _maybe(lam_gap_hi),
        "lam_red": _maybe(lam_red),
    }
