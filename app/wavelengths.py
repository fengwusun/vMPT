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

# Gap parameters (relative to fiducial wavelength span).
GAP_CENTER_REL: float = 0.50
GAP_WIDTH_REL: float = 0.10


def cutoffs(v2_arcsec: float, v3_arcsec: float, disperser: str, filt: str) -> dict:
    disperser = disperser.upper()
    filt = filt.upper()
    if disperser not in GRATING_RANGES or filt not in GRATING_RANGES[disperser]:
        raise ValueError(f"Unsupported (disperser, filter) = ({disperser}, {filt})")

    lam_min, lam_max = GRATING_RANGES[disperser][filt]
    span = lam_max - lam_min
    dlam_dv2 = span / V2_DISP_EXTENT

    # Wavelength shift induced by shutter V2 offset from the fiducial.
    shift = (v2_arcsec - MSA_V2_REF) * dlam_dv2

    lam_blue = lam_min + shift
    lam_red = lam_max + shift
    gap_half = 0.5 * GAP_WIDTH_REL * span
    lam_gap_center = lam_min + GAP_CENTER_REL * span + shift
    lam_gap_lo = lam_gap_center - gap_half
    lam_gap_hi = lam_gap_center + gap_half

    blue_cut = FILTER_BLUE_CUTOFF.get(filt, 0.0)

    def _maybe(lam: float, is_blue_edge: bool = False) -> float | None:
        if lam < blue_cut:
            # Blue edge clamps to filter cutoff if any red part of spectrum survives.
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
