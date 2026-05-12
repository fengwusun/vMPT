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

# Per-arcsec wavelength dispersion, calibrated so PRISM/CLEAR's ~4.7 μm range
# spans the full V2_DISP_EXTENT. Filters with narrower wavelength coverage
# (e.g. G140M/F070LP at 0.57 μm) get a *proportionally narrower* V2 overlap
# zone. This is a deliberate approximation for the spectral-conflict visual
# — see PLAN.md M5 for the placeholder notes.
LAMBDA_PER_ARCSEC: float = (5.30 - 0.60) / V2_DISP_EXTENT  # ~0.026 μm/arcsec


def v2_overlap_distance(disperser: str, filt: str) -> float:
    """Half-width V2 distance (arcsec) over which two shutters at the same
    s-row have overlapping wavelength coverage with the given grating/filter.

    For PRISM/CLEAR this is the full V2_DISP_EXTENT (~180″) — almost every
    same-row shutter conflicts. For narrow-filter combos (e.g.
    G140M/F070LP, ~0.57 μm) it shrinks to ~22″ — only nearby shutters share
    any wavelength.
    """
    disperser = disperser.upper()
    filt = filt.upper()
    if disperser not in GRATING_RANGES or filt not in GRATING_RANGES[disperser]:
        return V2_DISP_EXTENT
    lam_min, lam_max = GRATING_RANGES[disperser][filt]
    return (lam_max - lam_min) / LAMBDA_PER_ARCSEC

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
