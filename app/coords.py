"""V2/V3 <-> RA/Dec coordinate transforms ported from footprint_emerald.ipynb."""

import numpy as np
import pysiaf
from astropy import units as u
from astropy.coordinates import SkyCoord

_siaf = pysiaf.Siaf("NIRSpec")
_msa_ap = _siaf["NRS_FULL_MSA"]

MSA_V2_REF: float = float(_msa_ap.V2Ref)
MSA_V3_REF: float = float(_msa_ap.V3Ref)


def rot_matrix(rotation: float = 30.0) -> np.ndarray:
    theta = np.radians(rotation)
    c, s = np.cos(theta), np.sin(theta)
    return np.array(((c, -s), (s, c)))


def shutter_corners_v2v3(v2c: float, v3c: float, w: float = 0.20, h: float = 0.46) -> np.ndarray:
    # 138.5 deg is the MSA tilt within the V2/V3 frame.
    return np.array([v2c, v3c]) + np.dot(
        np.array([[-w / 2, w / 2, w / 2, -w / 2],
                  [-h / 2, -h / 2, h / 2, h / 2]]).T,
        rot_matrix(138.5),
    )


def v2v3_to_radec(coord_c: SkyCoord, pa_v3: float, corners_v2v3: np.ndarray) -> SkyCoord:
    offsets = corners_v2v3 - np.array([MSA_V2_REF, MSA_V3_REF])
    offsets = np.dot(offsets, rot_matrix(pa_v3))
    return coord_c.spherical_offsets_by(
        offsets.T[0] * u.arcsec, offsets.T[1] * u.arcsec
    )
