"""End-to-end smoke test on the two example fields, bypassing Bokeh widgets."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from vmpt.coords import shutter_corners_v2v3, v2v3_to_radec  # noqa: E402
from vmpt.empt_io import (  # noqa: E402
    OpenShutter,
    Pointing,
    write_observed_targets_cat,
    write_pointing_summary_txt,
    write_shutter_mask_csv,
)
from vmpt.image_io import load_fits, load_jpg_with_sidecar  # noqa: E402
from vmpt.msa import load_msa_grid, load_operability  # noqa: E402
from vmpt.wavelengths import cutoffs  # noqa: E402

EX_A370 = _ROOT / "example_a370" / "a370_f182m_f200w_f210m.fits"
EX_R0600_JPG = _ROOT / "example_r0600" / "JWST_F090W_F200W_F444W.jpg"
EX_R0600_WCS = _ROOT / "example_r0600" / "wcs.fits"


def _msa_outline_pixels(img, pointing: SkyCoord, pa_v3: float):
    """Return per-quadrant outline polygons in image-pixel coords."""
    v2_msa, v3_msa = load_msa_grid()
    polys = []
    for q in range(4):
        rows = [0, 0, 170, 170]
        cols = [0, 364, 364, 0]
        v2c = np.array([v2_msa[q, r, c] for r, c in zip(rows, cols)])
        v3c = np.array([v3_msa[q, r, c] for r, c in zip(rows, cols)])
        corners = np.stack([v2c, v3c], axis=1)
        sky = v2v3_to_radec(pointing, pa_v3, corners)
        x, y = skycoord_to_pixel(sky, img.wcs, origin=0)
        polys.append((np.asarray(x), np.asarray(y)))
    return polys


@pytest.mark.skipif(not EX_A370.exists(), reason="example_a370 missing")
def test_a370_msa_overlay_inside_image():
    img = load_fits(str(EX_A370))
    H, W = img.shape
    # Point the MSA at the image center
    center = img.wcs.pixel_to_world(W / 2, H / 2)
    polys = _msa_outline_pixels(img, center, pa_v3=0.0)

    # At least one quadrant's centroid should fall inside the image
    inside = 0
    for x, y in polys:
        cx, cy = float(x.mean()), float(y.mean())
        if 0 <= cx < W and 0 <= cy < H:
            inside += 1
    assert inside >= 1, "no MSA quadrant lands on the a370 cutout"


@pytest.mark.skipif(not (EX_R0600_JPG.exists() and EX_R0600_WCS.exists()), reason="example_r0600 missing")
def test_r0600_jpg_sidecar_overlay():
    img = load_jpg_with_sidecar(str(EX_R0600_JPG), str(EX_R0600_WCS), max_dim=4000)
    H, W = img.shape[:2]
    center = img.wcs.pixel_to_world(W / 2, H / 2)
    polys = _msa_outline_pixels(img, center, pa_v3=137.0)
    # MSA should cover ~3.5×3.5 arcmin → at 30 mas/pix downsampled = ~1700 pix.
    # All 4 quadrants should at least partially overlap a 4000² image.
    inside = sum(
        1 for x, y in polys
        if (np.any((x >= 0) & (x < W)) and np.any((y >= 0) & (y < H)))
    )
    assert inside == 4, f"only {inside}/4 quadrants overlap r0600 image"


@pytest.mark.skipif(not EX_A370.exists(), reason="example_a370 missing")
def test_export_bundle_round_trips(tmp_path):
    img = load_fits(str(EX_A370))
    operable, reason = load_operability()
    # Build a tiny mock open-shutter set
    open_shutters = [
        OpenShutter(q=2, s=86, d=200, target_id=123, role="target"),
        OpenShutter(q=2, s=87, d=200, target_id=123, role="sky"),
        OpenShutter(q=2, s=85, d=200, target_id=123, role="sky"),
    ]
    pointing = Pointing(ra_deg=40.0, dec_deg=-1.6, apa_v3_deg=120.0, pa_ap_deg=258.5)

    cat_path = tmp_path / "observed_targets.cat"
    summ_path = tmp_path / "pointing_summary.txt"
    csv_path = tmp_path / "shutter_mask.csv"

    write_observed_targets_cat(
        str(cat_path),
        [{"No_cat": 123, "Pr": 1, "ra_deg": 40.001, "dec_deg": -1.599}],
    )
    write_pointing_summary_txt(
        str(summ_path), pointing, "G395M", "F290LP",
        n_targets_total=1, n_targets_accepted=1,
    )
    write_shutter_mask_csv(str(csv_path), open_shutters, operable, reason)

    # All three files exist and have plausible sizes
    assert cat_path.exists() and cat_path.stat().st_size > 0
    assert summ_path.exists() and summ_path.stat().st_size > 0
    assert csv_path.exists() and csv_path.stat().st_size > 1000  # 731 lines × 683 chars ≈ 500 KB


def test_wavelengths_for_arbitrary_shutter():
    v2_msa, v3_msa = load_msa_grid()
    # PRISM/CLEAR: the on-detector spectrum nearly always spans the
    # full filter range (msaviz reports blue_edge 5-95 % ≈ [0.600,
    # 0.606] μm and red_edge ≈ [5.29, 5.30] μm across all shutters),
    # because PRISM dispersion compresses the red and expands the
    # blue. So endpoints should be essentially identical across
    # quadrants — what shifts per-shutter is the *gap* location
    # within the spectrum.
    c_q2 = cutoffs(v2_msa[1, 85, 200], v3_msa[1, 85, 200], "PRISM", "CLEAR")
    c_q4 = cutoffs(v2_msa[3, 85, 200], v3_msa[3, 85, 200], "PRISM", "CLEAR")
    assert c_q2["lam_blue"] == c_q4["lam_blue"]
    assert c_q2["lam_red"] == c_q4["lam_red"]
    # G395M/F290LP is filter-limited on blue (always 2.87) but
    # lam_red shifts — the gratings DO follow a roughly linear V2
    # dispersion, so different-V2 shutters see different red edges.
    g_q2 = cutoffs(v2_msa[1, 85, 200], v3_msa[1, 85, 200], "G395M", "F290LP")
    g_q4 = cutoffs(v2_msa[3, 85, 200], v3_msa[3, 85, 200], "G395M", "F290LP")
    assert g_q2["lam_blue"] == g_q4["lam_blue"] == 2.87
    assert abs(g_q2["lam_red"] - g_q4["lam_red"]) > 0.1
