"""NIRSpec MSA hand-picking planner — Bokeh server entry point.

Run:  bokeh serve app/ --show
"""
from __future__ import annotations

import base64
import sys
import tempfile
import traceback
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    CheckboxGroup,
    ColumnDataSource,
    Div,
    FileInput,
    HoverTool,
    Select,
    Slider,
    TapTool,
    TextInput,
)
from bokeh.plotting import figure

from app.catalog import Catalog, catalog_in_view, load_catalog
from app.coords import (
    MSA_V2_REF,
    MSA_V3_REF,
    V3_IDL_Y_ANGLE,
    rot_matrix,
    shutter_corners_v2v3,
    v2v3_to_radec,
)
from app.empt_io import (
    OpenShutter,
    Pointing,
    write_observed_targets_cat,
    write_pointing_summary_txt,
    write_shutter_mask_csv,
)
from app.image_io import LoadedImage, load_fits, load_jpg_with_sidecar, stretch_for_display
from app.msa import load_msa_grid, load_operability
from app.wavelengths import FILTER_BLUE_CUTOFF, GRATING_RANGES, cutoffs

# ---------------------------------------------------------------------------
# Singletons
# ---------------------------------------------------------------------------

V2_MSA, V3_MSA = load_msa_grid()           # (4, 171, 365)
OPERABLE, REASON = load_operability()      # (4, 171, 365) bool / int8

DISPERSERS = list(GRATING_RANGES.keys())   # PRISM, G140M, ...
FILTER_OPTIONS = ["CLEAR", "F070LP", "F100LP", "F170LP", "F290LP"]

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------

state: dict = {
    "image": None,
    "catalog": None,
    "tmp_sidecar_path": None,
    "open_shutters": {},  # (q,s,d) -> OpenShutter
    "history": [],        # stack of prior open_shutters snapshots (capped)
    "pa_v3": 0.0,
    "ra_deg": 0.0,
    "dec_deg": 0.0,
    "disperser": "PRISM",
    "filter": "CLEAR",
    "slitlet_height": 3,
    "snap_to_operable": True,
}


def _push_history() -> None:
    """Snapshot open_shutters for undo (cap at 50)."""
    state["history"].append(dict(state["open_shutters"]))
    if len(state["history"]) > 50:
        state["history"].pop(0)

# ---------------------------------------------------------------------------
# Bokeh widgets / glyphs
# ---------------------------------------------------------------------------

status = Div(text="Enter a file path below (or use the upload widgets for small files).", width=420)

# Path-based inputs are the primary way to load — no WebSocket size limit, no temp files.
fits_path_input = TextInput(title="FITS path (local)", value="", placeholder="/path/to/image.fits")
jpg_path_input = TextInput(title="JPG path (local)", value="", placeholder="/path/to/image.jpg")
sidecar_path_input = TextInput(title="Sidecar FITS path (WCS for JPG)", value="", placeholder="/path/to/wcs.fits")
catalog_path_input = TextInput(title="Catalog path (local)", value="", placeholder="/path/to/catalog.csv")

# Upload widgets work for small files but Bokeh's default WebSocket limit (~20 MB) will
# silently truncate larger ones. Start the server with --websocket-max-message-size if
# you want to use these for big files.
fits_input = FileInput(accept=".fits,.fit", title="FITS upload (small files)")
jpg_input = FileInput(accept=".jpg,.jpeg,.png", title="JPG upload (small files)")
sidecar_input = FileInput(accept=".fits", title="Sidecar FITS upload")
catalog_input = FileInput(accept=".csv,.cat,.txt,.fits", title="Catalog upload")

ra_input = TextInput(title="Pointing RA (deg)", value="")
dec_input = TextInput(title="Pointing Dec (deg)", value="")
pa_slider = Slider(title="APA_V3 (deg)", start=0.0, end=360.0, step=0.1, value=0.0)
pa_input = TextInput(title="APA_V3 (deg, exact)", value="0.0")

# Visibility window query (jwst_gtvt)
visibility_date_input = TextInput(
    title="Visibility date (YYYY-MM-DD)", value="",
    placeholder="leave blank for today",
)
visibility_btn = Button(label="Compute allowed APA_V3 (jwst_gtvt)", button_type="primary")
visibility_div = Div(text="<small>Allowed APA_V3 windows appear here.</small>", width=320)

disperser_select = Select(title="Disperser", options=DISPERSERS, value="PRISM")
filter_select = Select(title="Filter", options=FILTER_OPTIONS, value="CLEAR")

layers_box = CheckboxGroup(
    labels=["Show MSA outline", "Show operable shutters", "Apply operability", "Show targets"],
    active=[0, 1, 2, 3],
)
slitlet_select = Select(
    title="Slitlet height", options=["1", "3", "5"], value="3",
)

snap_box = CheckboxGroup(labels=["Snap target to nearest operable"], active=[0])

undo_btn = Button(label="Undo last", button_type="default")
clear_btn = Button(label="Clear open", button_type="warning")

export_dir_input = TextInput(title="Export dir", value=str(Path.cwd() / "exports"))
export_btn = Button(label="Export eMPT bundle", button_type="success")

# Glyph data sources
src_image = ColumnDataSource(data=dict(image=[], x=[], y=[], dw=[], dh=[]))
src_msa_outline = ColumnDataSource(data=dict(xs=[], ys=[]))
src_bg_shutters = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], reason=[])
)
src_open_shutters = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], target=[], lam_blue=[], lam_red=[])
)
src_targets = ColumnDataSource(data=dict(x=[], y=[], id=[], ra=[], dec=[], pr=[]))

fig = figure(
    width=900, height=900,
    tools="pan,wheel_zoom,box_zoom,reset,save,tap",
    match_aspect=True,
    output_backend="webgl",
    title="NIRSpec MSA planner",
    x_axis_label="pix x", y_axis_label="pix y",
)
fig.toolbar.active_scroll = fig.tools[1]  # wheel zoom enabled by default

img_glyph = fig.image_rgba(image="image", x="x", y="y", dw="dw", dh="dh", source=src_image)

msa_outline_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_msa_outline,
    line_color="dodgerblue", line_width=1.5, fill_alpha=0.0,
)
bg_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_bg_shutters,
    line_color="white", line_alpha=0.35, line_width=0.5,
    fill_color="white", fill_alpha=0.05,
)
open_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_open_shutters,
    line_color="#ff3333", line_width=1.5,
    fill_color="#ff8888", fill_alpha=0.35,
)
target_glyph = fig.scatter(
    x="x", y="y", source=src_targets,
    size=10, marker="circle", line_color="yellow", fill_alpha=0.0, line_width=1.5,
)

fig.add_tools(HoverTool(
    renderers=[bg_shutters_glyph],
    tooltips=[("Q,s,d", "@q,@s,@d"), ("op", "@reason")],
))
fig.add_tools(HoverTool(
    renderers=[open_shutters_glyph],
    tooltips=[("Q,s,d", "@q,@s,@d"),
              ("target", "@target"),
              ("λ blue/red", "@lam_blue{0.00} / @lam_red{0.00} μm")],
))
fig.add_tools(HoverTool(
    renderers=[target_glyph],
    tooltips=[("ID", "@id"), ("RA,Dec", "@ra{0.0000}, @dec{0.0000}"), ("Pr", "@pr")],
))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _set_status(msg: str, level: str = "info") -> None:
    color = {"info": "#222", "warn": "#a06000", "err": "#a00000", "ok": "#006020"}[level]
    status.text = f'<div style="color:{color}">{msg}</div>'


def _write_temp(b64_value: str, suffix: str) -> str:
    raw = base64.b64decode(b64_value)
    f = tempfile.NamedTemporaryFile(prefix="msa_planner_", suffix=suffix, delete=False)
    f.write(raw)
    f.close()
    return f.name


def _image_array_for_bokeh(img: LoadedImage, stretch: str = "asinh") -> np.ndarray:
    """Return a uint32 RGBA array oriented for Bokeh's image_rgba.

    Bokeh draws with image[0,0] at lower-left. FITS arrays are already that way.
    PIL JPEG arrays come in top-left origin; flip them vertically.
    """
    rgba = stretch_for_display(img.data, stretch=stretch)
    if img.mode == "jpg+sidecar":
        rgba = rgba[::-1]
    return rgba


def _world_to_pixel(coords: SkyCoord, wcs: WCS) -> tuple[np.ndarray, np.ndarray]:
    x, y = skycoord_to_pixel(coords, wcs, origin=0)
    return np.asarray(x), np.asarray(y)


def _pointing_skycoord() -> SkyCoord | None:
    try:
        ra = float(ra_input.value)
        dec = float(dec_input.value)
        return SkyCoord(ra, dec, unit=u.deg, frame="icrs")
    except (TypeError, ValueError):
        return None


def _shutter_polygons_in_view(
    fiducial: SkyCoord, pa_v3: float, wcs: WCS,
    apply_operability: bool,
    view_bbox_pix: tuple[float, float, float, float],
) -> dict:
    """Compute shutter polygons in image-pixel coords.

    Returns dict keyed by (q, s, d) -> {'xs': [...], 'ys': [...], 'reason': int}.
    Culls to shutters whose centers fall inside view_bbox_pix plus a 5″ margin.
    """
    # Convert *all* shutter centers to pixel once; cheap (250k coords).
    all_v2 = V2_MSA.reshape(-1)
    all_v3 = V3_MSA.reshape(-1)
    centers_v2v3 = np.stack([all_v2, all_v3], axis=1)
    sky_centers = v2v3_to_radec(fiducial, pa_v3, centers_v2v3)
    px, py = _world_to_pixel(sky_centers, wcs)

    xmin, xmax, ymin, ymax = view_bbox_pix
    margin = 50  # pixels
    in_view = (
        (px >= xmin - margin)
        & (px <= xmax + margin)
        & (py >= ymin - margin)
        & (py <= ymax + margin)
    )
    idx = np.where(in_view)[0]

    if apply_operability:
        flat_operable = OPERABLE.reshape(-1)
        idx = idx[flat_operable[idx]]

    if idx.size == 0:
        return {}

    # Build (M, 4, 2) corner array in V2/V3 for selected shutters
    sel_v2 = all_v2[idx]
    sel_v3 = all_v3[idx]
    # Stack: for each shutter, 4 corners
    M = idx.size
    corners_v2v3 = np.empty((M * 4, 2), dtype=np.float64)
    for k, (v2c, v3c) in enumerate(zip(sel_v2, sel_v3)):
        corners_v2v3[k * 4 : k * 4 + 4] = shutter_corners_v2v3(v2c, v3c)
    sky_corners = v2v3_to_radec(fiducial, pa_v3, corners_v2v3)
    cx, cy = _world_to_pixel(sky_corners, wcs)
    cx = cx.reshape(M, 4)
    cy = cy.reshape(M, 4)

    flat_reason = REASON.reshape(-1)
    # Recover (q, s, d) from flat index in (4, 171, 365)
    qs = idx // (171 * 365)
    rem = idx % (171 * 365)
    ss = rem // 365
    ds = rem % 365

    out = {}
    for k in range(M):
        out[(int(qs[k]) + 1, int(ss[k]) + 1, int(ds[k]) + 1)] = {
            "xs": cx[k].tolist(),
            "ys": cy[k].tolist(),
            "reason": int(flat_reason[idx[k]]),
        }
    return out


def _msa_outline_polygons(fiducial: SkyCoord, pa_v3: float, wcs: WCS) -> dict:
    """Trace each quadrant outline (using its 4 corner shutters)."""
    xs_all, ys_all = [], []
    for q in range(4):
        # 4 corners of the quadrant in V2/V3 (rows {0,170}, cols {0,364})
        rows = [0, 0, 170, 170]
        cols = [0, 364, 364, 0]
        v2c = np.array([V2_MSA[q, r, c] for r, c in zip(rows, cols)])
        v3c = np.array([V3_MSA[q, r, c] for r, c in zip(rows, cols)])
        corners = np.stack([v2c, v3c], axis=1)
        sky = v2v3_to_radec(fiducial, pa_v3, corners)
        x, y = _world_to_pixel(sky, wcs)
        xs_all.append([[x.tolist()]])
        ys_all.append([[y.tolist()]])
    return dict(xs=xs_all, ys=ys_all)


def _shutter_polys_to_cds(polys: dict) -> dict:
    """Convert shutter polygons dict to MultiPolygons CDS data."""
    xs, ys, qs, ss, ds, reasons = [], [], [], [], [], []
    for (q, s, d), p in polys.items():
        xs.append([[p["xs"]]])
        ys.append([[p["ys"]]])
        qs.append(q)
        ss.append(s)
        ds.append(d)
        reasons.append(p["reason"])
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds, reason=reasons)


def _open_shutters_cds_data() -> dict:
    """Recompute open-shutter polygons for the current pointing."""
    img = state["image"]
    fiducial = _pointing_skycoord()
    if img is None or fiducial is None or not state["open_shutters"]:
        return dict(xs=[], ys=[], q=[], s=[], d=[], target=[], lam_blue=[], lam_red=[])

    pa_v3 = state["pa_v3"]
    wcs = img.wcs
    xs, ys, qs, ss, ds, tgt, lam_b, lam_r = [], [], [], [], [], [], [], []
    for (q, s, d), sh in state["open_shutters"].items():
        v2c = V2_MSA[q - 1, s - 1, d - 1]
        v3c = V3_MSA[q - 1, s - 1, d - 1]
        corners_v2v3 = shutter_corners_v2v3(v2c, v3c)
        sky = v2v3_to_radec(fiducial, pa_v3, corners_v2v3)
        cx, cy = _world_to_pixel(sky, wcs)
        xs.append([[cx.tolist()]])
        ys.append([[cy.tolist()]])
        qs.append(q); ss.append(s); ds.append(d)
        tgt.append(str(sh.target_id) if sh.target_id is not None else "")
        cut = cutoffs(v2c, v3c, state["disperser"], state["filter"])
        lam_b.append(cut.get("lam_blue") if cut and cut.get("lam_blue") is not None else float("nan"))
        lam_r.append(cut.get("lam_red") if cut and cut.get("lam_red") is not None else float("nan"))
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds, target=tgt, lam_blue=lam_b, lam_red=lam_r)


# ---------------------------------------------------------------------------
# Refresh logic
# ---------------------------------------------------------------------------


def refresh_overlays() -> None:
    img: LoadedImage | None = state["image"]
    if img is None:
        _set_status("Load an image first.", "warn")
        return
    fiducial = _pointing_skycoord()
    if fiducial is None:
        _set_status("Enter a valid RA and Dec.", "warn")
        return

    pa_v3 = state["pa_v3"]
    wcs = img.wcs
    H, W = img.shape

    # MSA outline
    show_outline = 0 in layers_box.active
    if show_outline:
        src_msa_outline.data = _msa_outline_polygons(fiducial, pa_v3, wcs)
    else:
        src_msa_outline.data = dict(xs=[], ys=[])

    # Background shutters (operable filter optional)
    show_shutters = 1 in layers_box.active
    apply_op = 2 in layers_box.active
    if show_shutters:
        view_bbox = (0, W, 0, H)
        polys = _shutter_polygons_in_view(fiducial, pa_v3, wcs, apply_op, view_bbox)
        src_bg_shutters.data = _shutter_polys_to_cds(polys)
    else:
        src_bg_shutters.data = dict(xs=[], ys=[], q=[], s=[], d=[], reason=[])

    # Open shutters (always shown)
    src_open_shutters.data = _open_shutters_cds_data()

    # Targets
    show_targets = 3 in layers_box.active
    cat: Catalog | None = state["catalog"]
    if show_targets and cat is not None:
        coords = SkyCoord(cat.ra_deg, cat.dec_deg, unit=u.deg, frame="icrs")
        x, y = _world_to_pixel(coords, wcs)
        # Cull to image bounds
        mask = (x >= 0) & (x < W) & (y >= 0) & (y < H)
        src_targets.data = dict(
            x=x[mask].tolist(),
            y=y[mask].tolist(),
            id=[str(i) for i in np.asarray(cat.ids)[mask]],
            ra=cat.ra_deg[mask].tolist(),
            dec=cat.dec_deg[mask].tolist(),
            pr=cat.priority[mask].tolist(),
        )
    else:
        src_targets.data = dict(x=[], y=[], id=[], ra=[], dec=[], pr=[])

    n_op = len(state["open_shutters"])
    n_tgt_open = len({sh.target_id for sh in state["open_shutters"].values() if sh.target_id})
    _set_status(
        f"{n_op} open shutters covering {n_tgt_open} targets. "
        f"PA_V3 = {pa_v3:.2f}°.", "ok",
    )


def refresh_image_glyph() -> None:
    img = state["image"]
    if img is None:
        return
    arr = _image_array_for_bokeh(img)
    H, W = arr.shape
    src_image.data = dict(image=[arr], x=[0], y=[0], dw=[W], dh=[H])
    fig.x_range.start = 0
    fig.x_range.end = W
    fig.y_range.start = 0
    fig.y_range.end = H
    # Match figure aspect to image aspect so each FITS pixel is a square
    # screen pixel (no horizontal/vertical stretching). match_aspect=True
    # locks the data ratio through pan/zoom. Assumes square arcsec/pix in
    # the WCS, which is the case for typical JWST/HST drizzled mosaics.
    base = 900
    if W >= H:
        fig.width = base
        fig.height = max(300, int(base * H / W))
    else:
        fig.height = base
        fig.width = max(300, int(base * W / H))


# ---------------------------------------------------------------------------
# Callbacks: file inputs
# ---------------------------------------------------------------------------


def _set_image_and_recenter(img: LoadedImage, source_label: str) -> None:
    state["image"] = img
    H, W = img.shape[:2]
    try:
        center = img.wcs.pixel_to_world(W / 2, H / 2)
        ra_input.value = f"{center.ra.deg:.6f}"
        dec_input.value = f"{center.dec.deg:.6f}"
        state["ra_deg"] = center.ra.deg
        state["dec_deg"] = center.dec.deg
    except Exception:  # noqa: BLE001
        pass
    refresh_image_glyph()
    refresh_overlays()
    _set_status(f"Loaded {source_label} ({W}×{H}).", "ok")


def _load_fits_from_path(path: str) -> None:
    try:
        img = load_fits(path)
        _set_image_and_recenter(img, f"FITS {Path(path).name}")
    except Exception as e:  # noqa: BLE001
        _set_status(f"FITS load failed: {e}", "err")
        traceback.print_exc()


def _load_jpg_pair_from_paths(jpg_path: str, sidecar_path: str) -> None:
    try:
        img = load_jpg_with_sidecar(jpg_path, sidecar_path, max_dim=6000)
        _set_image_and_recenter(img, f"JPG+sidecar {Path(jpg_path).name}")
    except Exception as e:  # noqa: BLE001
        _set_status(f"JPG+sidecar load failed: {e}", "err")
        traceback.print_exc()


def _load_catalog_from_path(path: str) -> None:
    try:
        cat = load_catalog(path)
        state["catalog"] = cat
        refresh_overlays()
        _set_status(f"Catalog loaded: {len(cat.ra_deg)} targets from {Path(path).name}.", "ok")
    except Exception as e:  # noqa: BLE001
        _set_status(f"Catalog load failed: {e}", "err")
        traceback.print_exc()


# Path-based callbacks (primary input for a local tool)
def on_fits_path(attr, old, new):
    p = fits_path_input.value.strip()
    if p and Path(p).exists():
        _load_fits_from_path(p)
    elif p:
        _set_status(f"FITS path not found: {p}", "err")


def on_sidecar_path(attr, old, new):
    p = sidecar_path_input.value.strip()
    if not p:
        return
    if not Path(p).exists():
        _set_status(f"Sidecar path not found: {p}", "err")
        return
    state["tmp_sidecar_path"] = p
    jpg_p = jpg_path_input.value.strip()
    if jpg_p and Path(jpg_p).exists():
        _load_jpg_pair_from_paths(jpg_p, p)
    else:
        _set_status("Sidecar set. Now enter a JPG path.", "info")


def on_jpg_path(attr, old, new):
    jpg_p = jpg_path_input.value.strip()
    side_p = state.get("tmp_sidecar_path") or sidecar_path_input.value.strip()
    if not jpg_p:
        return
    if not Path(jpg_p).exists():
        _set_status(f"JPG path not found: {jpg_p}", "err")
        return
    if not side_p or not Path(side_p).exists():
        _set_status("Enter a sidecar FITS path first.", "warn")
        return
    _load_jpg_pair_from_paths(jpg_p, side_p)


def on_catalog_path(attr, old, new):
    p = catalog_path_input.value.strip()
    if p and Path(p).exists():
        _load_catalog_from_path(p)
    elif p:
        _set_status(f"Catalog path not found: {p}", "err")


# Upload-based callbacks (for small files only — Bokeh WebSocket limit applies)
def on_fits(attr, old, new):
    if not fits_input.value:
        return
    try:
        path = _write_temp(fits_input.value, suffix=".fits")
        _load_fits_from_path(path)
    except Exception as e:  # noqa: BLE001
        _set_status(f"FITS upload failed: {e}", "err")


def _try_load_jpg_upload():
    if not jpg_input.value or state["tmp_sidecar_path"] is None:
        return
    try:
        jpg_path = _write_temp(jpg_input.value, suffix=Path(jpg_input.filename).suffix or ".jpg")
        _load_jpg_pair_from_paths(jpg_path, state["tmp_sidecar_path"])
    except Exception as e:  # noqa: BLE001
        _set_status(f"JPG upload failed: {e}", "err")


def on_jpg(attr, old, new):
    _try_load_jpg_upload()


def on_sidecar(attr, old, new):
    if not sidecar_input.value:
        return
    try:
        state["tmp_sidecar_path"] = _write_temp(sidecar_input.value, suffix=".fits")
        _set_status("Sidecar uploaded. Now upload the JPG.", "info")
        _try_load_jpg_upload()
    except Exception as e:  # noqa: BLE001
        _set_status(f"Sidecar upload failed: {e}", "err")


def on_catalog(attr, old, new):
    if not catalog_input.value:
        return
    try:
        suffix = Path(catalog_input.filename).suffix or ".csv"
        path = _write_temp(catalog_input.value, suffix=suffix)
        _load_catalog_from_path(path)
    except Exception as e:  # noqa: BLE001
        _set_status(f"Catalog upload failed: {e}", "err")


# ---------------------------------------------------------------------------
# Callbacks: pointing & display
# ---------------------------------------------------------------------------


def on_pointing(attr, old, new):
    try:
        state["ra_deg"] = float(ra_input.value)
        state["dec_deg"] = float(dec_input.value)
    except (TypeError, ValueError):
        return
    refresh_overlays()


def on_pa_slider(attr, old, new):
    state["pa_v3"] = float(pa_slider.value)
    pa_input.value = f"{state['pa_v3']:.2f}"
    refresh_overlays()


def on_pa_text(attr, old, new):
    try:
        v = float(pa_input.value) % 360.0
        state["pa_v3"] = v
        pa_slider.value = v
    except (TypeError, ValueError):
        return
    refresh_overlays()


def on_layers(attr, old, new):
    refresh_overlays()


def on_disperser(attr, old, new):
    state["disperser"] = disperser_select.value
    refresh_overlays()


def on_filter(attr, old, new):
    state["filter"] = filter_select.value
    refresh_overlays()


def on_slitlet_height(attr, old, new):
    state["slitlet_height"] = int(slitlet_select.value)


def on_snap(attr, old, new):
    state["snap_to_operable"] = 0 in snap_box.active


# ---------------------------------------------------------------------------
# Hand-picking: tap callbacks and slitlet builder
# ---------------------------------------------------------------------------


def _sky_to_v2v3(sky: SkyCoord, fiducial: SkyCoord, pa_v3: float) -> tuple[float, float]:
    """Inverse of v2v3_to_radec: returns (V2, V3) in arcsec for a single sky point."""
    d_lon, d_lat = fiducial.spherical_offsets_to(sky)
    dx = d_lon.to_value(u.arcsec)
    dy = d_lat.to_value(u.arcsec)
    # Inverse of dot(offsets, rot_matrix(pa_v3)) is dot(offsets, rot_matrix(-pa_v3))
    rot = rot_matrix(-pa_v3)
    v2_off, v3_off = float(dx * rot[0, 0] + dy * rot[1, 0]), float(dx * rot[0, 1] + dy * rot[1, 1])
    return v2_off + MSA_V2_REF, v3_off + MSA_V3_REF


def _nearest_shutter(v2_target: float, v3_target: float,
                     require_operable: bool = True) -> tuple[int, int, int] | None:
    """Brute-force argmin over all shutters of squared V2/V3 distance."""
    dv2 = V2_MSA - v2_target
    dv3 = V3_MSA - v3_target
    d2 = dv2 * dv2 + dv3 * dv3
    if require_operable:
        d2 = np.where(OPERABLE, d2, np.inf)
    idx = int(np.argmin(d2))
    if not np.isfinite(d2.flat[idx]):
        return None
    q = idx // (171 * 365)
    rem = idx % (171 * 365)
    s = rem // 365
    d = rem % 365
    return (q + 1, s + 1, d + 1)


def _add_slitlet(q: int, s_center: int, d: int, target_id: str | None) -> int:
    """Add a slitlet of state['slitlet_height'] shutters centered on (q,s_center,d).

    Returns the number of shutters added (operable ones; skips failed)."""
    h = state["slitlet_height"]
    half = h // 2
    added = 0
    for offset in range(-half, half + 1):
        s = s_center + offset
        if not (1 <= s <= 171):
            continue
        if not OPERABLE[q - 1, s - 1, d - 1]:
            continue
        role = "target" if offset == 0 else "sky"
        state["open_shutters"][(q, s, d)] = OpenShutter(
            q=q, s=s, d=d, target_id=target_id, role=role,
        )
        added += 1
    return added


def on_target_tap(attr, old, new):
    if not new:
        return
    idx = int(new[0])
    src_targets.selected.indices = []  # clear selection so the next tap re-fires
    if idx >= len(src_targets.data["x"]):
        return
    img = state["image"]
    fiducial = _pointing_skycoord()
    if img is None or fiducial is None:
        return
    px = src_targets.data["x"][idx]
    py = src_targets.data["y"][idx]
    tgt_id = src_targets.data["id"][idx]
    sky_target = img.wcs.pixel_to_world(px, py)
    v2, v3 = _sky_to_v2v3(sky_target, fiducial, state["pa_v3"])
    nearest = _nearest_shutter(v2, v3, require_operable=state["snap_to_operable"])
    if nearest is None:
        _set_status(f"Target {tgt_id}: no operable shutter nearby.", "warn")
        return
    q, s, d = nearest
    _push_history()
    n = _add_slitlet(q, s, d, target_id=str(tgt_id))
    _set_status(f"Target {tgt_id} → slitlet ({q},{s},{d}), {n} shutters opened.", "ok")
    refresh_overlays()


def on_bg_shutter_tap(attr, old, new):
    if not new:
        return
    idx = int(new[0])
    src_bg_shutters.selected.indices = []
    if idx >= len(src_bg_shutters.data["q"]):
        return
    q = int(src_bg_shutters.data["q"][idx])
    s = int(src_bg_shutters.data["s"][idx])
    d = int(src_bg_shutters.data["d"][idx])
    key = (q, s, d)
    _push_history()
    if key in state["open_shutters"]:
        del state["open_shutters"][key]
        _set_status(f"Closed shutter ({q},{s},{d}).", "ok")
    else:
        state["open_shutters"][key] = OpenShutter(q=q, s=s, d=d, role="manual")
        _set_status(f"Opened shutter ({q},{s},{d}) manually.", "ok")
    refresh_overlays()


def on_open_shutter_tap(attr, old, new):
    if not new:
        return
    idx = int(new[0])
    src_open_shutters.selected.indices = []
    if idx >= len(src_open_shutters.data["q"]):
        return
    q = int(src_open_shutters.data["q"][idx])
    s = int(src_open_shutters.data["s"][idx])
    d = int(src_open_shutters.data["d"][idx])
    key = (q, s, d)
    _push_history()
    if key in state["open_shutters"]:
        sh = state["open_shutters"].pop(key)
        # Remove slitlet siblings (same target_id, same q, d, contiguous s)
        if sh.target_id:
            for k in list(state["open_shutters"].keys()):
                other = state["open_shutters"][k]
                if (other.target_id == sh.target_id
                        and other.q == q and other.d == d
                        and abs(other.s - s) <= state["slitlet_height"] // 2):
                    del state["open_shutters"][k]
        _set_status(f"Removed shutter ({q},{s},{d}) and slitlet siblings.", "ok")
    refresh_overlays()


def on_undo():
    if not state["history"]:
        _set_status("Nothing to undo.", "warn")
        return
    state["open_shutters"] = state["history"].pop()
    refresh_overlays()
    _set_status("Undid last action.", "ok")


def on_clear():
    if not state["open_shutters"]:
        return
    _push_history()
    state["open_shutters"] = {}
    refresh_overlays()
    _set_status("Cleared all open shutters.", "ok")


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


def on_export():
    if not state["open_shutters"]:
        _set_status("No open shutters to export.", "warn")
        return
    fiducial = _pointing_skycoord()
    if fiducial is None:
        _set_status("Set a valid RA/Dec before exporting.", "err")
        return
    out_dir = Path(export_dir_input.value).expanduser()
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    out_dir = out_dir / f"empt_bundle_{stamp}"
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        pass

    # observed_targets: one row per unique target_id with an open shutter
    cat = state["catalog"]
    targets_rows = []
    seen = set()
    for sh in state["open_shutters"].values():
        if sh.target_id and sh.target_id not in seen:
            seen.add(sh.target_id)
            # Look up RA/Dec from catalog if possible
            ra_d, dec_d = float("nan"), float("nan")
            pr = 1
            if cat is not None:
                ids_str = [str(i) for i in cat.ids]
                if str(sh.target_id) in ids_str:
                    k = ids_str.index(str(sh.target_id))
                    ra_d = float(cat.ra_deg[k])
                    dec_d = float(cat.dec_deg[k])
                    if np.isfinite(cat.priority[k]):
                        pr = int(cat.priority[k])
            try:
                target_no = int(sh.target_id)
            except (ValueError, TypeError):
                target_no = len(targets_rows) + 1
            targets_rows.append({
                "No_cat": target_no,
                "Pr": pr,
                "ra_deg": ra_d,
                "dec_deg": dec_d,
            })

    pa_v3 = state["pa_v3"]
    # PA_V3 - PA_AP = -V3IdlYAngle (mod 360); for NRS_FULL_MSA V3IdlYAngle ~ 138.5746°.
    pa_ap = (pa_v3 + V3_IDL_Y_ANGLE) % 360.0
    pointing = Pointing(
        ra_deg=float(ra_input.value),
        dec_deg=float(dec_input.value),
        apa_v3_deg=pa_v3,
        pa_ap_deg=pa_ap,
    )
    try:
        write_observed_targets_cat(str(out_dir / "observed_targets.cat"), targets_rows)
        write_pointing_summary_txt(
            str(out_dir / "pointing_summary.txt"),
            pointing, state["disperser"], state["filter"],
            n_targets_total=(len(cat.ra_deg) if cat is not None else 0),
            n_targets_accepted=len(targets_rows),
        )
        open_list = list(state["open_shutters"].values())
        write_shutter_mask_csv(
            str(out_dir / "shutter_mask.csv"),
            open_list, OPERABLE, REASON,
        )
        _set_status(
            f"Exported eMPT bundle to {out_dir} "
            f"({len(targets_rows)} targets, {len(open_list)} open shutters).",
            "ok",
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"Export failed: {e}", "err")
        traceback.print_exc()


# Tap wiring
src_targets.selected.on_change("indices", on_target_tap)
src_bg_shutters.selected.on_change("indices", on_bg_shutter_tap)
src_open_shutters.selected.on_change("indices", on_open_shutter_tap)
undo_btn.on_click(on_undo)
clear_btn.on_click(on_clear)
export_btn.on_click(on_export)
snap_box.on_change("active", on_snap)


# ---------------------------------------------------------------------------
# Visibility / APA_V3 constraints (jwst_gtvt)
# ---------------------------------------------------------------------------


def on_visibility():
    fiducial = _pointing_skycoord()
    if fiducial is None:
        _set_status("Set RA/Dec before computing visibility.", "warn")
        return
    try:
        from jwst_gtvt.jwst_tvt import Ephemeris
    except ImportError:
        _set_status("jwst_gtvt not installed (pip install jwst_gtvt).", "err")
        return

    try:
        _set_status("Querying jwst_gtvt ephemeris… (one-time, ~5 s)", "info")
        eph = Ephemeris()
        df = eph.get_fixed_target_positions(
            f"{fiducial.ra.deg:.6f}", f"{fiducial.dec.deg:.6f}"
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"jwst_gtvt query failed: {e}", "err")
        traceback.print_exc()
        return

    obs = df[df["in_FOR"]].copy()
    if len(obs) == 0:
        visibility_div.text = (
            "<small>Target never observable by JWST (likely too close to the Sun "
            "across the ephemeris range).</small>"
        )
        return

    # Find the requested-date row (or today)
    date_str = visibility_date_input.value.strip()
    if not date_str:
        date_str = datetime.utcnow().strftime("%Y-%m-%d")
    # Match by 'Calendar Date (TDB)' substring; gtvt format e.g. "A.D. 2026-Sep-18 00:00:00.0000"
    yyyy, mm, dd = date_str.split("-")
    month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][int(mm) - 1]
    needle = f"{yyyy}-{month_abbr}-{int(dd):02d}"
    matches = obs[obs["Calendar Date (TDB)"].str.contains(needle, na=False)]

    if len(matches) == 0:
        visibility_div.text = (
            f"<small><b>{date_str}: not observable.</b> "
            f"Observable windows in ephemeris range: "
            f"{len(obs)} days total.</small>"
        )
        return

    row = matches.iloc[0]
    nom = float(row["V3PA_nominal_angle"])
    lo = float(row["V3PA_min_pa_angle"])
    hi = float(row["V3PA_max_pa_angle"])
    # Allowed range may wrap mod 360
    visibility_div.text = (
        f"<small><b>{date_str}: observable.</b><br>"
        f"V3PA nominal = <b>{nom:.2f}°</b><br>"
        f"V3PA allowed: <b>{lo:.2f}° – {hi:.2f}°</b><br>"
        f"(window width {hi - lo:.1f}°)</small>"
    )
    # Snap slider into the allowed range and recolor band
    pa_slider.value = nom
    _set_status(
        f"jwst_gtvt: {date_str} V3PA ∈ [{lo:.1f}, {hi:.1f}]°, nominal {nom:.1f}°.",
        "ok",
    )


visibility_btn.on_click(on_visibility)


# ---------------------------------------------------------------------------
# Wiring
# ---------------------------------------------------------------------------

fits_input.on_change("value", on_fits)
jpg_input.on_change("value", on_jpg)
sidecar_input.on_change("value", on_sidecar)
catalog_input.on_change("value", on_catalog)
fits_path_input.on_change("value", on_fits_path)
jpg_path_input.on_change("value", on_jpg_path)
sidecar_path_input.on_change("value", on_sidecar_path)
catalog_path_input.on_change("value", on_catalog_path)
ra_input.on_change("value", on_pointing)
dec_input.on_change("value", on_pointing)
pa_slider.on_change("value", on_pa_slider)
pa_input.on_change("value", on_pa_text)
layers_box.on_change("active", on_layers)
disperser_select.on_change("value", on_disperser)
filter_select.on_change("value", on_filter)
slitlet_select.on_change("value", on_slitlet_height)

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

sidebar = column(
    Div(text="<h3>Image</h3><b>Paste a local path</b> (works for any size):"),
    fits_path_input,
    Div(text="<i>or</i> JPG + sidecar FITS by path:"),
    sidecar_path_input,
    jpg_path_input,
    Div(text="<small>Upload widgets below — small files only "
             "(Bokeh WS limit ~20 MB unless server started with "
             "<code>--websocket-max-message-size</code>):</small>"),
    fits_input,
    sidecar_input,
    jpg_input,
    Div(text="<h3>Catalog</h3>Path:"),
    catalog_path_input,
    catalog_input,
    Div(text="<h3>Pointing</h3>"),
    ra_input, dec_input,
    pa_slider, pa_input,
    visibility_date_input, visibility_btn, visibility_div,
    Div(text="<h3>Instrument</h3>"),
    disperser_select, filter_select,
    Div(text="<h3>Display</h3>"),
    layers_box,
    slitlet_select,
    snap_box,
    row(undo_btn, clear_btn),
    Div(text="<h3>Export</h3>"),
    export_dir_input, export_btn,
    status,
    width=340,
)

curdoc().add_root(row(sidebar, fig))
curdoc().title = "NIRSpec MSA planner"
