"""vMPT — visual MSA Planning Tool. Bokeh server entry point.

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

from bokeh.events import DoubleTap, RangesUpdate, Tap
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    CheckboxGroup,
    ColumnDataSource,
    CustomJSTickFormatter,
    Div,
    FileInput,
    HoverTool,
    Select,
    Slider,
    TextInput,
    WheelZoomTool,
)
from bokeh.plotting import figure

from app.catalog import Catalog, catalog_in_view, load_catalog
from app.coords import (
    MSA_V2_REF,
    MSA_V3_REF,
    V3_IDL_Y_ANGLE,
    fixed_slit_corners_v2v3,
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
from app.session_io import Session, export_session_json, import_session_json
from app.wavelengths import (
    FILTER_BLUE_CUTOFF,
    GRATING_RANGES,
    cutoffs,
    v2_overlap_distance,
)

# ---------------------------------------------------------------------------
# Singletons
# ---------------------------------------------------------------------------

V2_MSA, V3_MSA = load_msa_grid()           # (4, 171, 365)
OPERABLE, REASON = load_operability()      # (4, 171, 365) bool / int8

DISPERSERS = list(GRATING_RANGES.keys())   # PRISM, G140M, ...
FILTER_OPTIONS = ["CLEAR", "F070LP", "F100LP", "F170LP", "F290LP"]
# Canonical disperser/filter combinations available in NIRSpec MOS. The
# combined dropdown drives the wavelength tooltip; values match
# observation-mode names used by APT.
DISPERSER_FILTER_COMBOS: list[tuple[str, str]] = [
    ("PRISM",  "CLEAR"),
    ("G140M",  "F070LP"),
    ("G140M",  "F100LP"),
    ("G235M",  "F170LP"),
    ("G395M",  "F290LP"),
    ("G140H",  "F070LP"),
    ("G140H",  "F100LP"),
    ("G235H",  "F170LP"),
    ("G395H",  "F290LP"),
]
DISPERSER_FILTER_LABELS = [f"{d} / {f}" for d, f in DISPERSER_FILTER_COMBOS]

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------

state: dict = {
    "image": None,
    "catalog": None,
    "tmp_sidecar_path": None,
    "open_shutters": {},  # (q,s,d) -> OpenShutter
    "highlighted": set(), # set of (q,s,d) tuples — visual flag, not exported
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

# Catalog filters — hide-able. Numeric thresholds; leave blank/empty to skip.
catalog_priority_input = TextInput(
    title="Show priority class ≤ (blank = all)", value="", placeholder="e.g. 3",
)
catalog_mag_input = TextInput(
    title="Show mag ≤ (blank = all)", value="", placeholder="e.g. 28",
)

ra_input = TextInput(title="Pointing RA (deg)", value="")
dec_input = TextInput(title="Pointing Dec (deg)", value="")

# V3 PA = position angle of the JWST V3 axis on sky. This is what drives the
# V2/V3 -> RA/Dec math. APT/MPT's "NIRSpec PA" is the *aperture* PA (APA),
# which differs by the V3IdlYAngle of NRS_FULL_MSA (~138.57 deg). We show
# both, synchronized.
v3pa_slider = Slider(title="V3 PA (deg)", start=0.0, end=360.0, step=0.1, value=0.0)
v3pa_input = TextInput(title="V3 PA (deg, exact)", value="0.0")
apa_input = TextInput(
    title=f"NIRSpec APA (deg) — V3PA + {V3_IDL_Y_ANGLE:.3f}°",
    value=f"{V3_IDL_Y_ANGLE % 360.0:.2f}",
)
pa_help_div = Div(text=(
    "<small><b>V3 PA</b>: JWST V3 axis PA on sky (drives the overlay). "
    f"<b>NIRSpec APA</b>: aperture PA of NRS_FULL_MSA = V3PA + {V3_IDL_Y_ANGLE:.2f}° "
    "(mod 360). APT/MPT calls this NIRSpec's 'Aperture PA'. "
    "<a href='https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/"
    "jwst-position-angles-ranges-and-offsets' target='_blank'>JDox reference</a>.</small>"
), width=320)

# Visibility window query (jwst_gtvt)
visibility_date_input = TextInput(
    title="Visibility date (YYYY-MM-DD)", value="",
    placeholder="leave blank for today",
)
visibility_btn = Button(label="Compute allowed V3 PA (jwst_gtvt)", button_type="primary")
visibility_div = Div(text="<small>Allowed V3 PA windows appear here.</small>", width=320)

# Loading banner — prominent indicator at the top, hidden by default.
loading_banner = Div(
    text="", width=320, height=30, visible=False,
    styles=dict(background="#fff3cd", color="#664d03",
                padding="6px 10px", border="1px solid #ffecb5",
                **{"border-radius": "4px", "font-weight": "bold"}),
)

# Quick-help panel rendered on the right side of the figure.
help_toggle_btn = Button(label="Hide help", button_type="default", width=110)
help_div = Div(
    width=320,
    styles=dict(
        background="#f8f9fa", color="#212529",
        padding="10px 14px", border="1px solid #dee2e6",
        **{"border-radius": "6px"},
    ),
    text="""
<h3 style='margin-top:0'>Quick guide</h3>
<b>1. Load image</b>
<ul style='margin:4px 0 8px 18px'>
  <li>Paste a local <b>FITS path</b>, or</li>
  <li>Paste a <b>JPG path</b> + <b>sidecar FITS</b> path (sidecar supplies WCS).</li>
</ul>
<b>2. Optional: target catalog</b>
<ul style='margin:4px 0 8px 18px'>
  <li>CSV/ASCII/FITS with at least <code>ID, RA, DEC</code>.</li>
</ul>
<b>3. Set pointing</b>
<ul style='margin:4px 0 8px 18px'>
  <li><b>V3 PA</b> drives the math; <b>NIRSpec APA</b> = V3PA + 138.575° (mod 360).</li>
  <li>Click <b>Compute allowed V3 PA</b> after entering a date to query jwst_gtvt.</li>
  <li><b>Shift+click</b> anywhere on the image to move the pointing center there. The <span style='color:lime;font-weight:bold'>lime cross</span> shows the current pointing.</li>
</ul>
<b>4. Hand-pick shutters</b>
<ul style='margin:4px 0 8px 18px'>
  <li><b>Click anywhere</b> on the image → snaps to the nearest shutter and toggles it open.
  Click near a yellow target → opens a 3-shutter slitlet on that target.</li>
  <li><b>Double-click</b> a shutter → toggles its <span style='color:cyan;font-weight:bold;background:#222;padding:0 4px'>cyan highlight</span>
  (visual flag, not exported).</li>
  <li><span style='color:#ff3333;font-weight:bold'>Open shutters</span> are red-filled.
  <span style='color:#ff2222;font-weight:bold'>Stuck-open</span> shutters: red edge.
  Failed-closed shutters are hidden.</li>
  <li><span style='color:orange;font-weight:bold'>Orange-tinted</span> shutters share an s-row with
  an open shutter — their spectra would overlap on the detector.</li>
  <li><span style='color:gold;font-weight:bold'>Gold</span> polygons are the 5 NIRSpec fixed slits.</li>
</ul>
<b>5. Save / share session</b>
<ul style='margin:4px 0 8px 18px'>
  <li><b>Save session</b> → writes a JSON snapshot of pointing, PA, all
  open shutters, highlighted shutters, and the image/catalog paths.</li>
  <li><b>Load session</b> → restores the snapshot. Hand the JSON to a
  collaborator to continue picking where you left off.</li>
</ul>
<b>6. Export to APT</b>
<ul style='margin:4px 0 8px 18px'>
  <li>Set the export dir, click <b>Export eMPT bundle</b>.</li>
  <li>Produces: <code>observed_targets.cat</code> (APT source catalog),
  <code>pointing_summary.txt</code> (PA values to copy-paste), and
  <code>shutter_mask.csv</code> (per-nod MSA mask, format byte-compatible with eMPT).</li>
</ul>
<b>Interactions</b>
<ul style='margin:4px 0 8px 18px'>
  <li><b>Wheel</b>: zoom both axes equally.</li>
  <li><b>Drag</b>: pan the view.</li>
  <li><b>Box zoom</b>: select box-zoom icon, then drag.</li>
  <li><b>Reset</b>: toolbar reset icon.</li>
</ul>
<p style='margin:4px 0'>See <code>README.md</code> in the project for the full reference.</p>
""",
)


def on_help_toggle():
    help_div.visible = not help_div.visible
    help_toggle_btn.label = "Hide help" if help_div.visible else "Show help"


help_toggle_btn.on_click(on_help_toggle)
help_panel = column(help_toggle_btn, help_div, width=340)

disperser_filter_select = Select(
    title="Disperser / Filter",
    options=DISPERSER_FILTER_LABELS,
    value="PRISM / CLEAR",
)

layers_box = CheckboxGroup(
    labels=["Show MSA outline", "Show operable shutters", "Show targets"],
    # Operable shutters OFF by default — drawing 180k+ polygons is slow.
    # User can toggle on once zoomed in to a region of interest.
    active=[0, 2],
)
MAX_OPERABLE_RENDER = 8000  # hard cap to keep redraws fast even when enabled
slitlet_select = Select(
    title="Slitlet height", options=["1", "3", "5"], value="3",
)

snap_box = CheckboxGroup(labels=["Snap target to nearest operable"], active=[0])

undo_btn = Button(label="Undo last", button_type="default")
clear_btn = Button(label="Clear open", button_type="warning")

export_dir_input = TextInput(title="Export dir", value=str(Path.cwd() / "exports"))
export_btn = Button(label="Export eMPT bundle", button_type="success")

# Session save/load: round-trips the full picking state for collaborators.
session_save_path_input = TextInput(
    title="Session save path",
    value=str(Path.cwd() / "exports" / "session.json"),
    placeholder="/path/to/session.json",
)
session_save_btn = Button(label="Save session", button_type="primary")
session_load_path_input = TextInput(
    title="Session load path",
    value="",
    placeholder="/path/to/session.json",
)
session_load_btn = Button(label="Load session", button_type="primary")

# Example data quick-load buttons (onboarding).
example_a370_btn = Button(label="Load Abell 370 example", button_type="default")
example_r0600_btn = Button(label="Load RXCJ0600 example", button_type="default")

# Statistics panel — always-visible counts.
stats_div = Div(
    width=320,
    styles=dict(
        background="#eef5ff", color="#1a3b66",
        padding="6px 10px", border="1px solid #b8d4ff",
        **{"border-radius": "4px", "font-size": "12px"},
    ),
    text="<b>Stats:</b> waiting for image…",
)

# Glyph data sources
src_image = ColumnDataSource(data=dict(image=[], x=[], y=[], dw=[], dh=[]))
src_msa_outline = ColumnDataSource(data=dict(xs=[], ys=[]))
src_bg_shutters = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[])
)
src_stuck_open = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[])
)
src_open_shutters = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], target=[], lam_blue=[], lam_red=[])
)
src_highlighted = ColumnDataSource(data=dict(xs=[], ys=[], q=[], s=[], d=[]))
src_spec_overlap = ColumnDataSource(data=dict(xs=[], ys=[], q=[], s=[], d=[]))
src_fixed_slits = ColumnDataSource(data=dict(xs=[], ys=[], name=[]))
src_targets = ColumnDataSource(data=dict(x=[], y=[], id=[], ra=[], dec=[], pr=[]))
src_pointing_handle = ColumnDataSource(data=dict(x=[], y=[]))

fig = figure(
    width=900, height=900,
    # No "tap" tool: it auto-selects clicked glyphs, which causes Bokeh's
    # default nonselection-rendering to fade every *other* open shutter to
    # 20% alpha — making them look "pale" after a click. We still receive
    # mouse clicks via fig.on_event(Tap, on_tap) without a TapTool present.
    tools="pan,box_zoom,reset,save",
    match_aspect=True,
    output_backend="webgl",
    title="vMPT — visual MSA Planning Tool",
    x_axis_label="RA (deg)", y_axis_label="Dec (deg)",
)
# Custom WheelZoomTool: scroll always zooms both axes equally, even when the
# cursor is over an axis. Bokeh's default (zoom_on_axis=True) lets scrolling
# on an axis zoom *only* that axis, which would distort the image aspect.
wheel_zoom = WheelZoomTool(dimensions="both", zoom_on_axis=False)
fig.add_tools(wheel_zoom)
fig.toolbar.active_scroll = wheel_zoom

img_glyph = fig.image_rgba(image="image", x="x", y="y", dw="dw", dh="dh", source=src_image)

# Bottom-to-top render order: image, operable shutters, stuck-open shutters,
# open shutters (user picks), highlighted shutters, MSA outline, fixed slits,
# targets, pointing handle.
bg_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_bg_shutters,
    line_color="white", line_alpha=0.35, line_width=0.5,
    fill_color="white", fill_alpha=0.05,
)
stuck_open_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_stuck_open,
    line_color="#ff2222", line_alpha=0.9, line_width=1.0,
    fill_color="#ff2222", fill_alpha=0.10,
)
# Spectral-overlap shutters: any operable shutter in the same s row of the
# same quadrant as a currently-open shutter. If opened, their dispersed
# spectra would overlap on the detector (MPT-style spectral conflict).
spec_overlap_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap,
    line_color="orange", line_alpha=0.9, line_width=1.0,
    fill_color="orange", fill_alpha=0.25,
)
open_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_open_shutters,
    line_color="#ff3333", line_width=1.5,
    fill_color="#ff8888", fill_alpha=0.35,
)
highlighted_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_highlighted,
    line_color="cyan", line_width=2.0, fill_alpha=0.0,
)
msa_outline_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_msa_outline,
    line_color="dodgerblue", line_width=1.5, fill_alpha=0.0,
)
fixed_slits_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_fixed_slits,
    line_color="gold", line_width=2.0, fill_color="gold", fill_alpha=0.18,
)
target_glyph = fig.scatter(
    x="x", y="y", source=src_targets,
    size=10, marker="circle", line_color="yellow", fill_alpha=0.0, line_width=1.5,
)
pointing_handle_glyph = fig.scatter(
    x="x", y="y", source=src_pointing_handle,
    size=18, marker="cross", line_color="lime", fill_color="lime",
    line_width=3, fill_alpha=0.6,
)
# The lime cross is a *visual* indicator of the pointing center only.
# To move the pointing center, shift-click anywhere on the image — see
# on_tap() below. (Bokeh's PointDrawTool drag interaction was unreliable
# across versions, so we use a click-with-modifier instead.)

fig.add_tools(HoverTool(
    renderers=[bg_shutters_glyph],
    tooltips=[("Q,s,d", "@q,@s,@d"), ("state", "operable")],
))
fig.add_tools(HoverTool(
    renderers=[stuck_open_glyph],
    tooltips=[("Q,s,d", "@q,@s,@d"), ("state", "STUCK OPEN")],
))
fig.add_tools(HoverTool(
    renderers=[open_shutters_glyph],
    tooltips=[("Q,s,d", "@q,@s,@d"),
              ("target", "@target"),
              ("λ blue/red", "@lam_blue{0.00} / @lam_red{0.00} μm")],
))
fig.add_tools(HoverTool(
    renderers=[fixed_slits_glyph],
    tooltips=[("slit", "@name")],
))
fig.add_tools(HoverTool(
    renderers=[target_glyph],
    tooltips=[("ID", "@id"), ("RA,Dec", "@ra{0.0000}, @dec{0.0000}"), ("Pr", "@pr")],
))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


_STATUS_AUTOCLEAR_GENERATION = [0]


def _set_status(msg: str, level: str = "info", clear_after: float = 6.0) -> None:
    """Set the sidebar status line. If `clear_after` is positive, fade the
    message after that many seconds so the UI doesn't accumulate stale
    messages. Pass 0 to keep the message indefinitely.
    """
    color = {"info": "#222", "warn": "#a06000", "err": "#a00000", "ok": "#006020"}[level]
    status.text = f'<div style="color:{color}">{msg}</div>'
    _STATUS_AUTOCLEAR_GENERATION[0] += 1
    gen = _STATUS_AUTOCLEAR_GENERATION[0]
    if clear_after > 0:
        def _maybe_clear():
            # Only clear if no newer status has been set since.
            if _STATUS_AUTOCLEAR_GENERATION[0] == gen:
                status.text = '<div style="color:#888">Ready.</div>'
        try:
            curdoc().add_timeout_callback(_maybe_clear, int(clear_after * 1000))
        except Exception:  # noqa: BLE001
            pass


_LOADING_GENERATION = [0]


def _show_loading(msg: str) -> None:
    """Show the yellow loading banner. Pair with _hide_loading().

    Schedules a 60-second safety timeout so the banner never gets stuck
    if a callback path fails to call _hide_loading.
    """
    loading_banner.text = f"⏳ {msg}"
    loading_banner.visible = True
    _LOADING_GENERATION[0] += 1
    gen = _LOADING_GENERATION[0]
    def _safety_hide():
        if _LOADING_GENERATION[0] == gen and loading_banner.visible:
            loading_banner.visible = False
            loading_banner.text = ""
    try:
        curdoc().add_timeout_callback(_safety_hide, 60_000)
    except Exception:  # noqa: BLE001
        pass


def _hide_loading() -> None:
    _LOADING_GENERATION[0] += 1  # cancel any pending safety timer
    loading_banner.visible = False
    loading_banner.text = ""


def _deferred(fn, *args, **kwargs):
    """Run fn on the next Bokeh document tick so any prior UI state changes
    (e.g. _show_loading) flush to the browser first."""
    curdoc().add_next_tick_callback(lambda: fn(*args, **kwargs))


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


# --- Vectorized geometry cache --------------------------------------------
# Precomputed once at module load. Shutter centers in V2/V3 offset form
# (relative to MSA_V2_REF, MSA_V3_REF). A single corner template (4 x 2) holds
# the shutter-corner offsets in the V2/V3 frame (already rotated by the 138.5°
# MSA tilt — same for every shutter), so per-shutter corners are just
# (v2_off, v3_off) + template.
_V2_OFFSETS_ALL = V2_MSA.reshape(-1) - MSA_V2_REF        # (249660,)
_V3_OFFSETS_ALL = V3_MSA.reshape(-1) - MSA_V3_REF        # (249660,)
_SHUTTER_CORNER_TEMPLATE = shutter_corners_v2v3(0.0, 0.0)  # (4, 2) in V2/V3


def _compute_wcs_jacobian(wcs: WCS, fx: float, fy: float) -> np.ndarray:
    """Return the inverse Jacobian mapping (Δeast_arcsec, Δnorth_arcsec) →
    (Δpix_x, Δpix_y) evaluated at pixel (fx, fy).

    NOTE: `c0.spherical_offsets_to(c1)` returns the proper tangent-plane
    (east, north) offsets — i.e. the east offset is already cos(dec)-
    corrected. Using `(c1.ra - c0.ra)` directly gives the RA *angular*
    difference, which is 1/cos(dec) too large near the poles and produces
    visibly-mis-placed overlays at Dec |~| 20° or more (this was the cause
    of the r0600 "click opens wrong shutter" bug).
    """
    c0 = wcs.pixel_to_world(fx, fy)
    c1 = wcs.pixel_to_world(fx + 1, fy)
    c2 = wcs.pixel_to_world(fx, fy + 1)
    d_e_x, d_n_x = c0.spherical_offsets_to(c1)
    d_e_y, d_n_y = c0.spherical_offsets_to(c2)
    J = np.array([
        [d_e_x.to(u.arcsec).value, d_e_y.to(u.arcsec).value],
        [d_n_x.to(u.arcsec).value, d_n_y.to(u.arcsec).value],
    ])
    return np.linalg.inv(J)


def _project_v2v3_offsets_to_pixel(
    v2_offsets: np.ndarray,
    v3_offsets: np.ndarray,
    pa_v3: float,
    fid_pix: tuple[float, float],
    jinv: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized: rotate V2/V3 offsets by PA, then map to pixel via jinv."""
    R = rot_matrix(pa_v3)
    offsets = np.stack([v2_offsets, v3_offsets], axis=-1)  # (..., 2)
    sky = offsets @ R                                       # (..., 2) lon/lat arcsec
    pix = sky @ jinv.T                                      # (..., 2) pixel offsets
    return fid_pix[0] + pix[..., 0], fid_pix[1] + pix[..., 1]


def _pointing_skycoord() -> SkyCoord | None:
    try:
        ra = float(ra_input.value)
        dec = float(dec_input.value)
        return SkyCoord(ra, dec, unit=u.deg, frame="icrs")
    except (TypeError, ValueError):
        return None


_FLAT_REASON = REASON.reshape(-1)


def _project_indices_to_cds(
    idx: np.ndarray, pa_v3: float, fid_pix: tuple[float, float],
    jinv: np.ndarray,
) -> dict:
    """Given a numpy int array of flat shutter indices, return a CDS-shape
    dict {xs, ys, q, s, d}. Vectorized: no per-shutter Python loop."""
    if idx.size == 0:
        return dict(xs=[], ys=[], q=[], s=[], d=[])
    M = idx.size
    v2_centers = _V2_OFFSETS_ALL[idx]
    v3_centers = _V3_OFFSETS_ALL[idx]
    v2_corners = v2_centers[:, None] + _SHUTTER_CORNER_TEMPLATE[None, :, 0]
    v3_corners = v3_centers[:, None] + _SHUTTER_CORNER_TEMPLATE[None, :, 1]
    cx, cy = _project_v2v3_offsets_to_pixel(
        v2_corners.ravel(), v3_corners.ravel(), pa_v3, fid_pix, jinv,
    )
    cx = cx.reshape(M, 4)
    cy = cy.reshape(M, 4)
    qs = (idx // (171 * 365)) + 1
    ss = ((idx % (171 * 365)) // 365) + 1
    ds = (idx % 365) + 1
    # MultiPolygons expects xs/ys as list-of-list-of-list. Build via list
    # comprehensions over the M rows; this is the only Python loop now and
    # it's O(M) with no dict lookups.
    xs_list = [[[cx[k].tolist()]] for k in range(M)]
    ys_list = [[[cy[k].tolist()]] for k in range(M)]
    return dict(
        xs=xs_list, ys=ys_list,
        q=qs.astype(int).tolist(),
        s=ss.astype(int).tolist(),
        d=ds.astype(int).tolist(),
    )


def _in_view_mask(pa_v3: float, fid_pix: tuple[float, float],
                  jinv: np.ndarray,
                  view_bbox_pix: tuple[float, float, float, float]) -> np.ndarray:
    """Bool mask over all 250k shutter centers — True if inside the view bbox
    plus a 50-pixel margin."""
    px, py = _project_v2v3_offsets_to_pixel(
        _V2_OFFSETS_ALL, _V3_OFFSETS_ALL, pa_v3, fid_pix, jinv,
    )
    xmin, xmax, ymin, ymax = view_bbox_pix
    margin = 50.0
    return (
        (px >= xmin - margin) & (px <= xmax + margin)
        & (py >= ymin - margin) & (py <= ymax + margin)
    )


def _fixed_slit_polygons(pa_v3: float, fid_pix: tuple[float, float],
                         jinv: np.ndarray) -> dict:
    """Per-slit polygons in image-pixel coords for the 5 NIRSpec fixed slits."""
    xs_all, ys_all, names = [], [], []
    for name, corners_v2v3 in fixed_slit_corners_v2v3().items():
        v2_off = corners_v2v3[:, 0] - MSA_V2_REF
        v3_off = corners_v2v3[:, 1] - MSA_V3_REF
        x, y = _project_v2v3_offsets_to_pixel(v2_off, v3_off, pa_v3, fid_pix, jinv)
        xs_all.append([[x.tolist()]])
        ys_all.append([[y.tolist()]])
        names.append(name.replace("NRS_", "").replace("_SLIT", ""))
    return dict(xs=xs_all, ys=ys_all, name=names)


def _msa_outline_polygons(pa_v3: float, fid_pix: tuple[float, float],
                          jinv: np.ndarray) -> dict:
    """Trace each quadrant outline (using its 4 corner shutters)."""
    xs_all, ys_all = [], []
    rows = [0, 0, 170, 170]
    cols = [0, 364, 364, 0]
    for q in range(4):
        v2c = np.array([V2_MSA[q, r, c] for r, c in zip(rows, cols)])
        v3c = np.array([V3_MSA[q, r, c] for r, c in zip(rows, cols)])
        x, y = _project_v2v3_offsets_to_pixel(
            v2c - MSA_V2_REF, v3c - MSA_V3_REF, pa_v3, fid_pix, jinv,
        )
        xs_all.append([[x.tolist()]])
        ys_all.append([[y.tolist()]])
    return dict(xs=xs_all, ys=ys_all)


def _shutter_polys_to_cds(polys: dict, only_reason: int | None = None) -> dict:
    """Convert shutter polygons dict to MultiPolygons CDS data.

    If `only_reason` is set, include only shutters whose reason matches.
    Returned dict matches the CDS schema (no `reason` column — that's the filter).
    """
    xs, ys, qs, ss, ds = [], [], [], [], []
    for (q, s, d), p in polys.items():
        if only_reason is not None and p["reason"] != only_reason:
            continue
        xs.append([[p["xs"]]])
        ys.append([[p["ys"]]])
        qs.append(q)
        ss.append(s)
        ds.append(d)
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds)


def _polygons_for_shutter_keys(
    keys, pa_v3: float, fid_pix: tuple[float, float], jinv: np.ndarray,
) -> tuple[list, list, list, list, list]:
    """Vectorized: produce (xs, ys, q, s, d) lists for an iterable of (q,s,d).
    Empty inputs return empty lists.
    """
    keys = list(keys)
    if not keys:
        return [], [], [], [], []
    qs = np.array([k[0] for k in keys], dtype=np.int32)
    ss = np.array([k[1] for k in keys], dtype=np.int32)
    ds = np.array([k[2] for k in keys], dtype=np.int32)
    v2_centers = V2_MSA[qs - 1, ss - 1, ds - 1] - MSA_V2_REF
    v3_centers = V3_MSA[qs - 1, ss - 1, ds - 1] - MSA_V3_REF
    # Broadcast: (M, 4) per axis after adding corner template
    v2_corners = v2_centers[:, None] + _SHUTTER_CORNER_TEMPLATE[None, :, 0]
    v3_corners = v3_centers[:, None] + _SHUTTER_CORNER_TEMPLATE[None, :, 1]
    cx, cy = _project_v2v3_offsets_to_pixel(
        v2_corners.ravel(), v3_corners.ravel(), pa_v3, fid_pix, jinv,
    )
    M = len(keys)
    cx = cx.reshape(M, 4)
    cy = cy.reshape(M, 4)
    xs = [[[cx[k].tolist()]] for k in range(M)]
    ys = [[[cy[k].tolist()]] for k in range(M)]
    return xs, ys, qs.tolist(), ss.tolist(), ds.tolist()


def _highlighted_polygons(pa_v3: float, fid_pix: tuple[float, float],
                          jinv: np.ndarray) -> dict:
    """Polygons for shutters in state['highlighted'], vectorized."""
    if not state["highlighted"]:
        return dict(xs=[], ys=[], q=[], s=[], d=[])
    xs, ys, qs, ss, ds = _polygons_for_shutter_keys(
        state["highlighted"], pa_v3, fid_pix, jinv,
    )
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds)


def _open_shutters_cds_data(pa_v3: float, fid_pix: tuple[float, float],
                            jinv: np.ndarray) -> dict:
    """Open-shutter polygons + per-shutter target ID and wavelength cutoffs."""
    if not state["open_shutters"]:
        return dict(xs=[], ys=[], q=[], s=[], d=[], target=[], lam_blue=[], lam_red=[])
    keys = list(state["open_shutters"].keys())
    xs, ys, qs, ss, ds = _polygons_for_shutter_keys(keys, pa_v3, fid_pix, jinv)
    tgt: list[str] = []
    lam_b: list[float] = []
    lam_r: list[float] = []
    for (q, s, d) in keys:
        sh = state["open_shutters"][(q, s, d)]
        tgt.append(str(sh.target_id) if sh.target_id is not None else "")
        v2c = float(V2_MSA[q - 1, s - 1, d - 1])
        v3c = float(V3_MSA[q - 1, s - 1, d - 1])
        cut = cutoffs(v2c, v3c, state["disperser"], state["filter"])
        lam_b.append(cut.get("lam_blue") if cut and cut.get("lam_blue") is not None else float("nan"))
        lam_r.append(cut.get("lam_red") if cut and cut.get("lam_red") is not None else float("nan"))
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds, target=tgt, lam_blue=lam_b, lam_red=lam_r)


# ---------------------------------------------------------------------------
# Refresh logic
# ---------------------------------------------------------------------------


def refresh_overlays_light() -> None:
    """Cheap, real-time overlay update.

    Renders only the layers that are O(1) in polygon count: the 4 MSA
    quadrant outlines, the 5 fixed slits, and the pointing handle. Skips
    the per-shutter layers (operable, stuck-open, spectral-overlap, open,
    highlighted) and the target catalog. Intended for slider drags so the
    UI stays responsive; pair with refresh_overlays() on commit.
    """
    img: LoadedImage | None = state["image"]
    if img is None:
        return
    fiducial = _pointing_skycoord()
    if fiducial is None:
        return
    pa_v3 = state["pa_v3"]
    wcs = img.wcs
    fid_x, fid_y = _world_to_pixel(fiducial, wcs)
    fid_pix = (float(fid_x), float(fid_y))
    jinv = _compute_wcs_jacobian(wcs, fid_pix[0], fid_pix[1])

    # MSA outline (toggle still respected)
    if 0 in layers_box.active:
        src_msa_outline.data = _msa_outline_polygons(pa_v3, fid_pix, jinv)
    # Fixed slits always visible
    src_fixed_slits.data = _fixed_slit_polygons(pa_v3, fid_pix, jinv)
    src_pointing_handle.data = dict(x=[fid_pix[0]], y=[fid_pix[1]])


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

    # Compute the pointing pixel and the local WCS inverse-Jacobian once;
    # all overlay polygons reuse them via vectorized math.
    fid_x, fid_y = _world_to_pixel(fiducial, wcs)
    fid_pix = (float(fid_x), float(fid_y))
    jinv = _compute_wcs_jacobian(wcs, fid_pix[0], fid_pix[1])

    # MSA outline
    show_outline = 0 in layers_box.active
    if show_outline:
        src_msa_outline.data = _msa_outline_polygons(pa_v3, fid_pix, jinv)
    else:
        src_msa_outline.data = dict(xs=[], ys=[])

    # Cull to the current visible figure bounds (post-zoom), clipped to the
    # image. With the current view-bbox approach, all three shutter-flavour
    # layers (operable, stuck-open, spectral-overlap) reuse one mask.
    try:
        vx0 = max(0.0, min(float(fig.x_range.start), W))
        vx1 = max(0.0, min(float(fig.x_range.end), W))
        vy0 = max(0.0, min(float(fig.y_range.start), H))
        vy1 = max(0.0, min(float(fig.y_range.end), H))
        view_bbox = (min(vx0, vx1), max(vx0, vx1), min(vy0, vy1), max(vy0, vy1))
    except (TypeError, ValueError):
        view_bbox = (0.0, float(W), 0.0, float(H))
    in_view = _in_view_mask(pa_v3, fid_pix, jinv, view_bbox)

    # Stuck-open (always visible, very few): index-based projection.
    stuck_open_idx = np.where(in_view & (_FLAT_REASON == 2))[0]
    src_stuck_open.data = _project_indices_to_cds(stuck_open_idx, pa_v3, fid_pix, jinv)

    # Operable shutters: only if the layer is toggled on. Cap render count.
    show_shutters = 1 in layers_box.active
    if show_shutters:
        op_idx = np.where(in_view & (_FLAT_REASON == 0))[0]
        if op_idx.size > MAX_OPERABLE_RENDER:
            stride = max(1, op_idx.size // MAX_OPERABLE_RENDER)
            op_idx = op_idx[::stride]
        src_bg_shutters.data = _project_indices_to_cds(op_idx, pa_v3, fid_pix, jinv)
    else:
        src_bg_shutters.data = dict(xs=[], ys=[], q=[], s=[], d=[])

    # Spectral-overlap: operable shutters at the same s-row (same quadrant)
    # AND within the V2 distance over which their wavelength coverage overlaps
    # the open shutter's coverage. For PRISM/CLEAR this is essentially the
    # full s-row; for narrow filters (e.g. G140M/F070LP, 0.57 μm) it's only
    # ~20″ — only neighbouring shutters share any wavelength range.
    open_keys = state["open_shutters"].keys()
    if open_keys:
        v2_overlap = float(v2_overlap_distance(state["disperser"], state["filter"]))
        q_arr = np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) // (171 * 365)
        s_arr = (np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) % (171 * 365)) // 365
        is_overlap = np.zeros(_V2_OFFSETS_ALL.size, dtype=bool)
        for q_o, s_o, d_o in open_keys:
            v2_o = float(_V2_OFFSETS_ALL[(q_o - 1) * 171 * 365 + (s_o - 1) * 365 + (d_o - 1)] + MSA_V2_REF)
            same_row = (q_arr == q_o - 1) & (s_arr == s_o - 1)
            near_v2 = np.abs((_V2_OFFSETS_ALL + MSA_V2_REF) - v2_o) < v2_overlap
            is_overlap |= same_row & near_v2
        # Don't mark the open shutters themselves as overlap.
        open_flat = np.array(
            [(k[0]-1) * 171 * 365 + (k[1]-1) * 365 + (k[2]-1) for k in open_keys],
            dtype=np.int64,
        )
        is_overlap[open_flat] = False
        overlap_idx = np.where(in_view & (_FLAT_REASON == 0) & is_overlap)[0]
    else:
        overlap_idx = np.empty(0, dtype=np.int64)
    src_spec_overlap.data = _project_indices_to_cds(overlap_idx, pa_v3, fid_pix, jinv)

    # Open shutters (always shown)
    src_open_shutters.data = _open_shutters_cds_data(pa_v3, fid_pix, jinv)
    # Highlighted shutters
    src_highlighted.data = _highlighted_polygons(pa_v3, fid_pix, jinv)
    # Fixed slits always visible
    src_fixed_slits.data = _fixed_slit_polygons(pa_v3, fid_pix, jinv)
    # Pointing handle at the image-pixel of (RA, Dec)
    src_pointing_handle.data = dict(x=[fid_pix[0]], y=[fid_pix[1]])

    # Targets
    show_targets = 2 in layers_box.active
    cat: Catalog | None = state["catalog"]
    if show_targets and cat is not None:
        coords = SkyCoord(cat.ra_deg, cat.dec_deg, unit=u.deg, frame="icrs")
        x, y = _world_to_pixel(coords, wcs)
        mask = (x >= 0) & (x < W) & (y >= 0) & (y < H)
        # Apply optional catalog filters (priority class ≤, magnitude ≤).
        try:
            pr_cutoff = float(catalog_priority_input.value.strip())
        except (TypeError, ValueError):
            pr_cutoff = None
        try:
            mag_cutoff = float(catalog_mag_input.value.strip())
        except (TypeError, ValueError):
            mag_cutoff = None
        if pr_cutoff is not None:
            # Drop NaNs (treated as "no priority", excluded) and apply ≤
            mask &= np.where(np.isnan(cat.priority), False, cat.priority <= pr_cutoff)
        if mag_cutoff is not None:
            mask &= np.where(np.isnan(cat.mag), False, cat.mag <= mag_cutoff)
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

    # Update stats panel and status with the current configuration.
    n_op = len(state["open_shutters"])
    n_tgt_open = len({sh.target_id for sh in state["open_shutters"].values() if sh.target_id})
    n_hl = len(state["highlighted"])
    n_overlap = overlap_idx.size
    apa = (pa_v3 + V3_IDL_Y_ANGLE) % 360.0
    # Sexagesimal RA/Dec for the stats panel
    fid_sky = SkyCoord(fiducial.ra.deg, fiducial.dec.deg, unit=u.deg)
    ra_hms = fid_sky.ra.to_string(unit=u.hour, sep=":", precision=2, pad=True)
    dec_dms = fid_sky.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)
    stats_div.text = (
        f"<b>Pointing</b>: {fiducial.ra.deg:.5f}, {fiducial.dec.deg:.5f} "
        f"&nbsp;<span style='color:#666'>({ra_hms}, {dec_dms})</span><br>"
        f"<b>V3 PA</b> {pa_v3:.2f}° &nbsp;|&nbsp; <b>NIRSpec APA</b> {apa:.2f}°<br>"
        f"<b>Open</b>: {n_op} shutters across {n_tgt_open} targets &nbsp;|&nbsp; "
        f"<b>Highlight</b>: {n_hl}<br>"
        f"<b>Spectral conflicts</b>: {n_overlap} shutters &nbsp;|&nbsp; "
        f"<b>{state['disperser']}/{state['filter']}</b>"
    )
    _set_status(
        f"{n_op} open shutters covering {n_tgt_open} targets. "
        f"V3 PA = {pa_v3:.2f}°.", "ok",
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
    # Axis tick formatters: convert pixel ticks to RA/Dec degrees using the
    # WCS. Linear approximation around the image center — accurate at the
    # ~milliarcsec level for fields up to ~10 arcmin (so good for our use).
    try:
        H_half, W_half = H / 2.0, W / 2.0
        center = img.wcs.pixel_to_world(W_half, H_half)
        # Differentials in deg/pix at the center
        dx_ra, dx_dec = img.wcs.pixel_to_world(W_half + 1, H_half).ra.deg - center.ra.deg, \
                        img.wcs.pixel_to_world(W_half + 1, H_half).dec.deg - center.dec.deg
        dy_ra, dy_dec = img.wcs.pixel_to_world(W_half, H_half + 1).ra.deg - center.ra.deg, \
                        img.wcs.pixel_to_world(W_half, H_half + 1).dec.deg - center.dec.deg
        # Wrap RA jumps near 0/360 into the linear differential
        for _v in (dx_ra, dy_ra):
            if abs(_v) > 180:
                pass  # left as-is; small-field case won't hit this
        fig.xaxis.formatter = CustomJSTickFormatter(
            args=dict(ra0=center.ra.deg, dec0=center.dec.deg,
                      cx=W_half, cy=H_half,
                      dxra=dx_ra, dxdec=dx_dec, dyra=dy_ra, dydec=dy_dec),
            code="""
                // Linear approximation: tick x is on the central row.
                var ra = ra0 + (tick - cx) * dxra;
                return ra.toFixed(4);
            """,
        )
        fig.yaxis.formatter = CustomJSTickFormatter(
            args=dict(ra0=center.ra.deg, dec0=center.dec.deg,
                      cx=W_half, cy=H_half,
                      dxra=dx_ra, dxdec=dx_dec, dyra=dy_ra, dydec=dy_dec),
            code="""
                var dec = dec0 + (tick - cy) * dydec;
                return dec.toFixed(4);
            """,
        )
    except Exception:  # noqa: BLE001
        pass


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
    finally:
        _hide_loading()


def _load_jpg_pair_from_paths(jpg_path: str, sidecar_path: str) -> None:
    try:
        img = load_jpg_with_sidecar(jpg_path, sidecar_path, max_dim=6000)
        _set_image_and_recenter(img, f"JPG+sidecar {Path(jpg_path).name}")
    except Exception as e:  # noqa: BLE001
        _set_status(f"JPG+sidecar load failed: {e}", "err")
        traceback.print_exc()
    finally:
        _hide_loading()


def _load_catalog_from_path(path: str) -> None:
    try:
        cat = load_catalog(path)
        state["catalog"] = cat
        refresh_overlays()
        _set_status(f"Catalog loaded: {len(cat.ra_deg)} targets from {Path(path).name}.", "ok")
    except Exception as e:  # noqa: BLE001
        _set_status(f"Catalog load failed: {e}", "err")
        traceback.print_exc()
    finally:
        _hide_loading()


# Path-based callbacks (primary input for a local tool).
# Slow loads are deferred to the next tick so the loading banner renders first.
def on_fits_path(attr, old, new):
    p = fits_path_input.value.strip()
    if not p:
        return
    if not Path(p).exists():
        _set_status(f"FITS path not found: {p}", "err")
        return
    _show_loading(f"Loading FITS: {Path(p).name}…")
    _deferred(_load_fits_from_path, p)


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
        _show_loading(f"Loading JPG + WCS sidecar…")
        _deferred(_load_jpg_pair_from_paths, jpg_p, p)
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
    _show_loading(f"Loading JPG: {Path(jpg_p).name}… (large JPGs take 5–10 s)")
    _deferred(_load_jpg_pair_from_paths, jpg_p, side_p)


def on_catalog_path(attr, old, new):
    p = catalog_path_input.value.strip()
    if not p:
        return
    if not Path(p).exists():
        _set_status(f"Catalog path not found: {p}", "err")
        return
    _show_loading(f"Loading catalog: {Path(p).name}…")
    _deferred(_load_catalog_from_path, p)


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


def _sync_pa_widgets(v3pa: float, source: str | None = None) -> None:
    """Update PA widgets from a V3 PA value (mod 360).

    `source` names the widget the change came from so we *don't* write back
    to it (avoids clobbering the user's in-flight typed text).
    """
    state["_syncing_pa"] = True
    try:
        state["pa_v3"] = v3pa % 360.0
        if source != "slider":
            v3pa_slider.value = state["pa_v3"]
        if source != "v3pa_text":
            v3pa_input.value = f"{state['pa_v3']:.2f}"
        if source != "apa_text":
            apa = (state["pa_v3"] + V3_IDL_Y_ANGLE) % 360.0
            apa_input.value = f"{apa:.2f}"
    finally:
        state["_syncing_pa"] = False


def on_v3pa_slider(attr, old, new):
    """Fires continuously during slider drag. Light refresh only — just the
    4 MSA quadrant outlines and the pointing handle follow the value
    in real time. The full shutter overlay is recomputed on slider
    release via on_v3pa_slider_done(value_throttled)."""
    if state.get("_syncing_pa"):
        return
    _sync_pa_widgets(float(v3pa_slider.value), source="slider")
    refresh_overlays_light()


def on_v3pa_slider_done(attr, old, new):
    """Fires once when slider drag ends. Full refresh."""
    if state.get("_syncing_pa"):
        return
    refresh_overlays()


def on_v3pa_text(attr, old, new):
    if state.get("_syncing_pa"):
        return
    try:
        v = float(v3pa_input.value)
    except (TypeError, ValueError):
        return
    _sync_pa_widgets(v, source="v3pa_text")
    refresh_overlays()


def on_apa_text(attr, old, new):
    if state.get("_syncing_pa"):
        return
    try:
        apa = float(apa_input.value)
    except (TypeError, ValueError):
        return
    _sync_pa_widgets(apa - V3_IDL_Y_ANGLE, source="apa_text")
    refresh_overlays()


def on_layers(attr, old, new):
    refresh_overlays()


def on_disperser_filter(attr, old, new):
    label = disperser_filter_select.value
    try:
        d_part, f_part = [s.strip() for s in label.split("/")]
    except ValueError:
        return
    state["disperser"] = d_part
    state["filter"] = f_part
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


def _data_distance_for_screen_pixels(n_pix: float = 15.0) -> float:
    """Convert a screen-pixel distance to a rough data-coord distance at the
    current zoom level, for proximity tests. Uses the x-range / fig.width."""
    try:
        rng = float(fig.x_range.end) - float(fig.x_range.start)
        w = max(int(getattr(fig, "width", 900) or 900), 1)
        return n_pix * (rng / w)
    except Exception:  # noqa: BLE001
        return 30.0


def _open_or_toggle_slitlet_at(q: int, s: int, d: int, target_id: str | None) -> None:
    """Open a slitlet centered at (q, s, d). If it's already open (same key),
    remove it plus its slitlet siblings instead.
    """
    key = (q, s, d)
    _push_history()
    if key in state["open_shutters"]:
        sh = state["open_shutters"].pop(key)
        if sh.target_id:
            for k in list(state["open_shutters"].keys()):
                other = state["open_shutters"][k]
                if (other.target_id == sh.target_id
                        and other.q == q and other.d == d
                        and abs(other.s - s) <= state["slitlet_height"] // 2):
                    del state["open_shutters"][k]
        _set_status(f"Closed shutter ({q},{s},{d}) and slitlet siblings.", "ok")
    elif target_id is not None:
        n = _add_slitlet(q, s, d, target_id=target_id)
        _set_status(f"Target {target_id} → slitlet ({q},{s},{d}), {n} shutters opened.", "ok")
    else:
        state["open_shutters"][key] = OpenShutter(q=q, s=s, d=d, role="manual")
        _set_status(f"Opened shutter ({q},{s},{d}) manually.", "ok")
    refresh_overlays()


def on_tap(event):
    """Unified single-tap handler.

    Snap-to-nearest model so the user doesn't have to land precisely on a
    polygon. Modes:
      - Shift+click: move the pointing center to the click location.
      - Click near a yellow target: open a 3-shutter slitlet there.
      - Click anywhere else: snap to the nearest shutter and toggle it.
    """
    img = state["image"]
    fiducial = _pointing_skycoord()
    if img is None or fiducial is None:
        return
    x_data, y_data = float(event.x), float(event.y)

    # 0) Shift-click → move pointing center. Bokeh 3.x: event.modifiers is a
    # KeyModifiers model with .shift attribute; older versions exposed it as
    # a dict or as a top-level event attribute.
    mods = getattr(event, "modifiers", None)
    shift_held = False
    if mods is not None:
        if isinstance(mods, dict):
            shift_held = bool(mods.get("shift"))
        else:
            shift_held = bool(getattr(mods, "shift", False))
    if not shift_held:
        shift_held = bool(getattr(event, "shift", False))
    if shift_held:
        _move_pointing_to(x_data, y_data)
        return

    # 1) Near a target marker? Open a slitlet for it.
    tx = np.asarray(src_targets.data.get("x") or [], dtype=float)
    ty = np.asarray(src_targets.data.get("y") or [], dtype=float)
    if tx.size > 0:
        d2 = (tx - x_data) ** 2 + (ty - y_data) ** 2
        i = int(np.argmin(d2))
        thr = _data_distance_for_screen_pixels(15.0)
        if d2[i] < thr * thr:
            tgt_id = str(src_targets.data["id"][i])
            try:
                sky_target = img.wcs.pixel_to_world(float(tx[i]), float(ty[i]))
                v2, v3 = _sky_to_v2v3(sky_target, fiducial, state["pa_v3"])
                nearest = _nearest_shutter(v2, v3, require_operable=state["snap_to_operable"])
            except Exception:  # noqa: BLE001
                nearest = None
            if nearest is None:
                _set_status(f"Target {tgt_id}: no operable shutter nearby.", "warn")
                return
            q, s, d = nearest
            _open_or_toggle_slitlet_at(q, s, d, target_id=tgt_id)
            return

    # 2) Otherwise: snap to nearest operable shutter and toggle it manually.
    # Non-operable shutters (stuck-open / failed-closed) must never be opened.
    try:
        sky = img.wcs.pixel_to_world(x_data, y_data)
        v2, v3 = _sky_to_v2v3(sky, fiducial, state["pa_v3"])
        nearest = _nearest_shutter(v2, v3, require_operable=True)
    except Exception:  # noqa: BLE001
        return
    if nearest is None:
        _set_status("No operable shutter near that click.", "warn")
        return
    q, s, d = nearest
    _open_or_toggle_slitlet_at(q, s, d, target_id=None)


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
        # Also drop a session.json next to the APT files so the same
        # directory can be re-loaded later via Session → Load to restore
        # picks + pointing + everything.
        export_session_json(_build_current_session(), str(out_dir / "session.json"))
        _set_status(
            f"Exported eMPT bundle + session.json to {out_dir} "
            f"({len(targets_rows)} targets, {len(open_list)} open shutters).",
            "ok", clear_after=15,
        )
        # Pre-fill the Session-load input so the user can re-load with one click.
        session_load_path_input.value = str(out_dir / "session.json")
    except Exception as e:  # noqa: BLE001
        _set_status(f"Export failed: {e}", "err")
        traceback.print_exc()


# Tap wiring
fig.on_event(Tap, on_tap)


def on_ranges_update(event):
    """Fires when the user finishes pan or zoom. Re-cull shutters to the new
    view so panning into a fresh region brings its shutters into view."""
    refresh_overlays()


fig.on_event(RangesUpdate, on_ranges_update)
undo_btn.on_click(on_undo)
clear_btn.on_click(on_clear)
export_btn.on_click(on_export)


# ---------------------------------------------------------------------------
# Session save/load
# ---------------------------------------------------------------------------


def _build_current_session() -> Session:
    return Session(
        pointing_ra_deg=float(state["ra_deg"]),
        pointing_dec_deg=float(state["dec_deg"]),
        pa_v3_deg=float(state["pa_v3"]),
        disperser=state["disperser"],
        filter_name=state["filter"],
        slitlet_height=int(state["slitlet_height"]),
        open_shutters=list(state["open_shutters"].values()),
        highlighted=list(state["highlighted"]),
        image_path=(state["image"].source_path if state["image"] else None),
        catalog_path=(state["catalog"].source_path if state["catalog"] else None),
    )


def on_session_save():
    path = session_save_path_input.value.strip()
    if not path:
        _set_status("Set a session save path first.", "warn")
        return
    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        export_session_json(_build_current_session(), str(p))
        _set_status(f"Session saved → {p}", "ok", clear_after=10)
    except Exception as e:  # noqa: BLE001
        _set_status(f"Session save failed: {e}", "err")
        traceback.print_exc()


def on_session_load():
    path = session_load_path_input.value.strip()
    if not path:
        _set_status("Enter a session load path first.", "warn")
        return
    p = Path(path).expanduser()
    if not p.exists():
        _set_status(f"Session file not found: {p}", "err")
        return
    try:
        sess = import_session_json(str(p))
    except ValueError as e:
        _set_status(f"Session load failed: {e}", "err")
        return
    # Apply session to the live UI. Load image/catalog first if paths still exist.
    if sess.image_path and Path(sess.image_path).exists() and (
        state["image"] is None
        or getattr(state["image"], "source_path", None) != sess.image_path
    ):
        fits_path_input.value = sess.image_path  # triggers on_fits_path
    if sess.catalog_path and Path(sess.catalog_path).exists() and (
        state["catalog"] is None
        or getattr(state["catalog"], "source_path", None) != sess.catalog_path
    ):
        catalog_path_input.value = sess.catalog_path
    # Pointing & PA — set via the inputs so the standard on_change handlers run
    ra_input.value = f"{sess.pointing_ra_deg:.6f}"
    dec_input.value = f"{sess.pointing_dec_deg:.6f}"
    _sync_pa_widgets(sess.pa_v3_deg)
    combo_label = f"{sess.disperser} / {sess.filter_name}"
    if combo_label in DISPERSER_FILTER_LABELS:
        disperser_filter_select.value = combo_label
    else:
        # Fall back to direct state mutation if the combo isn't in the list
        state["disperser"] = sess.disperser
        state["filter"] = sess.filter_name
    slitlet_select.value = str(sess.slitlet_height)
    state["slitlet_height"] = int(sess.slitlet_height)
    # Restore the open-shutter dict and highlighted set
    state["open_shutters"] = {
        (sh.q, sh.s, sh.d): sh for sh in sess.open_shutters
    }
    state["highlighted"] = set(sess.highlighted)
    state["history"] = []
    refresh_overlays()
    _set_status(
        f"Session loaded: {len(state['open_shutters'])} open shutters, "
        f"{len(state['highlighted'])} highlighted.", "ok", clear_after=10,
    )


session_save_btn.on_click(on_session_save)
session_load_btn.on_click(on_session_load)


# ---------------------------------------------------------------------------
# Example data quick-load
# ---------------------------------------------------------------------------


_EX_DIR = Path(__file__).resolve().parent.parent

def on_example_a370():
    p = _EX_DIR / "example_a370" / "a370_f182m_f200w_f210m.fits"
    if not p.exists():
        _set_status(f"Example missing: {p}", "err")
        return
    fits_path_input.value = str(p)


def on_example_r0600():
    jpg = _EX_DIR / "example_r0600" / "JWST_F090W_F200W_F444W.jpg"
    wcs = _EX_DIR / "example_r0600" / "wcs.fits"
    if not (jpg.exists() and wcs.exists()):
        _set_status("Example r0600 files missing.", "err")
        return
    # Set sidecar first so the JPG callback finds it.
    sidecar_path_input.value = str(wcs)
    jpg_path_input.value = str(jpg)


example_a370_btn.on_click(on_example_a370)
example_r0600_btn.on_click(on_example_r0600)
snap_box.on_change("active", on_snap)


# ---------------------------------------------------------------------------
# Drag-pointing handle and double-tap highlight
# ---------------------------------------------------------------------------


def _move_pointing_to(x_data: float, y_data: float) -> None:
    """Move the pointing center to the given image-pixel position."""
    img = state["image"]
    if img is None:
        return
    try:
        sky = img.wcs.pixel_to_world(float(x_data), float(y_data))
        ra_input.value = f"{sky.ra.deg:.6f}"
        dec_input.value = f"{sky.dec.deg:.6f}"
        _set_status(
            f"Pointing → RA={sky.ra.deg:.5f}°, Dec={sky.dec.deg:.5f}°.", "ok",
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"Move-pointing failed: {e}", "err")




def on_double_tap(event):
    """Double-click anywhere → toggle highlight on the nearest shutter."""
    img = state["image"]
    fiducial = _pointing_skycoord()
    if img is None or fiducial is None:
        return
    try:
        sky = img.wcs.pixel_to_world(float(event.x), float(event.y))
    except Exception:  # noqa: BLE001
        return
    v2, v3 = _sky_to_v2v3(sky, fiducial, state["pa_v3"])
    nearest = _nearest_shutter(v2, v3, require_operable=False)
    if nearest is None:
        return
    key = nearest
    if key in state["highlighted"]:
        state["highlighted"].remove(key)
        _set_status(f"Un-highlighted shutter {key}.", "ok")
    else:
        state["highlighted"].add(key)
        _set_status(f"Highlighted shutter {key} (cyan).", "ok")
    refresh_overlays()


fig.on_event(DoubleTap, on_double_tap)


# ---------------------------------------------------------------------------
# Visibility / APA_V3 constraints (jwst_gtvt)
# ---------------------------------------------------------------------------


def on_visibility():
    fiducial = _pointing_skycoord()
    if fiducial is None:
        _set_status("Set RA/Dec before computing visibility.", "warn")
        return
    try:
        from jwst_gtvt.jwst_tvt import Ephemeris  # noqa: F401
    except ImportError:
        _set_status("jwst_gtvt not installed (pip install jwst_gtvt).", "err")
        return
    _show_loading("Querying jwst_gtvt ephemeris… (first call ~5–8 s, cached after).")
    _deferred(_do_visibility_query, fiducial)


def _do_visibility_query(fiducial):
    try:
        from jwst_gtvt.jwst_tvt import Ephemeris
        eph = Ephemeris()
        df = eph.get_fixed_target_positions(
            f"{fiducial.ra.deg:.6f}", f"{fiducial.dec.deg:.6f}"
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"jwst_gtvt query failed: {e}", "err")
        traceback.print_exc()
        _hide_loading()
        return
    _hide_loading()
    _render_visibility(df)


def _render_visibility(df):

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
    _sync_pa_widgets(nom)
    refresh_overlays()
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
catalog_priority_input.on_change("value", lambda a, o, n: refresh_overlays())
catalog_mag_input.on_change("value", lambda a, o, n: refresh_overlays())
ra_input.on_change("value", on_pointing)
dec_input.on_change("value", on_pointing)
v3pa_slider.on_change("value", on_v3pa_slider)
v3pa_slider.on_change("value_throttled", on_v3pa_slider_done)
v3pa_input.on_change("value", on_v3pa_text)
apa_input.on_change("value", on_apa_text)
layers_box.on_change("active", on_layers)
disperser_filter_select.on_change("value", on_disperser_filter)
slitlet_select.on_change("value", on_slitlet_height)

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

sidebar = column(
    loading_banner,
    stats_div,
    Div(text="<h3>Image</h3>Try an example field:"),
    row(example_a370_btn, example_r0600_btn),
    Div(text="<small><b>or</b> paste a local path:</small>"),
    fits_path_input,
    Div(text="<i>or</i> JPG + sidecar FITS by path:"),
    sidecar_path_input,
    jpg_path_input,
    Div(text="<small>Upload widgets below — small files only "
             "(unless server started with <code>--websocket-max-message-size</code>):</small>"),
    fits_input,
    sidecar_input,
    jpg_input,
    Div(text="<h3>Catalog</h3>Path:"),
    catalog_path_input,
    catalog_input,
    catalog_priority_input,
    catalog_mag_input,
    Div(text="<h3>Pointing</h3>"),
    ra_input, dec_input,
    v3pa_slider, v3pa_input, apa_input, pa_help_div,
    visibility_date_input, visibility_btn, visibility_div,
    Div(text="<h3>Instrument</h3>"),
    disperser_filter_select,
    Div(text="<h3>Display</h3>"),
    layers_box,
    slitlet_select,
    snap_box,
    row(undo_btn, clear_btn),
    Div(text="<h3>Session (collaborate)</h3>"),
    session_save_path_input, session_save_btn,
    session_load_path_input, session_load_btn,
    Div(text="<h3>Export to APT</h3>"),
    export_dir_input, export_btn,
    status,
    width=340,
)

curdoc().add_root(row(sidebar, fig, help_panel))
curdoc().title = "vMPT — visual MSA Planning Tool"
