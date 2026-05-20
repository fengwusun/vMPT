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
    HoverTool,
    Range1d,
    Select,
    Slider,
    TabPanel,
    Tabs,
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
    write_mpt_catalog,
    write_observed_targets_cat,
    write_pointing_summary_txt,
    write_shutter_mask_csv,
)
from app.image_io import LoadedImage, load_fits, load_jpg_with_sidecar, stretch_for_display
from app.msa import load_msa_grid, load_operability
from app.mpt_io import (
    MPTPlan,
    download_apt_program,
    list_mpt_plans_in_aptx,
    parse_mpt_json,
    parse_mpt_json_in_aptx,
    parse_shutter_csv,
)
from app.session_io import (
    EMPT_OBSERVED_FILENAME,
    EMPT_POINTING_FILENAME,
    EMPT_SHUTTER_MASK_FILENAME,
    MPT_CATALOG_FILENAME,
    MPT_PLAN_FILENAME,
    Session,
    WORKSPACE_FILENAME,
    export_session_json,
    import_session_json,
)
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
    # Cache: (q,s,d) → [catalog source id, …] for sources falling inside the
    # shutter footprint at the current pointing + PA. Rebuilt whenever
    # pointing / PA / catalog changes.
    "shutter_to_catids": {},
}


def _push_history() -> None:
    """Snapshot open_shutters for undo (cap at 50)."""
    state["history"].append(dict(state["open_shutters"]))
    if len(state["history"]) > 50:
        state["history"].pop(0)

# ---------------------------------------------------------------------------
# Layout constants — declared early so they're available everywhere
# (help_panel and sidebar reference these during widget construction
# before the figure block where the figure dims are used).
# ---------------------------------------------------------------------------
FIG_W_HINT = 900     # initial canvas width hint (Bokeh stretches it)
FIG_H_HINT = 800     # initial canvas height hint
SIDEBAR_W = 340      # left tab panel (Image / Aim / Pick / MPT)
HELPPANEL_W = 340    # right help panel (Quick guide + rotating tip)


# ---------------------------------------------------------------------------
# Bokeh widgets / glyphs
# ---------------------------------------------------------------------------

status = Div(
    text="Enter a file path below (or use the upload widgets for small files).",
    # Stretch to the sidebar's width (340 px) so the status box can't
    # overflow horizontally — previously the fixed width=420 was wider
    # than the sidebar and the wrapped status text would render on top
    # of the MPT-tab content below it.
    sizing_mode="stretch_width",
    styles={
        "padding": "4px 8px",
        "font-size": "11.5px",
        "line-height": "1.35",
        "border-top": "1px solid #e0e6f0",
        "box-sizing": "border-box",
        "overflow-wrap": "anywhere",
        "word-break": "break-word",
    },
)

# Path-based inputs are the primary way to load — no WebSocket size limit, no temp files.
fits_path_input = TextInput(title="FITS path (local)", value="", placeholder="/path/to/image.fits")
jpg_path_input = TextInput(title="JPG path (local)", value="", placeholder="/path/to/image.jpg")
sidecar_path_input = TextInput(title="Sidecar FITS path (WCS for JPG)", value="", placeholder="/path/to/wcs.fits")
catalog_path_input = TextInput(title="Catalog path (local)", value="", placeholder="/path/to/catalog.csv")

# Upload widgets work for small files but Bokeh's default WebSocket limit (~20 MB) will
# silently truncate larger ones. Start the server with --websocket-max-message-size if
# you want to use these for big files.
# "Browse…" buttons paired with each path TextInput. Clicking one opens
# a native file picker (via a tkinter subprocess) and writes the chosen
# path into the text input — which triggers the existing on_<path>_path
# callback. No upload, no WebSocket size limit.
fits_browse_btn = Button(label="Browse…", button_type="default", width=80)
jpg_browse_btn = Button(label="Browse…", button_type="default", width=80)
sidecar_browse_btn = Button(label="Browse…", button_type="default", width=80)
catalog_browse_btn = Button(label="Browse…", button_type="default", width=80)
mpt_json_browse_btn = Button(label="Browse…", button_type="default", width=80)
mpt_csv_browse_btn = Button(label="Browse…", button_type="default", width=80)
apt_path_browse_btn = Button(label="Browse…", button_type="default", width=80)
session_save_browse_btn = Button(label="Browse…", button_type="default", width=80)
session_load_browse_btn = Button(label="Browse…", button_type="default", width=80)
export_dir_browse_btn = Button(label="Browse…", button_type="default", width=80)

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

# Full-page loading overlay — a centered translucent backdrop with an
# animated spinner. The widget itself is a zero-size Bokeh Div, but its
# inner HTML uses `position: fixed` to escape the layout and cover the
# whole viewport.
loading_banner = Div(
    text="", width=0, height=0, visible=False,
    # The outer Bokeh container must be invisible — only the position-fixed
    # inner overlay should render.
    styles=dict(background="transparent", border="none", padding="0",
                margin="0"),
)

# Quick-help panel rendered on the right side of the figure.
help_toggle_btn = Button(label="Hide help", button_type="default", width=110)

# ── Rotating tip strip ───────────────────────────────────────────────────────
# The help panel shows a rotating one-liner tip on top, and the full
# reference below it. Tips fade in/out every 15 s so the help feels alive
# without being annoying. See `_TIPS` for content; rotation is wired by
# `_advance_tip` registered as a periodic callback further down.
_TIPS = [
    ("🎯", "Pick mode", "Click anywhere on the image — vMPT snaps to the nearest operable shutter and opens an <b>N-shutter slitlet</b> (set N=1/2/3/5 in the <b>Pick</b> tab)."),
    ("✋", "Move the pointing", "<b>Shift + click</b> anywhere on the image to recentre the pointing on that spot. The <span style='color:#2e9b3f;font-weight:600'>lime cross</span> marks the current pointing."),
    ("🔁", "Toggle a slitlet", "Click an already-open shutter to close it. Its slitlet siblings come down with it."),
    ("🎨", "Cyan flag", "Double-click a shutter to toggle a <span style='color:#0aa;font-weight:600'>cyan highlight</span> — a visual flag for your own review. It's not exported."),
    ("🔭", "Pick a roll", "In the <b>Aim</b> tab, enter a visibility date and click <b>Compute allowed V3 PA</b>. jwst_gtvt reports the valid window for the date."),
    ("🌈", "Wavelength check", "Hover any open shutter to see its λ<sub>blue</sub> / λ<sub>red</sub> and the NRS1 / NRS2 detector-gap range for the current disperser."),
    ("⚠️", "Orange = collision", "Orange-tinted shutters share a dispersed-y row with an open or stuck-open shutter — opening them would put two spectra on the same detector pixels."),
    ("🪞", "Cross-quadrant", "Spec-overlap correctly pairs Q1↔Q3 (NRS1) and Q2↔Q4 (NRS2). A pick in Q1 will never light up Q2 or Q4."),
    ("💎", "Catalog match", "Open a shutter with a catalog source inside it — vMPT auto-tags the slitlet with that source's ID. Status bar names the match."),
    ("📤", "Export bundle", "<b>MPT</b> tab → <b>Export eMPT bundle</b> writes a folder with <code>MPT_plan.json</code>, an APT-importable <code>.cat</code> target list, and the eMPT pipeline's three files."),
    ("⏪", "Undo", "<b>Pick</b> tab → <b>Undo last</b> reverts the most recent slitlet open/close action. History is 50 deep."),
    ("📐", "Slitlet sizes", "N=2 means clicked-shutter + one row of lower-y on the detector. N=3/5 are centred on the click. Switch any time in <b>Pick</b>."),
    ("🛰️", "Two ways to load APT", "<b>MPT</b> tab → either point at a local <code>.aptx</code>, or just type a JWST program ID (e.g. <code>1208</code>) and vMPT pulls it from STScI."),
]

tip_div = Div(
    sizing_mode="stretch_width",
    styles=dict(
        background="linear-gradient(135deg, #fef9e7 0%, #fff5d6 100%)",
        color="#5a4a00",
        padding="10px 14px",
        border="1px solid #f0d990",
        **{"border-radius": "6px",
            "transition": "opacity 350ms ease",
            "font-size": "12.5px",
            "line-height": "1.45",
            "min-height": "62px",
        },
    ),
    text="",  # populated below
)


def _render_tip(idx: int) -> str:
    """One-tip card. Each render carries an explicit @keyframes block so a
    fresh DOM swap (Bokeh redoes the inner HTML when `text` changes) restarts
    the fade-in animation."""
    emoji, header, body = _TIPS[idx % len(_TIPS)]
    return (
        '<style>'
        '@keyframes vmpt-tip-fadein { '
        '  0% { opacity: 0; transform: translateY(4px); } '
        '  100% { opacity: 1; transform: translateY(0); } '
        '}'
        '</style>'
        '<div style="display: flex; gap: 10px; align-items: flex-start; '
        '            animation: vmpt-tip-fadein 350ms ease-out;">'
        f'  <div style="font-size: 22px; line-height: 1; flex-shrink: 0;">{emoji}</div>'
        '  <div>'
        '    <div style="font-weight: 700; color: #8a6300; letter-spacing: 0.3px; '
        '                margin-bottom: 3px; font-size: 11.5px; text-transform: uppercase;">'
        f'      Tip · {header}'
        '    </div>'
        f'    <div>{body}</div>'
        '  </div>'
        '</div>'
    )


# Index lives in module-state so the periodic callback can advance it.
_tip_state = {"idx": 0}
tip_div.text = _render_tip(_tip_state["idx"])

help_div = Div(
    # Note: the OUTER help_panel column constrains width — match it here
    # with box-sizing so padding stays inside the 340-px envelope and
    # the long inline tokens like <code>vMPT_workspace.json</code> can
    # wrap rather than overflow horizontally.
    sizing_mode="stretch_width",
    styles={
        "background": "#f8f9fa",
        "color": "#212529",
        "padding": "8px 10px",
        "border": "1px solid #dee2e6",
        "border-radius": "6px",
        "box-sizing": "border-box",
        "font-size": "12px",
        "line-height": "1.4",
        "overflow-wrap": "anywhere",
        "word-break": "break-word",
    },
    text="""
<style>
  /* Local typography for the Quick guide — keeps the nested lists from
     pushing content off the right edge of the 340-px help panel. */
  .vmpt-help h3 { margin: 0 0 6px 0; font-size: 14px; }
  .vmpt-help b  { color: #1a3b66; }
  .vmpt-help ul { margin: 2px 0 6px 0; padding-left: 16px; }
  .vmpt-help ul ul { margin: 1px 0 1px 0; padding-left: 12px; }
  .vmpt-help li { margin: 1px 0; }
  .vmpt-help code { font-size: 11px; padding: 0 2px;
                    background: #ececec; border-radius: 2px;
                    word-break: break-all; }
</style>
<div class="vmpt-help">
<h3>Quick guide</h3>
<b>1. Load an image</b>
<ul>
  <li>One-click <b>Load Abell 370 example</b> or <b>Load RXCJ0600 example</b> from the <b>Image</b> tab — fastest.</li>
  <li>Or paste a local <b>FITS</b> path (with WCS), or a <b>JPG + sidecar FITS</b> pair.</li>
</ul>
<b>2. Optional: target catalog</b>
<ul>
  <li>CSV / ASCII / FITS with at least <code>ID, RA, DEC</code>.</li>
  <li>Targets render as yellow circles. A shutter containing a catalog source auto-tags the slitlet on click.</li>
</ul>
<b>3. Aim the MSA</b>
<ul>
  <li><b>V3 PA</b> drives the math; <b>NIRSpec APA</b> = V3 PA + 138.575° (mod 360).</li>
  <li><b>Shift + click</b> to move pointing. The <span style='color:#2e9b3f;font-weight:600'>lime cross</span> marks it.</li>
  <li>Type a date in <b>Visibility</b> → <b>Compute allowed V3 PA</b> to query jwst_gtvt.</li>
</ul>
<b>4. Hand-pick shutters</b>
<ul>
  <li>Pick the <b>N-shutter slitlet</b> size (1/2/3/5) in <b>Pick</b>.</li>
  <li><b>Click</b> → opens N-shutter slitlet at the nearest operable shutter. Click an open shutter to close the slitlet.</li>
  <li><b>Double-click</b> → toggles <span style='color:#0aa;font-weight:600;background:#222;padding:0 4px'>cyan highlight</span> (visual flag, not exported).</li>
  <li>Layers (Pick tab → <b>Layers</b>):
    <ul>
      <li><span style='background:silver;padding:0 4px'>silver</span> = operable</li>
      <li><span style='color:#d63d3d;font-weight:700'>red fill</span> = your picks</li>
      <li><span style='color:#b30000;font-weight:700'>dark red</span> = stuck-open</li>
      <li><span style='color:#e26a00;font-weight:700'>orange</span> = spec-overlap warning</li>
      <li><span style='color:gold;font-weight:700'>gold</span> = fixed slits</li>
      <li><span style='color:#ddd200;font-weight:700'>yellow ○</span> = catalog target · <span style='color:#2e9b3f;font-weight:700'>green ○</span> = matched to an open shutter</li>
    </ul>
  </li>
</ul>
<b>5. Save / share / export</b>
<ul>
  <li><b>MPT</b> tab → <b>Save session</b> writes a bundle.</li>
  <li><b>Load session</b> — point at <code>MPT_plan.json</code> or <code>vMPT_workspace.json</code>; the sibling auto-loads.</li>
  <li><b>Export eMPT bundle</b> writes a timestamped folder:
    <ul>
      <li><code>MPT_plan.json</code> + <code>&lt;catalog&gt;.cat</code> → APT MPT</li>
      <li><code>vMPT_workspace.json</code> → vMPT round-trip</li>
      <li><code>eMPT_*</code> three files → eMPT pipeline</li>
    </ul>
  </li>
</ul>
<b>Interactions</b>
<ul>
  <li><b>Wheel</b>: zoom · <b>Drag</b>: pan · <b>Box zoom</b>: toolbar → drag</li>
  <li><b>Reset</b>: toolbar · <b>Undo</b>: Pick → <b>Undo last</b></li>
</ul>
</div>
<p style='margin:4px 0'>Full reference in <code>README.md</code> · file roles in <code>CONTEXT.md</code>.</p>
""",
)


def on_help_toggle():
    help_div.visible = not help_div.visible
    tip_div.visible = not tip_div.visible
    help_toggle_btn.label = "Hide help" if help_div.visible else "Show help"


help_toggle_btn.on_click(on_help_toggle)
help_panel = column(
    help_toggle_btn, tip_div, help_div,
    width=HELPPANEL_W,
    # The Quick guide is long. Make the help panel scroll vertically
    # within whatever height it gets in the page layout, so users on
    # smaller screens can still reach the bottom of the guide.
    height_policy="max",
    styles={
        "overflow-y": "auto",
        "overflow-x": "hidden",   # never horizontally scroll the guide
        "max-height": "100vh",
        # Long inline tokens like vMPT_workspace.json shouldn't push the
        # panel wider — break them anywhere if needed.
        "overflow-wrap": "anywhere",
        "word-break": "break-word",
        "box-sizing": "border-box",
    },
)

disperser_filter_select = Select(
    title="Disperser / Filter",
    options=DISPERSER_FILTER_LABELS,
    value="PRISM / CLEAR",
)

layers_box = CheckboxGroup(
    labels=["Show MSA outline", "Show operable shutters", "Show catalog targets"],
    # All three layers on by default. The operable-shutter layer is now
    # filtered to *unaffected* shutters only (excludes user-opens,
    # stuck-opens, and spec-overlap rows), keeping the polygon count
    # manageable so it works at typical zoom levels.
    active=[0, 1, 2],
)
MAX_OPERABLE_RENDER = 10000  # cap for operable-shutter silver-edge layer.
                              # Below the cap we draw every shutter (no
                              # stride); above it we skip rendering — user
                              # zooms in further to see all silver edges.
# Slitlet selector. Each click opens N shutters at the picked column:
#   N=1: just the click
#   N=2: click + the shutter one row below (s-1)
#   N=3: click ±1 (centred)
#   N=5: click ±2 (centred)
slitlet_select = Select(
    title="N-shutter slitlet", options=["1", "2", "3", "5"], value="3",
)

snap_box = CheckboxGroup(labels=["Snap target to nearest operable"], active=[0])

# ── Overlay-appearance picker ────────────────────────────────────────────
# One dropdown + two sliders (alpha, stroke) that retarget themselves at
# whichever layer the user selects. The mapping from each layer name to
# its glyph properties is in `_OVERLAY_LAYER_CONFIG` further down; the
# slider on_change callbacks read that config to decide which glyph
# attribute to mutate.
overlay_layer_select = Select(
    title="Adjust layer",
    options=[
        "Operable shutters",
        "Overlapping shutters",
        "Picked shutters",
        "Stuck open",
        "Catalog sources",
    ],
    value="Operable shutters",
)
overlay_alpha_slider = Slider(
    start=0.0, end=1.0, step=0.05, value=0.20,
    title="Alpha",
    width=SIDEBAR_W - 40,
)
overlay_stroke_slider = Slider(
    start=0.0, end=3.0, step=0.05, value=0.75,
    title="Stroke (px)",
    width=SIDEBAR_W - 40,
)

undo_btn = Button(label="Undo last", button_type="default")
clear_btn = Button(label="Clear open", button_type="warning")

export_dir_input = TextInput(title="Export dir", value=str(Path.cwd() / "exports"))
export_btn = Button(label="Export eMPT bundle", button_type="success")

# Session save/load: round-trips the full picking state for collaborators.
session_save_path_input = TextInput(
    title="Session save path",
    value=str(Path.cwd() / "exports" / MPT_PLAN_FILENAME),
    placeholder=f"/path/to/{MPT_PLAN_FILENAME}",
)
session_save_btn = Button(label="Save session", button_type="primary")
session_load_path_input = TextInput(
    title="Session load path",
    value="",
    placeholder=f"/path/to/{MPT_PLAN_FILENAME} (or {WORKSPACE_FILENAME})",
)
session_load_btn = Button(label="Load session", button_type="primary")

# Example data quick-load buttons (onboarding).
example_a370_btn = Button(label="Load Abell 370 example", button_type="default")
example_r0600_btn = Button(label="Load RXCJ0600 example", button_type="default")

# Import an existing APT plan: either an MPT JSON (preferred — has pointing,
# PA, multiple plans, target IDs) or a shutter CSV (just the open mask).
mpt_json_path_input = TextInput(
    title="APT/MPT plan JSON path",
    value="",
    placeholder="/path/to/plan.json",
)
mpt_plan_select = Select(
    title="Plan from JSON",
    options=[],
    value="",
    visible=False,
)
mpt_load_btn = Button(
    label="Load plan from JSON", button_type="primary", disabled=True,
)
mpt_csv_path_input = TextInput(
    title="Shutter CSV path (open-mask only)",
    value="",
    placeholder="/path/to/shutter_mask.csv",
)
mpt_csv_load_btn = Button(label="Load shutter CSV", button_type="primary")

# Direct APT support: either an .aptx file path on disk, or a program
# ID we can download from STScI.
apt_path_input = TextInput(
    title="APT (.aptx) path on disk",
    value="",
    placeholder="/path/to/1208.aptx",
)
apt_program_input = TextInput(
    title="… or JWST program ID (fetched from STScI)",
    value="",
    placeholder="e.g. 1208",
)
apt_fetch_btn = Button(label="Fetch / open .aptx", button_type="primary")
apt_plan_select = Select(
    title="Plan inside .aptx",
    options=[],
    value="",
    visible=False,
)
apt_load_btn = Button(label="Load selected plan", button_type="primary", disabled=True)

# Status bar — wide strip above the figure showing pointing / PA /
# disperser / open-shutter count / spec-conflict count. Always visible
# regardless of which sidebar tab the user is on; designed so a planner
# can glance up from picking and confirm their state in one read.
stats_div = Div(
    sizing_mode="stretch_width", height=44,
    styles=dict(
        background="linear-gradient(180deg, #f7fbff 0%, #eaf2ff 100%)",
        color="#1a3b66",
        padding="6px 14px",
        border="1px solid #c2d6f0",
        **{"border-radius": "6px",
            "font-size": "13px",
            "line-height": "1.35",
            "box-shadow": "0 1px 2px rgba(0,0,0,0.04)",
        },
    ),
    text="<i>Loading vMPT… pick an example from the Image tab to begin.</i>",
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
    data=dict(xs=[], ys=[], q=[], s=[], d=[], target=[],
              lam_blue=[], lam_red=[], lam_gap_lo=[], lam_gap_hi=[],
              gap_label=[])
)
src_highlighted = ColumnDataSource(data=dict(xs=[], ys=[], q=[], s=[], d=[]))
src_spec_overlap = ColumnDataSource(data=dict(xs=[], ys=[], q=[], s=[], d=[]))
src_fixed_slits = ColumnDataSource(data=dict(xs=[], ys=[], name=[]))
src_targets = ColumnDataSource(data=dict(
    x=[], y=[], id=[], ra=[], dec=[], pr=[],
    # Per-target rendering: matched sources (their id is the target_id of
    # at least one open shutter) flip to green with a thicker line so the
    # user can see at a glance which catalogue entries they've already
    # placed in slitlets.
    line_color=[], line_width=[],
))
src_pointing_handle = ColumnDataSource(data=dict(x=[], y=[]))

# Canvas size adapts to the browser window: the figure stretches both
# axes to fill the centre cell of the layout, while `match_aspect=True`
# locks the DATA aspect to 1:1 (image pixels are square on screen,
# letterboxed if the canvas aspect ≠ image aspect).
#
# We keep the *initial* width/height as a hint for the layout engine,
# but the actual size is set by `sizing_mode='stretch_both'` against
# the parent row's stretch_both. On a laptop (1280–1440 px wide) the
# canvas shrinks so the sidebars stay visible without horizontal
# scrolling; on a 24"+ monitor it grows up to the row's natural size.
#
# Sidebar / help-panel are fixed-width — the canvas absorbs the rest.
# (FIG_W_HINT / FIG_H_HINT / SIDEBAR_W / HELPPANEL_W are defined earlier
# in the file because help_panel references HELPPANEL_W during its
# construction before this point.)
fig = figure(
    width=FIG_W_HINT, height=FIG_H_HINT,
    sizing_mode="stretch_both",
    match_aspect=True,
    # IMPORTANT: leave x_range / y_range at the default DataRange1d.
    # Per Bokeh docs match_aspect=True only works with DataRange1d —
    # switching to explicit Range1d silently breaks the aspect lock and
    # rectangular images (e.g. 2200×2500 FITS) get stretched horizontally
    # by ~factor canvas_aspect/data_aspect. The earlier Range1d "fix to
    # prevent auto-zoom-out after a click" is reverted; that potential
    # zoom-out is harmless in practice (the spec-overlap polygons stay
    # within the already-visible MSA footprint at the relevant zooms).
    #
    # No "tap" tool: it auto-selects clicked glyphs, which causes Bokeh's
    # default nonselection-rendering to fade every *other* open shutter
    # to 20% alpha. Mouse clicks come via fig.on_event(Tap, on_tap).
    tools="pan,box_zoom,reset,save",
    output_backend="webgl",
    title=None,  # vertical-space saver; identity lives in the top status bar.
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
    line_color="silver", line_alpha=0.20, line_width=0.75,
    fill_alpha=0.0,
)
stuck_open_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_stuck_open,
    # Thicker, darker outline to distinguish from user-opened shutters
    # (which are red with 1.5 px lines). Stuck-open is permanent state,
    # not a user pick — the chunkier border makes that obvious.
    line_color="#b30000", line_alpha=1.0, line_width=2.5,
    fill_color="#ff2222", fill_alpha=0.15,
)
# Spectral-overlap shutters: any operable shutter in the same s row of the
# same quadrant as a currently-open shutter. If opened, their dispersed
# spectra would overlap on the detector (MPT-style spectral conflict).
spec_overlap_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap,
    # No edge — overlap shutters are fill-only, alpha 0.1 per conflict so
    # the colour intensifies where multiple open shutters' spectra
    # contribute (alpha compositing stacks).
    line_alpha=0.0, line_width=0,
    fill_color="orange", fill_alpha=0.10,
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
    size=10, marker="circle",
    line_color="line_color",   # field-driven: yellow normally, green when matched
    line_width="line_width",
    fill_alpha=0.0,
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

# Tooltips: single-line HTML strings. The Bokeh wrapper `.bk-tooltip`
# is also tightened via CSS in templates/index.html so the popup is a
# thin pill rather than a paragraph block.
_TIP_BASE_STYLE = (
    "font-family: Calibri, Helvetica, Arial, sans-serif; "
    "font-size: 11px; line-height: 1.25; padding: 0 2px; "
    "white-space: nowrap;"
)
fig.add_tools(HoverTool(
    renderers=[bg_shutters_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#555;">'
        f'  Q@q · s@s · d@d &nbsp;<span style="color:#aaa">operable</span>'
        f'</div>'
    ),
))
fig.add_tools(HoverTool(
    renderers=[stuck_open_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#b30000; font-weight:600;">'
        f'  Q@q · s@s · d@d &nbsp;<span style="color:#7a0000">STUCK OPEN</span>'
        f'</div>'
    ),
))
fig.add_tools(HoverTool(
    renderers=[open_shutters_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#1a3b66;">'
        f'<b>Q@q s@s d@d</b>'
        f'<span style="color:#888"> · </span>@target'
        f'<span style="color:#888"> · </span>'
        f'<b>@lam_blue{{0.00}}</b><span style="color:#888">–</span>'
        f'<b>@lam_red{{0.00}}</b><span style="color:#888"> μm · gap </span>'
        f'<b>@gap_label</b>'
        f'</div>'
    ),
))
fig.add_tools(HoverTool(
    renderers=[fixed_slits_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#8a6300;">'
        f'  <b>@name</b>'
        f'</div>'
    ),
))
fig.add_tools(HoverTool(
    renderers=[target_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#1a3b66;">'
        f'  <b>@id</b>'
        f'  <span style="color:#888;"> · </span>@ra{{0.0000}}, @dec{{0.0000}}'
        f'  <span style="color:#888;"> · Pr </span>@pr'
        f'</div>'
    ),
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


def _loading_overlay_html(msg: str) -> str:
    """HTML for the full-page loading overlay (escapes the Bokeh layout via
    position: fixed, so it covers the whole viewport)."""
    safe_msg = (msg or "Loading…").replace("<", "&lt;").replace(">", "&gt;")
    return (
        "<style>"
        "@keyframes vmpt-spin { to { transform: rotate(360deg); } }"
        "@keyframes vmpt-fadein { from { opacity: 0; } to { opacity: 1; } }"
        "</style>"
        '<div style="position: fixed; top: 0; left: 0; right: 0; bottom: 0;'
        "            z-index: 9999;"
        "            display: flex; flex-direction: column;"
        "            align-items: center; justify-content: center;"
        "            background: rgba(8, 16, 32, 0.55);"
        "            backdrop-filter: blur(2px);"
        "            -webkit-backdrop-filter: blur(2px);"
        '            animation: vmpt-fadein 180ms ease-out;">'
        '  <div style="width: 84px; height: 84px;'
        "             border: 7px solid rgba(255,255,255,0.18);"
        "             border-top-color: #FFB400;"
        "             border-right-color: #FFB400;"
        "             border-radius: 50%;"
        '             animation: vmpt-spin 0.9s linear infinite;"></div>'
        '  <div style="margin-top: 22px;'
        "             color: white;"
        "             font-family: Calibri, Helvetica, Arial, sans-serif;"
        "             font-size: 16px;"
        "             font-weight: 600;"
        "             letter-spacing: 0.3px;"
        '             text-shadow: 0 1px 4px rgba(0,0,0,0.5);">'
        f"    {safe_msg}"
        "  </div>"
        "</div>"
    )


def _show_loading(msg: str) -> None:
    """Show the full-page spinner overlay. Pair with _hide_loading().

    Schedules a 60-second safety timeout so the overlay never gets stuck
    if a callback path fails to call _hide_loading.
    """
    loading_banner.text = _loading_overlay_html(msg)
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
        return dict(xs=[], ys=[], q=[], s=[], d=[], target=[],
                    lam_blue=[], lam_red=[], lam_gap_lo=[], lam_gap_hi=[],
                    gap_label=[])
    keys = list(state["open_shutters"].keys())
    xs, ys, qs, ss, ds = _polygons_for_shutter_keys(keys, pa_v3, fid_pix, jinv)
    tgt: list[str] = []
    lam_b: list[float] = []
    lam_r: list[float] = []
    lam_glo: list[float] = []
    lam_ghi: list[float] = []
    gap_lbl: list[str] = []
    for (q, s, d) in keys:
        sh = state["open_shutters"][(q, s, d)]
        tgt.append(str(sh.target_id) if sh.target_id is not None else "")
        v2c = float(V2_MSA[q - 1, s - 1, d - 1])
        v3c = float(V3_MSA[q - 1, s - 1, d - 1])
        cut = cutoffs(v2c, v3c, state["disperser"], state["filter"])
        b = cut.get("lam_blue"); r = cut.get("lam_red")
        glo = cut.get("lam_gap_lo"); ghi = cut.get("lam_gap_hi")
        lam_b.append(b if b is not None else float("nan"))
        lam_r.append(r if r is not None else float("nan"))
        lam_glo.append(glo if glo is not None else float("nan"))
        lam_ghi.append(ghi if ghi is not None else float("nan"))
        if glo is None or ghi is None:
            gap_lbl.append("(no gap on this spectrum)")
        else:
            gap_lbl.append(f"{glo:.2f} – {ghi:.2f} μm")
    return dict(xs=xs, ys=ys, q=qs, s=ss, d=ds, target=tgt,
                lam_blue=lam_b, lam_red=lam_r,
                lam_gap_lo=lam_glo, lam_gap_hi=lam_ghi, gap_label=gap_lbl)


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

    # ── Operable-shutter (silver-edge) layer is built AFTER spec-overlap
    # below, because it filters out user-opens + spec-overlap + stuck-open
    # so the silver edge cleanly highlights the *unaffected, ready-to-pick*
    # shutters. We compute the spec-overlap set first, then this layer.

    # Spectral-overlap: operable shutters whose spectrum collides with any
    # dispersed shutter's spectrum on the detector. Two shutters share a
    # detector y-row iff (a) they sit on the SAME detector half (Q1/Q3 →
    # NRS1, Q2/Q4 → NRS2; the MSA's s-axis tiles into a single column on
    # each detector half, see eMPT's CSV layout) and (b) their s indices
    # are within SHVAL_S_TOLERANCE. They share dispersed pixels iff their
    # V2 separation is below v2_overlap_distance(disperser, filter).
    #
    # Cross-quadrant pairings other than Q1↔Q3 and Q2↔Q4 image onto
    # *different* detectors and therefore never overlap, even when their
    # V2 windows would otherwise pass the distance check (relevant for the
    # H gratings, whose V2 half-extent is ~500″).
    #
    # Sources of dispersion: every user-opened shutter AND every shutter
    # known to be stuck open (REASON == 2). Stuck-opens always disperse
    # light onto the detector, so their spec-overlap rows must light up
    # even when the user hasn't picked them.
    NRS1_QUADS = {1, 3}
    NRS2_QUADS = {2, 4}
    # |Δs| ≤ 1 → only the open shutter's row and its immediate neighbour
    # rows above/below disperse onto overlapping detector pixels. Wider
    # tolerances over-paint inconvenient (eMPT uses shval ≈ s exactly).
    SHVAL_S_TOLERANCE = 1
    user_opens = list(state["open_shutters"].keys())
    stuck_flat = np.where(_FLAT_REASON == 2)[0]
    stuck_keys = [
        (int(f // (171 * 365)) + 1,
         int((f % (171 * 365)) // 365) + 1,
         int(f % 365) + 1)
        for f in stuck_flat
    ]
    dispersion_sources = user_opens + stuck_keys
    if dispersion_sources:
        v2_overlap = float(v2_overlap_distance(state["disperser"], state["filter"]))
        s_arr = (np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) % (171 * 365)) // 365
        q_arr = np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) // (171 * 365) + 1
        in_view_op = in_view & (_FLAT_REASON == 0)
        chunks: list[np.ndarray] = []
        for q_o, s_o, d_o in dispersion_sources:
            open_flat = (q_o - 1) * 171 * 365 + (s_o - 1) * 365 + (d_o - 1)
            v2_o = float(_V2_OFFSETS_ALL[open_flat] + MSA_V2_REF)
            # Only the matching detector half can share a y-row.
            partners = NRS1_QUADS if q_o in NRS1_QUADS else NRS2_QUADS
            same_det = np.isin(q_arr, list(partners))
            same_row = np.abs(s_arr - (s_o - 1)) <= SHVAL_S_TOLERANCE
            near_v2 = np.abs((_V2_OFFSETS_ALL + MSA_V2_REF) - v2_o) < v2_overlap
            idx_this = np.where(in_view_op & same_det & same_row & near_v2)[0]
            idx_this = idx_this[idx_this != open_flat]
            if idx_this.size:
                chunks.append(idx_this)
        overlap_idx = (
            np.unique(np.concatenate(chunks)) if chunks else np.empty(0, dtype=np.int64)
        )
    else:
        overlap_idx = np.empty(0, dtype=np.int64)
    src_spec_overlap.data = _project_indices_to_cds(overlap_idx, pa_v3, fid_pix, jinv)

    # ── Unaffected operable (silver-edge) layer.
    # Show only shutters that are operable, in view, NOT currently open
    # (red), NOT stuck-open (dark red), and NOT in the spec-overlap set
    # (orange). The silver edges then act as a "click here, ready" hint
    # rather than overlapping the heavier coloured layers underneath.
    show_shutters = 1 in layers_box.active
    if show_shutters:
        op_mask = in_view & (_FLAT_REASON == 0)
        if state["open_shutters"]:
            open_flat = np.array([
                (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
                for (q, s, d) in state["open_shutters"].keys()
            ], dtype=np.int64)
            op_mask[open_flat] = False
        if overlap_idx.size:
            op_mask[overlap_idx] = False
        op_idx = np.where(op_mask)[0]
        if op_idx.size > MAX_OPERABLE_RENDER:
            # Above the cap: blank rather than stride-sample (a sparse
            # grid looks broken). User zooms in to see them.
            src_bg_shutters.data = dict(xs=[], ys=[], q=[], s=[], d=[])
        else:
            src_bg_shutters.data = _project_indices_to_cds(op_idx, pa_v3, fid_pix, jinv)
    else:
        src_bg_shutters.data = dict(xs=[], ys=[], q=[], s=[], d=[])

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
        ids = [str(i) for i in np.asarray(cat.ids)[mask]]
        # Highlight catalog entries that already correspond to an open
        # shutter's target_id: green ring, thicker line.
        matched_target_ids = {
            str(sh.target_id) for sh in state["open_shutters"].values()
            if sh.target_id is not None
        }
        line_colors = [
            "#2e9b3f" if tid in matched_target_ids else "yellow"
            for tid in ids
        ]
        line_widths = [
            2.5 if tid in matched_target_ids else 1.5
            for tid in ids
        ]
        src_targets.data = dict(
            x=x[mask].tolist(),
            y=y[mask].tolist(),
            id=ids,
            ra=cat.ra_deg[mask].tolist(),
            dec=cat.dec_deg[mask].tolist(),
            pr=cat.priority[mask].tolist(),
            line_color=line_colors,
            line_width=line_widths,
        )
    else:
        src_targets.data = dict(
            x=[], y=[], id=[], ra=[], dec=[], pr=[],
            line_color=[], line_width=[],
        )

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
    # Wide single-row status bar above the figure. Each cell is a
    # labelled key/value chip; spec-conflict count is colour-coded
    # (green=0, yellow=<100, orange=<1000, red=>=1000).
    if n_overlap == 0:
        oc = "#16803c"   # green
    elif n_overlap < 100:
        oc = "#9a7400"   # amber
    elif n_overlap < 1000:
        oc = "#cc6a00"   # orange
    else:
        oc = "#b3261e"   # red
    img_label = state["image"].source_path.split("/")[-1] if state["image"] else "—"
    cell = (
        "display:inline-block; padding:0 12px; "
        "border-right:1px solid #c2d6f0; line-height:32px; vertical-align:middle;"
    )
    label_style = "color:#5a6b85; font-size:11px; text-transform:uppercase; letter-spacing:0.5px;"
    val_style   = "color:#1a3b66; font-weight:600; font-size:14px; margin-left:4px;"
    stats_div.text = (
        f'<div style="display:flex; flex-wrap:wrap; align-items:center; gap:0;">'
        f'  <span style="{cell}">'
        f'    <span style="{label_style}">Image</span>'
        f'    <span style="{val_style}">{img_label}</span>'
        f'  </span>'
        f'  <span style="{cell}">'
        f'    <span style="{label_style}">RA · Dec</span>'
        f'    <span style="{val_style}">{fiducial.ra.deg:.5f} · {fiducial.dec.deg:.5f}</span>'
        f'    <span style="color:#7c8aa0; font-size:11px; margin-left:6px;">'
        f'      ({ra_hms} · {dec_dms})'
        f'    </span>'
        f'  </span>'
        f'  <span style="{cell}">'
        f'    <span style="{label_style}">V3 PA</span>'
        f'    <span style="{val_style}">{pa_v3:.2f}°</span>'
        f'    <span style="{label_style}; margin-left:10px;">APA</span>'
        f'    <span style="{val_style}">{apa:.2f}°</span>'
        f'  </span>'
        f'  <span style="{cell}">'
        f'    <span style="{label_style}">Disperser</span>'
        f'    <span style="{val_style}">{state["disperser"]} / {state["filter"]}</span>'
        f'  </span>'
        f'  <span style="{cell}">'
        f'    <span style="{label_style}">Open</span>'
        f'    <span style="{val_style}">{n_op}</span>'
        f'    <span style="color:#7c8aa0; font-size:11px; margin-left:4px;">'
        f'      across {n_tgt_open} target{"s" if n_tgt_open != 1 else ""}'
        f'    </span>'
        f'  </span>'
        f'  <span style="{cell} border-right:none;">'
        f'    <span style="{label_style}">Conflicts</span>'
        f'    <span style="color:{oc}; font-weight:700; font-size:14px; margin-left:4px;">{n_overlap}</span>'
        f'    <span style="color:#7c8aa0; font-size:11px; margin-left:4px;">shutters</span>'
        f'  </span>'
        f'</div>'
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
    # Re-fit figure ranges to the new image's pixel extent. We use a
    # single atomic update per axis (via `.update()`) so match_aspect's
    # aspect-locking pass sees a coherent change and re-evaluates the
    # constraint. Setting `.start` and `.end` separately on DataRange1d
    # can leave the ranges in an intermediate state where match_aspect
    # doesn't fire — that was the "switch from A370 to a JPG and the
    # aspect comes out wrong" bug.
    fig.x_range.update(start=0, end=W)
    fig.y_range.update(start=0, end=H)
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


def _set_image_and_recenter(
    img: LoadedImage, source_label: str, force_recenter: bool = False,
) -> None:
    """Install a freshly-loaded image as the active canvas.

    `force_recenter=True` (used by the "Load Abell 370 / RXCJ0600 example"
    buttons) overrides the preserve-existing-pointing guard. The guard
    exists so that loading an image AFTER an APT plan doesn't clobber
    the plan's pointing; but when the user is explicitly switching
    examples, they expect the pointing to follow. We also auto-recenter
    if the existing pointing is more than 30 arcmin away from the new
    image (then no MSA would land on the image anyway).
    """
    state["image"] = img
    H, W = img.shape[:2]
    has_pointing = bool(
        (ra_input.value or "").strip() and (dec_input.value or "").strip()
    )
    auto_recentered = False
    try:
        center = img.wcs.pixel_to_world(W / 2, H / 2)
    except Exception:  # noqa: BLE001
        center = None
    if has_pointing and not force_recenter and center is not None:
        # Detect "pointing is on a different sky" — if the current
        # pointing is more than 30' from the new image's centre, no MSA
        # would land on it. Recenter rather than leave the user staring
        # at an image with no overlay.
        try:
            cur = SkyCoord(float(ra_input.value), float(dec_input.value),
                           unit=u.deg, frame="icrs")
            if cur.separation(center).to_value(u.arcmin) > 30.0:
                force_recenter = True
                auto_recentered = True
        except (ValueError, TypeError):
            pass
    if (force_recenter or not has_pointing) and center is not None:
        try:
            ra_input.value = f"{center.ra.deg:.6f}"
            dec_input.value = f"{center.dec.deg:.6f}"
            state["ra_deg"] = center.ra.deg
            state["dec_deg"] = center.dec.deg
        except Exception:  # noqa: BLE001
            pass
    refresh_image_glyph()
    refresh_overlays()
    if force_recenter and auto_recentered:
        suffix = " — auto-recentered (previous pointing was outside this field)"
    elif force_recenter:
        suffix = ""
    elif has_pointing:
        suffix = " (kept existing pointing)"
    else:
        suffix = ""
    _set_status(f"Loaded {source_label} ({W}×{H}).{suffix}", "ok")


def _load_fits_from_path(path: str, force_recenter: bool = False) -> None:
    try:
        img = load_fits(path)
        _set_image_and_recenter(
            img, f"FITS {Path(path).name}", force_recenter=force_recenter,
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"FITS load failed: {e}", "err")
        traceback.print_exc()
    finally:
        _hide_loading()


def _load_jpg_pair_from_paths(
    jpg_path: str, sidecar_path: str, force_recenter: bool = False,
) -> None:
    try:
        img = load_jpg_with_sidecar(jpg_path, sidecar_path, max_dim=6000)
        _set_image_and_recenter(
            img, f"JPG+sidecar {Path(jpg_path).name}",
            force_recenter=force_recenter,
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"JPG+sidecar load failed: {e}", "err")
        traceback.print_exc()
    finally:
        _hide_loading()


def _load_catalog_from_path(path: str) -> None:
    try:
        cat = load_catalog(path)
        state["catalog"] = cat
        _rebuild_shutter_catalog_index()
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


# ---------------------------------------------------------------------------
# Callbacks: pointing & display
# ---------------------------------------------------------------------------


def on_pointing(attr, old, new):
    try:
        state["ra_deg"] = float(ra_input.value)
        state["dec_deg"] = float(dec_input.value)
    except (TypeError, ValueError):
        return
    _show_loading("Updating pointing…")
    def _do():
        try:
            _rebuild_shutter_catalog_index()
            refresh_overlays()
        finally:
            _hide_loading()
    _deferred(_do)


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
    _show_loading("Updating V3 PA…")
    def _do():
        try:
            _rebuild_shutter_catalog_index()
            refresh_overlays()
        finally:
            _hide_loading()
    _deferred(_do)


def on_v3pa_text(attr, old, new):
    if state.get("_syncing_pa"):
        return
    try:
        v = float(v3pa_input.value)
    except (TypeError, ValueError):
        return
    _sync_pa_widgets(v, source="v3pa_text")
    _show_loading("Updating V3 PA…")
    def _do():
        try:
            _rebuild_shutter_catalog_index()
            refresh_overlays()
        finally:
            _hide_loading()
    _deferred(_do)


def on_apa_text(attr, old, new):
    if state.get("_syncing_pa"):
        return
    try:
        apa = float(apa_input.value)
    except (TypeError, ValueError):
        return
    _sync_pa_widgets(apa - V3_IDL_Y_ANGLE, source="apa_text")
    _show_loading("Updating APA…")
    def _do():
        try:
            _rebuild_shutter_catalog_index()
            refresh_overlays()
        finally:
            _hide_loading()
    _deferred(_do)


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
    _show_loading(f"Recomputing for {d_part} / {f_part}…")
    def _do():
        try:
            refresh_overlays()
        finally:
            _hide_loading()
    _deferred(_do)


def on_slitlet_height(attr, old, new):
    state["slitlet_height"] = int(slitlet_select.value)


def on_snap(attr, old, new):
    state["snap_to_operable"] = 0 in snap_box.active


# ---------------------------------------------------------------------------
# Hand-picking: tap callbacks and slitlet builder
# ---------------------------------------------------------------------------


def _v2v3_to_radec(v2: float, v3: float, fiducial: SkyCoord, pa_v3: float) -> tuple[float, float]:
    """Forward map a single (V2, V3) arcsec coord onto (RA, Dec) degrees."""
    offsets = np.array([[v2 - MSA_V2_REF, v3 - MSA_V3_REF]])
    rotated = np.dot(offsets, rot_matrix(pa_v3))
    sky = fiducial.spherical_offsets_by(
        rotated[0, 0] * u.arcsec, rotated[0, 1] * u.arcsec
    )
    return float(sky.ra.deg), float(sky.dec.deg)


def _sky_to_v2v3(sky: SkyCoord, fiducial: SkyCoord, pa_v3: float) -> tuple[float, float]:
    """Inverse of v2v3_to_radec: returns (V2, V3) in arcsec for a single sky point."""
    d_lon, d_lat = fiducial.spherical_offsets_to(sky)
    dx = d_lon.to_value(u.arcsec)
    dy = d_lat.to_value(u.arcsec)
    # Inverse of dot(offsets, rot_matrix(pa_v3)) is dot(offsets, rot_matrix(-pa_v3))
    rot = rot_matrix(-pa_v3)
    v2_off, v3_off = float(dx * rot[0, 0] + dy * rot[1, 0]), float(dx * rot[0, 1] + dy * rot[1, 1])
    return v2_off + MSA_V2_REF, v3_off + MSA_V3_REF


def _rebuild_shutter_catalog_index() -> None:
    """Vectorised: for every catalog source visible at the current pointing
    + PA, find the nearest operable shutter and bucket the source id under
    that shutter. The result lives in state['shutter_to_catids'] as
    {(q,s,d) → [source_id, …]}. Cleared and rebuilt whenever the
    pointing, PA, or catalog changes.

    Matching uses APT's "Unconstrained" Source Centering rule
    (https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/
    nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-planner):
    a source still matches the shutter even if its centre falls behind
    the opaque bars separating shutters. So the matching cell is the
    full MSA **pitch** (≈0.27″ × 0.53″) — the V2/V3 Voronoi cell of the
    shutter — not just the narrower open aperture (0.20″ × 0.46″).
    """
    state["shutter_to_catids"] = {}
    cat = state.get("catalog")
    fiducial = _pointing_skycoord()
    if cat is None or fiducial is None:
        return
    if len(cat.ra_deg) == 0:
        return
    pa_v3 = state["pa_v3"]
    sky = SkyCoord(cat.ra_deg, cat.dec_deg, unit=u.deg, frame="icrs")
    d_lon, d_lat = fiducial.spherical_offsets_to(sky)
    dx = d_lon.to_value(u.arcsec)
    dy = d_lat.to_value(u.arcsec)
    rot = rot_matrix(-pa_v3)
    v2_arr = dx * rot[0, 0] + dy * rot[1, 0] + MSA_V2_REF
    v3_arr = dx * rot[0, 1] + dy * rot[1, 1] + MSA_V3_REF
    # Half-pitch in V2 / V3 (the full shutter "cell" including the bars).
    # MSA pitch is ≈0.27″ × 0.53″; half = ≈0.135″ × ≈0.265″.
    SHUTTER_HALF_PITCH_V2 = 0.135
    SHUTTER_HALF_PITCH_V3 = 0.265
    bucket: dict[tuple[int, int, int], list] = {}
    for i in range(len(v2_arr)):
        dv2 = V2_MSA - v2_arr[i]
        dv3 = V3_MSA - v3_arr[i]
        d2 = dv2 * dv2 + dv3 * dv3
        idx = int(np.argmin(d2))
        q = idx // (171 * 365)
        rem = idx % (171 * 365)
        s = rem // 365
        d = rem % 365
        if not OPERABLE[q, s, d]:
            continue
        # "Unconstrained" footprint check: the source must lie within the
        # full shutter pitch (i.e. closer to this shutter centre than to
        # any neighbouring shutter). Reject anything outside the pitch
        # box — that means the source is entirely outside the MSA grid,
        # not just sitting on a bar.
        if (abs(dv2.flat[idx]) > SHUTTER_HALF_PITCH_V2
                or abs(dv3.flat[idx]) > SHUTTER_HALF_PITCH_V3):
            continue
        key = (q + 1, s + 1, d + 1)
        bucket.setdefault(key, []).append(cat.ids[i])
    state["shutter_to_catids"] = bucket


def _shutter_source_id(q: int, s: int, d: int) -> str | None:
    """Return the catalog source id (as a string) that falls inside this
    shutter, or None if none does. If multiple sources land in the same
    shutter we return the first one (catalog order)."""
    bucket = state.get("shutter_to_catids") or {}
    ids = bucket.get((int(q), int(s), int(d)))
    if not ids:
        return None
    return str(ids[0])


def _nearest_shutter(v2_target: float, v3_target: float,
                     require_operable: bool = True,
                     max_dist_arcsec: float | None = None) -> tuple[int, int, int] | None:
    """Brute-force argmin over all shutters of squared V2/V3 distance.

    Returns None if the closest matching shutter is farther than
    `max_dist_arcsec` from (v2_target, v3_target). Useful so a click in
    empty sky doesn't open a faraway shutter.
    """
    dv2 = V2_MSA - v2_target
    dv3 = V3_MSA - v3_target
    d2 = dv2 * dv2 + dv3 * dv3
    if require_operable:
        d2 = np.where(OPERABLE, d2, np.inf)
    idx = int(np.argmin(d2))
    min_d2 = float(d2.flat[idx])
    if not np.isfinite(min_d2):
        return None
    if max_dist_arcsec is not None and min_d2 > max_dist_arcsec * max_dist_arcsec:
        return None
    q = idx // (171 * 365)
    rem = idx % (171 * 365)
    s = rem // 365
    d = rem % 365
    return (q + 1, s + 1, d + 1)


def _slitlet_offsets(n: int) -> list[int]:
    """Return the list of s-offsets (relative to the clicked shutter) for an
    N-shutter slitlet:
      • N=1 → [0]
      • N=2 → [-1, 0]      (clicked shutter + one row lower y in detector)
      • N=3 → [-1, 0, +1]  (centred)
      • N=5 → [-2,-1,0,+1,+2]
    Anything else falls back to a centred (or near-centred) layout.
    """
    if n <= 1:
        return [0]
    if n == 2:
        return [-1, 0]
    half = n // 2
    return list(range(-half, n - half))


def _add_slitlet(q: int, s_click: int, d: int, target_id: str | None) -> int:
    """Open an N-shutter slitlet at column (q, d) anchored on the clicked
    shutter s=s_click. N comes from state['slitlet_height']. The clicked
    shutter is always opened (and marked 'target' if it's an offset of 0
    in the layout from `_slitlet_offsets`); siblings are 'sky'.

    If `target_id` is None we check whether any catalog source falls inside
    each opened shutter's footprint and adopt it as the shutter's target
    id (so source-IDs propagate automatically into the APT export). All
    shutters in this slitlet inherit the same id (the first one found
    when scanning the layout).

    Returns the number of shutters added (operable ones; failed-closed are skipped).
    """
    n = int(state["slitlet_height"])
    offsets = _slitlet_offsets(n)
    # If the caller didn't already know a target_id, pick the first
    # catalog source landing inside any opened shutter of this slitlet.
    if target_id is None:
        for offset in offsets:
            s_try = s_click + offset
            if not (1 <= s_try <= 171):
                continue
            if not OPERABLE[q - 1, s_try - 1, d - 1]:
                continue
            tid = _shutter_source_id(q, s_try, d)
            if tid is not None:
                target_id = tid
                break
    added = 0
    for offset in offsets:
        s = s_click + offset
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
    """Open an N-shutter slitlet anchored on (q, s, d) — N comes from
    state['slitlet_height']. If the clicked shutter is already open we
    remove it plus its slitlet siblings instead.

    `target_id` may be:
      • a known catalog id (e.g. user clicked a target marker)
      • None → `_add_slitlet` looks up a catalog source whose footprint
        falls inside any of the opened shutters and tags the slitlet
        with that id automatically.
    """
    key = (q, s, d)
    _push_history()
    if key in state["open_shutters"]:
        sh = state["open_shutters"].pop(key)
        # Remove slitlet siblings — siblings share (q, d) and a small Δs.
        # Match by either target_id (if present on both) or position
        # alone (so manual single shutters and their N-shutter siblings
        # still come down together).
        n = int(state["slitlet_height"])
        half = max(1, n)  # close up to N rows around the click (≥ 1)
        for k in list(state["open_shutters"].keys()):
            other = state["open_shutters"][k]
            if other.q != q or other.d != d:
                continue
            if abs(other.s - s) > half:
                continue
            if other.target_id != sh.target_id:
                # Don't sweep up a different target's slitlet that happens
                # to sit nearby — only same-target (or both None) siblings.
                continue
            del state["open_shutters"][k]
        _set_status(f"Closed shutter ({q},{s},{d}) and slitlet siblings.", "ok")
    else:
        n_added = _add_slitlet(q, s, d, target_id=target_id)
        # Re-read the freshly opened shutter to surface the auto-matched
        # target id (if any) in the status line.
        opened = state["open_shutters"].get(key)
        auto_id = (opened.target_id if opened else None)
        if auto_id and auto_id != target_id:
            _set_status(
                f"Opened {n_added}-shutter slitlet at ({q},{s},{d}); "
                f"matched catalog source {auto_id}.", "ok"
            )
        elif target_id is not None:
            _set_status(
                f"Target {target_id} → slitlet ({q},{s},{d}), "
                f"{n_added} shutters opened.", "ok"
            )
        else:
            _set_status(
                f"Opened {n_added}-shutter slitlet at ({q},{s},{d}).", "ok"
            )
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

    # 2) Otherwise: snap to nearest shutter and toggle. Clicks in empty
    # sky (far from any shutter) or on stuck-open / failed-closed shutters
    # do nothing — silent — so the user can pan and click around without
    # accidentally opening faraway shutters.
    try:
        sky = img.wcs.pixel_to_world(x_data, y_data)
        v2, v3 = _sky_to_v2v3(sky, fiducial, state["pa_v3"])
        # Threshold = ~half a shutter diagonal. Shutters are 0.20" wide x
        # 0.46" tall in V2/V3; half-diagonal ~ 0.25". Use 0.5" to be a bit
        # forgiving on aim.
        nearest = _nearest_shutter(v2, v3, require_operable=False,
                                   max_dist_arcsec=0.5)
    except Exception:  # noqa: BLE001
        return
    if nearest is None:
        return  # clicked far from any shutter
    q, s, d = nearest
    # If the user is closing an already-open shutter, that's allowed even
    # if the shutter is somehow flagged non-operable (shouldn't happen,
    # but be defensive). Otherwise refuse to open non-operable shutters.
    if (q, s, d) not in state["open_shutters"] and not OPERABLE[q - 1, s - 1, d - 1]:
        return  # clicked on a stuck-open or failed-closed shutter

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

    # ─── Build the observed_targets list ───────────────────────────────────
    # One row per "real" target_id picked up from the catalog (looked up by
    # shutter-footprint at pick time, or re-checked here so newly-loaded
    # catalogs are picked up). For each contiguous slitlet of open shutters
    # in the same (q, d) column that has NO real catalog source, fake one
    # entry positioned at the centre shutter's sky location. Stuck-open
    # shutters never get faked entries.
    cat = state["catalog"]
    targets_rows = []
    real_ids_seen: set[str] = set()
    used_target_nos: set[int] = set()

    def _cat_lookup(tid: str) -> tuple[float, float, int, str] | None:
        """Look up a catalog row by id. Returns (ra, dec, weight, label).
        `label` is the catalog's `label`/`name` column value when
        available, else the literal string "real" (so the output
        catalog's Label column always distinguishes real vs synth)."""
        if cat is None:
            return None
        ids_str = [str(i) for i in cat.ids]
        if tid not in ids_str:
            return None
        k = ids_str.index(tid)
        pr = int(cat.priority[k]) if np.isfinite(cat.priority[k]) else 1
        label_val = ""
        try:
            label_val = str(cat.label[k]) if cat.label is not None else ""
        except (AttributeError, IndexError, TypeError):
            label_val = ""
        label_val = label_val.strip()
        if not label_val:
            label_val = "real"
        return float(cat.ra_deg[k]), float(cat.dec_deg[k]), pr, label_val

    def _push_target_row(
        tid: str | None, ra_d: float, dec_d: float, pr: int, label: str,
    ) -> int:
        """Append a row and return the assigned No_cat."""
        try:
            target_no = int(tid) if tid is not None else None
        except (ValueError, TypeError):
            target_no = None
        if target_no is None or target_no in used_target_nos:
            # Generate a fresh sequential number that doesn't collide.
            target_no = max(used_target_nos, default=0) + 1
            while target_no in used_target_nos:
                target_no += 1
        used_target_nos.add(target_no)
        targets_rows.append({
            "No_cat": target_no,
            "Pr": pr,
            "ra_deg": ra_d,
            "dec_deg": dec_d,
            "label": label,
        })
        return target_no

    # Step 1: real catalog sources tied to user picks. For each open
    # user-shutter that has a target_id OR sits inside a catalog source's
    # footprint, register the source.
    for (q, s, d), sh in state["open_shutters"].items():
        tid = sh.target_id or _shutter_source_id(q, s, d)
        if not tid or tid in real_ids_seen:
            continue
        info = _cat_lookup(str(tid))
        if info is None:
            continue  # tid we don't know — leave to be re-faked from geometry
        ra_d, dec_d, pr, label_val = info
        real_ids_seen.add(str(tid))
        _push_target_row(str(tid), ra_d, dec_d, pr, label=label_val)

    # Step 2: group open shutters into per-(q,d) consecutive-s runs (slitlets)
    # and fake an entry for any run that has no real source attached.
    fiducial = _pointing_skycoord()
    col_to_s: dict[tuple[int, int], list[tuple[int, str | None]]] = {}
    for (q, s, d), sh in state["open_shutters"].items():
        col_to_s.setdefault((q, d), []).append((s, sh.target_id))
    for (q, d), rows in col_to_s.items():
        rows.sort(key=lambda t: t[0])
        # Split into consecutive-s runs
        runs: list[list[tuple[int, str | None]]] = []
        run: list[tuple[int, str | None]] = []
        for s, tid in rows:
            if run and s == run[-1][0] + 1:
                run.append((s, tid))
            else:
                if run:
                    runs.append(run)
                run = [(s, tid)]
        if run:
            runs.append(run)
        # For each run, check if it already has a real catalog source.
        for run in runs:
            has_real = any(
                (tid and str(tid) in real_ids_seen) or
                (_shutter_source_id(q, s, d) and
                 str(_shutter_source_id(q, s, d)) in real_ids_seen)
                for s, tid in run
            )
            if has_real:
                continue
            # Fake an entry at the middle shutter's RA/Dec; mark it as
            # synthesized so the output catalog's Label column tells the
            # user (and APT) which rows weren't in their input catalog.
            mid = run[len(run) // 2]
            s_mid = mid[0]
            v2 = float(V2_MSA[q - 1, s_mid - 1, d - 1])
            v3 = float(V3_MSA[q - 1, s_mid - 1, d - 1])
            if fiducial is not None:
                ra_d, dec_d = _v2v3_to_radec(v2, v3, fiducial, state["pa_v3"])
            else:
                ra_d, dec_d = float("nan"), float("nan")
            no_cat = _push_target_row(
                None, ra_d, dec_d, pr=5, label="vMPT_synth",
            )
            # Tag every shutter in the run with this fake id so later
            # exporters (MPT plan primaryIds) see consistent target IDs.
            fake_id = str(no_cat)
            for s, _ in run:
                cur = state["open_shutters"].get((q, s, d))
                if cur is not None and (cur.target_id is None or cur.target_id == ""):
                    state["open_shutters"][(q, s, d)] = OpenShutter(
                        q=cur.q, s=cur.s, d=cur.d,
                        target_id=fake_id,
                        role=cur.role if cur.role != "manual" else "target" if s == s_mid else "sky",
                    )

    pa_v3 = state["pa_v3"]
    # PA_V3 - PA_AP = -V3IdlYAngle (mod 360); for NRS_FULL_MSA V3IdlYAngle ~ 138.5746°.
    pa_ap = (pa_v3 + V3_IDL_Y_ANGLE) % 360.0
    pointing = Pointing(
        ra_deg=float(ra_input.value),
        dec_deg=float(dec_input.value),
        apa_v3_deg=pa_v3,
        pa_ap_deg=pa_ap,
    )
    # Resolve the MPT catalog filename: align the .cat basename with the
    # catalog.name we write into the plan JSON. If the user loaded a
    # specific catalog, mirror its basename; otherwise use the default.
    if cat is not None and getattr(cat, "source_path", None):
        mpt_catalog_name = Path(cat.source_path).stem + ".cat"
    else:
        mpt_catalog_name = MPT_CATALOG_FILENAME

    try:
        # 1) APT-importable primaries catalog (the file APT's Target List
        # importer wants). The plan JSON's catalog.name matches its stem.
        write_mpt_catalog(str(out_dir / mpt_catalog_name), targets_rows)
        # 2) eMPT-style outputs (observed targets, pointing, shutter mask)
        write_observed_targets_cat(str(out_dir / EMPT_OBSERVED_FILENAME), targets_rows)
        write_pointing_summary_txt(
            str(out_dir / EMPT_POINTING_FILENAME),
            pointing, state["disperser"], state["filter"],
            n_targets_total=(len(cat.ra_deg) if cat is not None else 0),
            n_targets_accepted=len(targets_rows),
        )
        open_list = list(state["open_shutters"].values())
        write_shutter_mask_csv(
            str(out_dir / EMPT_SHUTTER_MASK_FILENAME),
            open_list, OPERABLE, REASON,
        )
        # 3) MPT plan + vMPT workspace sidecar so the bundle round-trips
        # cleanly: Session → Load on either file restores everything.
        export_session_json(_build_current_session(), str(out_dir / MPT_PLAN_FILENAME))
        _set_status(
            f"Wrote bundle to {out_dir} — {MPT_PLAN_FILENAME} + "
            f"{mpt_catalog_name} (APT import) + {WORKSPACE_FILENAME} "
            f"(vMPT state) + 3 eMPT_* files. "
            f"{len(targets_rows)} targets, {len(open_list)} open shutters.",
            "ok", clear_after=18,
        )
        # Pre-fill the Session-load input so the user can re-load with one click.
        session_load_path_input.value = str(out_dir / MPT_PLAN_FILENAME)
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
    img = state["image"]
    return Session(
        pointing_ra_deg=float(state["ra_deg"]),
        pointing_dec_deg=float(state["dec_deg"]),
        pa_v3_deg=float(state["pa_v3"]),
        disperser=state["disperser"],
        filter_name=state["filter"],
        slitlet_height=int(state["slitlet_height"]),
        open_shutters=list(state["open_shutters"].values()),
        highlighted=list(state["highlighted"]),
        image_path=(img.source_path if img else None),
        wcs_sidecar_path=(getattr(img, "wcs_sidecar_path", None) if img else None),
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


def _resolve_jpg_sidecar(jpg_path: Path, recorded: Optional[str]) -> Optional[str]:
    """Find the WCS sidecar for a JPG session image. Prefer the path recorded
    in the session JSON; fall back to a sibling .fits in the same directory."""
    if recorded and Path(recorded).exists():
        return recorded
    candidates = sorted(jpg_path.parent.glob("*.fits"))
    return str(candidates[0]) if candidates else None


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

    # Apply pointing/PA/disperser/shutters FIRST so they survive even if the
    # image fails to load. We push them straight into state and the widgets
    # rather than via on_change handlers so a missing image doesn't abort
    # the session restore midway.
    state["ra_deg"] = float(sess.pointing_ra_deg)
    state["dec_deg"] = float(sess.pointing_dec_deg)
    ra_input.value = f"{sess.pointing_ra_deg:.6f}"
    dec_input.value = f"{sess.pointing_dec_deg:.6f}"
    _sync_pa_widgets(sess.pa_v3_deg)
    combo_label = f"{sess.disperser} / {sess.filter_name}"
    if combo_label in DISPERSER_FILTER_LABELS:
        disperser_filter_select.value = combo_label
    else:
        state["disperser"] = sess.disperser
        state["filter"] = sess.filter_name
    slitlet_select.value = str(sess.slitlet_height)
    state["slitlet_height"] = int(sess.slitlet_height)
    state["open_shutters"] = {
        (sh.q, sh.s, sh.d): sh for sh in sess.open_shutters
    }
    state["highlighted"] = set(sess.highlighted)
    state["history"] = []

    # Now try to load the image. Route by extension so a JPG session goes
    # through the JPG+sidecar loader, not load_fits. Tolerate failures —
    # the user keeps the picks and can re-load the image manually.
    image_note = ""
    img_path = Path(sess.image_path) if sess.image_path else None
    if img_path and img_path.exists():
        ext = img_path.suffix.lower()
        if ext in (".jpg", ".jpeg", ".png"):
            sidecar = _resolve_jpg_sidecar(img_path, sess.wcs_sidecar_path)
            if sidecar is None:
                image_note = (
                    f" Image '{img_path.name}' needs a WCS sidecar FITS — "
                    f"none found alongside it; load the image manually."
                )
            else:
                sidecar_path_input.value = sidecar
                jpg_path_input.value = str(img_path)  # triggers JPG load
        elif ext in (".fits", ".fit", ".fts"):
            fits_path_input.value = str(img_path)  # triggers FITS load
        else:
            image_note = f" Unknown image extension {ext!r}; load manually."
    elif sess.image_path:
        image_note = f" Image not found at {sess.image_path}; load manually."

    if sess.catalog_path and Path(sess.catalog_path).exists() and (
        state["catalog"] is None
        or getattr(state["catalog"], "source_path", None) != sess.catalog_path
    ):
        catalog_path_input.value = sess.catalog_path

    refresh_overlays()
    if state["image"] is None and not image_note:
        image_note = " Load an image to see the overlay."
    _set_status(
        f"Session loaded: {len(state['open_shutters'])} open shutters, "
        f"{len(state['highlighted'])} highlighted.{image_note}",
        "warn" if image_note else "ok", clear_after=14,
    )


session_save_btn.on_click(on_session_save)
session_load_btn.on_click(on_session_load)


# ---------------------------------------------------------------------------
# Example data quick-load
# ---------------------------------------------------------------------------


_EX_DIR = Path(__file__).resolve().parent.parent

def _reset_pointing_inputs() -> None:
    """Clear RA/Dec inputs so the next image-load auto-recenters."""
    ra_input.value = ""
    dec_input.value = ""


def on_example_a370():
    p = _EX_DIR / "example_a370" / "a370_f182m_f200w_f210m.fits"
    if not p.exists():
        _set_status(f"Example missing: {p}", "err")
        return
    # Example buttons are a hard reset: clear pointing so the new image
    # auto-recenters (otherwise loading R0600 after A370 leaves pointing
    # at A370 and the MSA disappears off-screen).
    _reset_pointing_inputs()
    fits_path_input.value = str(p)


def on_example_r0600():
    jpg = _EX_DIR / "example_r0600" / "JWST_F090W_F200W_F444W.jpg"
    wcs = _EX_DIR / "example_r0600" / "wcs.fits"
    if not (jpg.exists() and wcs.exists()):
        _set_status("Example r0600 files missing.", "err")
        return
    _reset_pointing_inputs()
    sidecar_path_input.value = str(wcs)
    jpg_path_input.value = str(jpg)


example_a370_btn.on_click(on_example_a370)
example_r0600_btn.on_click(on_example_r0600)


# ---------------------------------------------------------------------------
# MPT JSON / shutter CSV import (load an existing APT MSA plan)
# ---------------------------------------------------------------------------


# Cache the parsed plan list keyed on file path so the Select callback
# doesn't re-parse on every plan-switch.
_mpt_plans_cache: dict[str, list] = {}


def _parse_and_cache_mpt(path: str) -> list:
    plans = parse_mpt_json(path)
    _mpt_plans_cache[path] = plans
    return plans


def on_mpt_json_path(attr, old, new):
    path = mpt_json_path_input.value.strip()
    if not path:
        mpt_plan_select.options = []
        mpt_plan_select.visible = False
        mpt_load_btn.disabled = True
        return
    if not Path(path).exists():
        _set_status(f"MPT JSON path not found: {path}", "err")
        return
    try:
        plans = _parse_and_cache_mpt(path)
    except ValueError as e:
        _set_status(f"MPT JSON parse failed: {e}", "err")
        return
    labels = [f"{i+1}. {p.name} ({p.n_open_shutters} open)"
              for i, p in enumerate(plans)]
    mpt_plan_select.options = labels
    mpt_plan_select.value = labels[0] if labels else ""
    mpt_plan_select.visible = True
    mpt_load_btn.disabled = not labels
    _set_status(f"Found {len(plans)} plan(s) in MPT JSON. Pick one and click "
                f"'Load plan from JSON'.", "ok")


def on_mpt_load():
    path = mpt_json_path_input.value.strip()
    if not path or path not in _mpt_plans_cache:
        _set_status("Set the MPT JSON path first.", "warn")
        return
    plans = _mpt_plans_cache[path]
    sel = mpt_plan_select.value
    if not sel:
        _set_status("Pick a plan from the dropdown.", "warn")
        return
    idx = mpt_plan_select.options.index(sel)
    plan = plans[idx]
    _apply_plan(plan)


def _apply_plan(plan) -> None:
    """Apply an MPTPlan to the live UI state: pointing, V3 PA, disperser,
    and the unfolded open-shutter set."""
    _push_history()
    if plan.ra_deg is not None and plan.dec_deg is not None:
        ra_input.value = f"{plan.ra_deg:.6f}"
        dec_input.value = f"{plan.dec_deg:.6f}"
        state["ra_deg"] = plan.ra_deg
        state["dec_deg"] = plan.dec_deg
    _sync_pa_widgets(plan.v3_pa_deg)
    if plan.grating and plan.filter_name:
        combo = f"{plan.grating} / {plan.filter_name}"
        if combo in DISPERSER_FILTER_LABELS:
            disperser_filter_select.value = combo
    state["open_shutters"] = {
        (sh.q, sh.s, sh.d): sh for sh in plan.to_open_shutters()
    }
    n_open = len(state["open_shutters"])
    img: LoadedImage | None = state["image"]
    if img is None:
        _set_status(
            f"Loaded plan '{plan.name}': {n_open} open shutters at "
            f"APA={plan.aperture_pa_deg:.2f}°, V3 PA={plan.v3_pa_deg:.2f}°. "
            f"Load an image (Image tab) to see and edit the overlay.",
            "warn", clear_after=20,
        )
        return
    refresh_overlays()
    _set_status(
        f"Loaded plan '{plan.name}': {n_open} open shutters, "
        f"APA={plan.aperture_pa_deg:.2f}°, V3 PA={plan.v3_pa_deg:.2f}°.",
        "ok", clear_after=12,
    )


def on_mpt_csv_load():
    path = mpt_csv_path_input.value.strip()
    if not path:
        _set_status("Set the shutter CSV path first.", "warn")
        return
    if not Path(path).exists():
        _set_status(f"Shutter CSV not found: {path}", "err")
        return
    try:
        opens = parse_shutter_csv(path)
    except ValueError as e:
        _set_status(f"Shutter CSV parse failed: {e}", "err")
        return
    _push_history()
    state["open_shutters"] = {(sh.q, sh.s, sh.d): sh for sh in opens}
    refresh_overlays()
    _set_status(
        f"Loaded shutter CSV: {len(opens)} open shutters (no slitlet "
        f"grouping or target IDs — CSV doesn't carry that info).",
        "ok", clear_after=12,
    )


mpt_json_path_input.on_change("value", on_mpt_json_path)
mpt_load_btn.on_click(on_mpt_load)
mpt_csv_load_btn.on_click(on_mpt_csv_load)


# .aptx handling: either an on-disk file or a program ID. After we have
# the archive, list its embedded MPT plans into apt_plan_select.

_apt_state: dict = {"aptx_path": None}


def on_apt_fetch():
    aptx_path = apt_path_input.value.strip()
    pid = apt_program_input.value.strip()
    if not aptx_path and not pid:
        _set_status("Set an .aptx path or a program ID first.", "warn")
        return
    if aptx_path:
        if not Path(aptx_path).exists():
            _set_status(f".aptx not found: {aptx_path}", "err")
            return
        path = aptx_path
    else:
        _show_loading(f"Downloading APT {pid} from STScI…")
        try:
            path = download_apt_program(pid)
        except ValueError as e:
            _hide_loading()
            _set_status(f"Download failed: {e}", "err")
            return
        _hide_loading()
    try:
        plans = list_mpt_plans_in_aptx(path)
    except Exception as e:  # noqa: BLE001
        _set_status(f"Could not read .aptx: {e}", "err")
        return
    if not plans:
        _set_status(f".aptx contained no MPT-format plans.", "warn")
        return
    _apt_state["aptx_path"] = path
    apt_plan_select.options = plans
    apt_plan_select.value = plans[0]
    apt_plan_select.visible = True
    apt_load_btn.disabled = False
    _set_status(
        f"Found {len(plans)} MSA plan(s) in {Path(path).name}. "
        f"Pick one and click 'Load selected plan'.", "ok", clear_after=15,
    )


def on_apt_load():
    aptx_path = _apt_state.get("aptx_path")
    member = apt_plan_select.value
    if not aptx_path or not member:
        _set_status("Open an .aptx first, then pick a plan.", "warn")
        return
    try:
        plans = parse_mpt_json_in_aptx(aptx_path, member)
    except ValueError as e:
        _set_status(f"Plan parse failed: {e}", "err")
        return
    if not plans:
        _set_status(f"No plans inside {member}.", "warn")
        return
    # MPT JSON inside an .aptx member may have multiple configs (sub-plans).
    # For now, load the first; future enhancement: a second dropdown.
    _apply_plan(plans[0])
    if len(plans) > 1:
        _set_status(
            f"Loaded first of {len(plans)} configs in {member}. "
            f"(Future: select sub-config from a dropdown.)",
            "ok", clear_after=20,
        )


apt_fetch_btn.on_click(on_apt_fetch)
apt_load_btn.on_click(on_apt_load)


# ---------------------------------------------------------------------------
# Native file picker (replaces Bokeh's FileInput upload widgets)
# ---------------------------------------------------------------------------


def _pick_file_dialog(
    title: str,
    filetypes: list[tuple[str, str]] | None = None,
    mode: str = "open",
) -> str | None:
    """Pop a native file dialog and return the chosen path (or None).

    Runs tkinter in a subprocess so it can't interfere with the Bokeh
    server's Tornado event loop, and so any tkinter init weirdness is
    isolated. Blocks the calling Bokeh callback while the dialog is
    open — typical, and the user expects it.

    `mode` is "open", "save", or "directory".
    """
    import json as _json
    import subprocess as _sp
    import sys as _sys
    ft_repr = repr(filetypes or [])
    script = (
        "import json, tkinter as tk\n"
        "from tkinter import filedialog\n"
        "root = tk.Tk()\n"
        "root.withdraw()\n"
        "root.attributes('-topmost', True)\n"
        f"title = {title!r}\n"
        f"ftypes = {ft_repr}\n"
        f"mode = {mode!r}\n"
        "if mode == 'directory':\n"
        "    p = filedialog.askdirectory(title=title, mustexist=True)\n"
        "elif mode == 'save':\n"
        "    p = filedialog.asksaveasfilename(title=title, filetypes=ftypes)\n"
        "else:\n"
        "    p = filedialog.askopenfilename(title=title, filetypes=ftypes)\n"
        "print(json.dumps(p or ''))\n"
        "root.destroy()\n"
    )
    try:
        out = _sp.run(
            [_sys.executable, "-c", script],
            capture_output=True, text=True, timeout=300,
        )
    except _sp.SubprocessError as e:
        _set_status(f"File picker failed: {e}", "err")
        return None
    last = out.stdout.strip().splitlines()[-1] if out.stdout.strip() else ""
    try:
        result = _json.loads(last) if last else ""
    except Exception:  # noqa: BLE001
        result = ""
    return result or None


_FITS_TYPES = [("FITS", "*.fits *.fit"), ("All files", "*")]
_JPG_TYPES = [("Image", "*.jpg *.jpeg *.png"), ("All files", "*")]
_CATALOG_TYPES = [("Catalog", "*.csv *.cat *.txt *.fits *.ascii"),
                  ("All files", "*")]
_JSON_TYPES = [("JSON", "*.json"), ("All files", "*")]
_CSV_TYPES = [("CSV", "*.csv"), ("All files", "*")]
_APTX_TYPES = [("APT archive", "*.aptx *.zip"), ("All files", "*")]


def _bind_browse(button, target_input, title, filetypes, mode="open"):
    """Wire one Browse… button to write its picked path into a TextInput."""
    def _cb():
        p = _pick_file_dialog(title, filetypes, mode)
        if p:
            target_input.value = p
    button.on_click(_cb)


_bind_browse(fits_browse_btn, fits_path_input, "Select FITS image", _FITS_TYPES)
_bind_browse(jpg_browse_btn, jpg_path_input, "Select JPG/PNG image", _JPG_TYPES)
_bind_browse(sidecar_browse_btn, sidecar_path_input, "Select sidecar FITS (WCS)", _FITS_TYPES)
_bind_browse(catalog_browse_btn, catalog_path_input, "Select catalog file", _CATALOG_TYPES)
_bind_browse(mpt_json_browse_btn, mpt_json_path_input, "Select MPT JSON plan", _JSON_TYPES)
_bind_browse(mpt_csv_browse_btn, mpt_csv_path_input, "Select shutter mask CSV", _CSV_TYPES)
_bind_browse(apt_path_browse_btn, apt_path_input, "Select APT (.aptx)", _APTX_TYPES)
_bind_browse(session_save_browse_btn, session_save_path_input,
             "Save session as…", _JSON_TYPES, mode="save")
_bind_browse(session_load_browse_btn, session_load_path_input,
             "Load session JSON", _JSON_TYPES)
_bind_browse(export_dir_browse_btn, export_dir_input,
             "Pick export directory", None, mode="directory")
snap_box.on_change("active", on_snap)


# ── Overlay-appearance picker callbacks ──────────────────────────────────
# Each layer entry describes which glyph property the alpha and stroke
# sliders write to, the slider ranges/steps, and the layer's default
# values. When the user changes the dropdown, the sliders are
# reconfigured to that layer's profile.
#
# `stroke_attr` is normally "line_width" but for catalog markers we use
# "size" (marker diameter in screen px) since stroke is data-driven there.
_OVERLAY_LAYER_CONFIG = {
    "Operable shutters":     {
        "glyph": bg_shutters_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
    },
    "Overlapping shutters":  {
        "glyph": spec_overlap_glyph,
        "alpha_attr": "fill_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 0.8, 0.02), "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
        # When the user moves stroke off 0, also reveal line_alpha so the
        # outline actually shows.
        "stroke_extra": lambda v: setattr(
            spec_overlap_glyph.glyph, "line_alpha", 0.6 if v > 0 else 0.0
        ),
    },
    "Picked shutters":       {
        "glyph": open_shutters_glyph,
        "alpha_attr": "fill_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 4.0, 0.1),
        "stroke_label": "Stroke (px)",
    },
    "Stuck open":            {
        "glyph": stuck_open_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 5.0, 0.1),
        "stroke_label": "Stroke (px)",
    },
    "Catalog sources":       {
        "glyph": target_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "size",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (4.0, 30.0, 1.0),
        "stroke_label": "Marker size (px)",
    },
}

# Guard flag — programmatic slider updates from _on_overlay_layer
# shouldn't cascade back through _on_overlay_alpha / _on_overlay_stroke
# and stomp on the glyph properties we just read out.
_overlay_syncing = {"flag": False}


def _on_overlay_layer(attr, old, new):
    cfg = _OVERLAY_LAYER_CONFIG.get(new)
    if cfg is None:
        return
    _overlay_syncing["flag"] = True
    try:
        a_lo, a_hi, a_step = cfg["alpha_range"]
        s_lo, s_hi, s_step = cfg["stroke_range"]
        overlay_alpha_slider.update(start=a_lo, end=a_hi, step=a_step)
        overlay_stroke_slider.update(
            start=s_lo, end=s_hi, step=s_step, title=cfg["stroke_label"],
        )
        # Read current values from the glyph and reflect them.
        cur_alpha = getattr(cfg["glyph"].glyph, cfg["alpha_attr"])
        cur_stroke = getattr(cfg["glyph"].glyph, cfg["stroke_attr"])
        # Bokeh string-field references aren't numbers — fall back to a
        # sensible default if the property is data-driven.
        try:
            overlay_alpha_slider.value = max(a_lo, min(a_hi, float(cur_alpha)))
        except (TypeError, ValueError):
            overlay_alpha_slider.value = a_lo
        try:
            overlay_stroke_slider.value = max(s_lo, min(s_hi, float(cur_stroke)))
        except (TypeError, ValueError):
            overlay_stroke_slider.value = s_lo
    finally:
        _overlay_syncing["flag"] = False


def _on_overlay_alpha(attr, old, new):
    if _overlay_syncing["flag"]:
        return
    cfg = _OVERLAY_LAYER_CONFIG.get(overlay_layer_select.value)
    if cfg is None:
        return
    setattr(cfg["glyph"].glyph, cfg["alpha_attr"], float(new))


def _on_overlay_stroke(attr, old, new):
    if _overlay_syncing["flag"]:
        return
    cfg = _OVERLAY_LAYER_CONFIG.get(overlay_layer_select.value)
    if cfg is None:
        return
    setattr(cfg["glyph"].glyph, cfg["stroke_attr"], float(new))
    extra = cfg.get("stroke_extra")
    if extra is not None:
        extra(float(new))


overlay_layer_select.on_change("value", _on_overlay_layer)
overlay_alpha_slider.on_change("value", _on_overlay_alpha)
overlay_stroke_slider.on_change("value", _on_overlay_stroke)


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

# SIDEBAR_W / HELPPANEL_W are defined above (near the figure init) and
# referenced by the figure-sizing comment. Don't redeclare here.

# Compact sidebar organised as tabs around the natural workflow:
# Load (image + catalog) → Aim (pointing + PA + visibility) →
# Pick (instrument, layers, slitlet, filters, undo/clear) →
# Save (session save/load, APT export).

image_tab = TabPanel(title="Image", child=column(
    Div(text="<b>Image</b> — try an example:"),
    row(example_a370_btn, example_r0600_btn),
    Div(text="<small><b>or</b> a local FITS:</small>"),
    row(fits_path_input, fits_browse_btn),
    Div(text="<small><b>or</b> JPG + sidecar FITS:</small>"),
    row(sidecar_path_input, sidecar_browse_btn),
    row(jpg_path_input, jpg_browse_btn),
    Div(text="<b>Catalog</b> <small>(CSV / ASCII / FITS with ID, RA, DEC)</small>"),
    row(catalog_path_input, catalog_browse_btn),
    catalog_priority_input,
    catalog_mag_input,
    width=SIDEBAR_W - 20,
))

aim_tab = TabPanel(title="Aim", child=column(
    Div(text="<b>Pointing center</b>"),
    ra_input, dec_input,
    Div(text="<b>Rotation</b>"),
    v3pa_slider, v3pa_input, apa_input, pa_help_div,
    Div(text="<b>Visibility window</b>"),
    visibility_date_input, visibility_btn, visibility_div,
    width=SIDEBAR_W - 20,
))

pick_tab = TabPanel(title="Pick", child=column(
    Div(text="<b>Disperser / Filter</b>"),
    disperser_filter_select,
    Div(text="<b>Layers</b>"),
    layers_box,
    Div(text="<b>Slitlet</b>"),
    slitlet_select,
    snap_box,
    Div(text="<b>Overlay appearance</b>"),
    overlay_layer_select,
    overlay_alpha_slider,
    overlay_stroke_slider,
    Div(text="<b>Actions</b>"),
    row(undo_btn, clear_btn),
    width=SIDEBAR_W - 20,
))

# MPT tab — everything to do with APT / MPT plans: import (from JSON,
# shutter CSV, or .aptx archive / program ID), the session save/load
# round-trip for collaboration, and the eMPT export bundle.
mpt_tab = TabPanel(title="MPT", child=column(
    Div(text="<b>Import a plan</b> — from a single MPT JSON:"),
    row(mpt_json_path_input, mpt_json_browse_btn),
    mpt_plan_select,
    mpt_load_btn,
    Div(text="<small><i>or</i> a shutter CSV (open mask only):</small>"),
    row(mpt_csv_path_input, mpt_csv_browse_btn),
    mpt_csv_load_btn,
    Div(text="<small><i>or</i> straight from an .aptx file or program ID:</small>"),
    row(apt_path_input, apt_path_browse_btn),
    apt_program_input,
    apt_fetch_btn,
    apt_plan_select,
    apt_load_btn,
    Div(text="<b>Save / share session</b>"),
    row(session_save_path_input, session_save_browse_btn),
    session_save_btn,
    row(session_load_path_input, session_load_browse_btn),
    session_load_btn,
    Div(text=f"<b>Export to APT</b> (eMPT bundle + {MPT_PLAN_FILENAME})"),
    row(export_dir_input, export_dir_browse_btn),
    export_btn,
    width=SIDEBAR_W - 20,
))

sidebar_tabs = Tabs(
    tabs=[image_tab, aim_tab, pick_tab, mpt_tab],
    width=SIDEBAR_W,
)

# Sidebar (no stats — stats now live above the figure as a wide bar).
# The MPT tab content (import + save/load + export) is taller than the
# viewport on most laptops, so the sidebar scrolls internally on overflow.
sidebar = column(
    loading_banner,
    sidebar_tabs,
    status,
    width=SIDEBAR_W,
    height_policy="max",
    styles={"overflow-y": "auto", "max-height": "100vh"},
)

# Figure column: status bar on top (stretches to canvas width) and the
# stretch-both canvas below. The whole column fills the centre cell of
# the root row, which itself stretches to fit the browser window.
#
# Layout shape on a 1440-wide laptop:
#   [ sidebar 340 ][ status bar | canvas (stretches) ][ help 340 ]
#                                      ~760 px wide
# On a 1920-wide monitor the canvas just gets bigger (~1240 px).
# `match_aspect=True` on the figure keeps data pixels 1:1 regardless.
figure_column = column(
    stats_div,
    fig,
    sizing_mode="stretch_both",
)

curdoc().add_root(row(
    sidebar, figure_column, help_panel,
    sizing_mode="stretch_both",
))
curdoc().title = "vMPT — visual MSA Planning Tool"


# Rotate the tip card every 15 s. _render_tip carries its own @keyframes
# fade-in block; setting the Bokeh Div's `text` swaps the DOM so the
# animation restarts cleanly without us needing to hand-splice strings.
def _advance_tip() -> None:
    _tip_state["idx"] = (_tip_state["idx"] + 1) % len(_TIPS)
    tip_div.text = _render_tip(_tip_state["idx"])


curdoc().add_periodic_callback(_advance_tip, 15_000)


# ──────────────────────────────────────────────────────────────────────────
# GitHub version-check (non-blocking)
# ──────────────────────────────────────────────────────────────────────────
# On startup, in a background thread, ask GitHub for the latest commit
# on `main`. If the local checkout's HEAD differs, show a dismissible
# notification overlay. Failures (no network, no git, rate-limit) are
# silent — version-checking is a courtesy, not a hard dependency.
import json as _json
import subprocess as _subprocess
import threading as _threading
import urllib.error as _urllib_error
import urllib.request as _urllib_request

_GITHUB_OWNER = "fengwusun"
_GITHUB_REPO = "vMPT"
_GITHUB_WEB_URL = f"https://github.com/{_GITHUB_OWNER}/{_GITHUB_REPO}"


def _github_compare_url(local_sha: str) -> str:
    """GitHub's compare endpoint tells us whether the user's local HEAD
    is ahead/behind/diverged from main, plus the list of new commits.
    We use the API JSON variant (not the HTML page) for status detection."""
    return (
        f"https://api.github.com/repos/{_GITHUB_OWNER}/"
        f"{_GITHUB_REPO}/compare/{local_sha}...main"
    )


def _local_git_head() -> str | None:
    """Return the SHA of the local checkout's HEAD, or None if not in a
    git working tree."""
    try:
        result = _subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(Path(__file__).resolve().parent.parent),
            capture_output=True, text=True, timeout=2, check=False,
        )
    except (FileNotFoundError, _subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    sha = result.stdout.strip()
    return sha if len(sha) == 40 else None


# The notification widgets — built up-front, hidden by default; the
# background thread reveals them via a doc tick callback when an update
# is detected. CSS `position: fixed` makes them float above the layout
# so they don't displace the canvas.
update_div = Div(
    text="",
    width=360,
    styles={
        "background": "linear-gradient(135deg, #fff3cd 0%, #ffe9a0 100%)",
        "border": "1px solid #ddb84c",
        "border-radius": "8px",
        "padding": "12px 14px 6px 14px",
        "box-shadow": "0 4px 16px rgba(0,0,0,0.15)",
        "color": "#5a4500",
        "font-size": "12.5px",
        "line-height": "1.45",
    },
)
update_dismiss_btn = Button(
    label="Dismiss", button_type="default", width=80, height=28,
)
update_box = column(
    update_div, update_dismiss_btn,
    visible=False,
    styles={
        "position": "fixed",
        "top": "70px",
        "right": "20px",
        "z-index": "99999",
        "max-width": "380px",
    },
)


def _on_update_dismiss() -> None:
    update_box.visible = False


update_dismiss_btn.on_click(_on_update_dismiss)
curdoc().add_root(update_box)


def _on_update_available(
    local_sha: str, remote_sha: str, commit_msg: str, n_commits: int,
) -> None:
    """Runs on the Bokeh document thread. Reveals the update notification."""
    short_local = local_sha[:7]
    short_remote = remote_sha[:7] if len(remote_sha) >= 7 else remote_sha
    safe_msg = (commit_msg or "").splitlines()[0][:140]
    safe_msg = safe_msg.replace("<", "&lt;").replace(">", "&gt;")
    commits_label = (
        f"{n_commits} new commit{'s' if n_commits != 1 else ''} on GitHub"
        if n_commits else "Newer commits on GitHub"
    )
    update_div.text = (
        '<div style="display:flex; align-items:flex-start; gap:10px;">'
        '  <div style="font-size:22px; line-height:1; flex-shrink:0;">📦</div>'
        '  <div style="min-width:0;">'
        '    <div style="font-weight:700; color:#7a5d00; font-size:13px; '
        '                margin-bottom:4px; letter-spacing:0.2px;">'
        '      vMPT update available'
        '    </div>'
        '    <div style="font-size:11.5px; color:#6f5b1f; margin-bottom:6px;">'
        f'      {commits_label}: <code>{short_local}</code> → '
        f'      <code>{short_remote}</code>'
        '    </div>'
        f'    <div style="font-size:11.5px; color:#6f5b1f; margin-bottom:8px; '
        '                 font-style:italic;">'
        f'      Latest: {safe_msg}'
        '    </div>'
        f'    <a href="{_GITHUB_WEB_URL}/compare/{short_local}...main" '
        '       target="_blank" rel="noopener" '
        '       style="color:#7a5d00; font-weight:600; text-decoration:none; '
        '              font-size:11.5px; border-bottom:1px solid #7a5d0033;">'
        '      Review changes on GitHub →'
        '    </a>'
        '    <div style="font-size:10.5px; color:#8a7530; margin-top:6px;">'
        '      Run <code>git pull</code> and restart the server to update.'
        '    </div>'
        '  </div>'
        '</div>'
    )
    update_box.visible = True


def _check_for_updates_blocking(session_doc) -> None:
    """Run in a background thread. `session_doc` is the user-session doc
    captured BEFORE the thread started — curdoc() is thread-local, so
    re-querying it from this thread returns a different (irrelevant)
    document and the next-tick callback never fires on the user's UI."""
    local = _local_git_head()
    if not local:
        return  # not a git checkout — nothing to compare against
    try:
        req = _urllib_request.Request(
            _github_compare_url(local),
            headers={
                "Accept": "application/vnd.github+json",
                "User-Agent": "vMPT-version-check",
            },
        )
        with _urllib_request.urlopen(req, timeout=6) as resp:
            data = _json.loads(resp.read().decode("utf-8"))
    except (_urllib_error.URLError, _urllib_error.HTTPError, TimeoutError,
            OSError, _json.JSONDecodeError, UnicodeDecodeError):
        return  # offline / rate-limited / unparseable — silently skip
    # GitHub's compare/BASE...HEAD reports `status` from HEAD's
    # perspective relative to BASE. We call compare(<local_sha>...main),
    # so BASE=local, HEAD=main. The meaningful statuses for us:
    #   "ahead"     → main has commits the local checkout doesn't (USER NEEDS UPDATE)
    #   "behind"    → local has commits not on main (user is ahead of remote)
    #   "identical" → up to date
    #   "diverged"  → both sides have unique commits (fork / WIP branch)
    # We only nag on "ahead". Leave "behind"/"diverged"/"identical" alone.
    status = (data.get("status") or "").strip().lower()
    if status != "ahead":
        return
    remote_head = ((data.get("commits") or [{}])[-1].get("sha") or "").strip()
    head_commit = (data.get("commits") or [{}])[-1].get("commit") or {}
    commit_msg = (head_commit.get("message") or "").strip()
    n_commits = len(data.get("commits") or [])
    try:
        session_doc.add_next_tick_callback(
            lambda: _on_update_available(
                local, remote_head or "main", commit_msg, n_commits,
            )
        )
    except Exception:  # noqa: BLE001
        pass


# Only run the version-check when we're being served by Bokeh — not
# during pytest imports or REPL imports. Capture the session doc here
# (on the main session thread); the worker thread can't query curdoc()
# itself because curdoc() is thread-local.
_session_doc = curdoc()
if _session_doc.session_context is not None:
    _threading.Thread(
        target=_check_for_updates_blocking,
        args=(_session_doc,),
        daemon=True,
    ).start()
