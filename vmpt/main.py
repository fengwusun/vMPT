"""vMPT — visual MSA Planning Tool. Bokeh server entry point.

Run:  bokeh serve vmpt/ --show   (or `vmpt` after `pip install jwst-vmpt`)
"""
from __future__ import annotations

import base64
import os
import re
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

from bokeh.events import MouseLeave, MouseMove, RangesUpdate, Tap
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.events import DocumentReady
from bokeh.models import (
    Button,
    CheckboxGroup,
    ColumnDataSource,
    CustomJS,
    CustomJSTickFormatter,
    DataTable,
    Div,
    GlobalInlineStyleSheet,
    HoverTool,
    InlineStyleSheet,
    HTMLTemplateFormatter,
    MultiChoice,
    NumberEditor,
    NumberFormatter,
    Range1d,
    RadioGroup,
    Select,
    Slider,
    Spinner,
    StringEditor,
    TableColumn,
    TabPanel,
    Tabs,
    TextInput,
    Toggle,
    WheelZoomTool,
)
from bokeh.plotting import figure

from vmpt.catalog import (
    Catalog,
    catalog_in_view,
    evaluate_catalog_condition,
    load_catalog,
)
from vmpt.coords import (
    MSA_V2_REF,
    MSA_V3_REF,
    V3_IDL_Y_ANGLE,
    fixed_slit_corners_v2v3,
    rot_matrix,
    shutter_corners_v2v3,
    v2v3_to_radec,
)
from vmpt.empt_io import (
    OpenShutter,
    Pointing,
    write_mpt_catalog,
    write_observed_targets_cat,
    write_pointing_summary_txt,
    write_shutter_mask_csv,
)
from vmpt.image_io import LoadedImage, load_fits, load_jpg_with_sidecar, stretch_for_display
from vmpt.msa import ensure_current_operability, load_msa_grid, load_operability
from vmpt.mpt_io import (
    MPTPlan,
    download_apt_program,
    list_mpt_plans_in_aptx,
    parse_mpt_json,
    parse_mpt_json_in_aptx,
    parse_shutter_csv,
)
from vmpt.session_io import (
    EMPT_OBSERVED_FILENAME,
    EMPT_POINTING_FILENAME,
    EMPT_SHUTTER_MASK_FILENAME,
    MPT_CATALOG_FILENAME,
    MPT_PLAN_FILENAME,
    Session,
    WORKSPACE_FILENAME,
    apt_catalog_basename,
    apt_plan_basename,
    bundle_readme_text,
    export_session_json,
    import_session_json,
)
from vmpt.wavelengths import (
    FILTER_BLUE_CUTOFF,
    GRATING_RANGES,
    cutoffs,
    grating_adjacency_min_colsep,
    v2_overlap_distance,
)


def _resubscribe_late_event_handlers(*models) -> None:
    """Re-enrol models in their document's event dispatcher.

    Bokeh only subscribes a model to ``Document.callbacks._subscribed_models``
    when the model is *attached* to the document (``_attach_document`` runs
    ``_update_event_callbacks``). A handler wired with ``on_event`` /
    ``on_click`` *after* the model's layout root was added to ``curdoc()`` is
    recorded on the model's ``_event_callbacks`` but never enrolled with the
    document, so the Bokeh server silently drops its events (e.g. a button
    that "does nothing" when clicked). Call this right after any such late
    wiring. Idempotent — ``subscribe`` is a ``set.add`` — and a no-op when
    the model isn't attached yet (its eventual attach will subscribe it).
    """
    for model in models:
        try:
            if getattr(model, "document", None) is not None:
                model._update_event_callbacks()
        except Exception:  # noqa: BLE001 — never let a wiring helper crash startup
            pass


# ---------------------------------------------------------------------------
# Singletons
# ---------------------------------------------------------------------------

V2_MSA, V3_MSA = load_msa_grid()           # (4, 171, 365)
# Each time vMPT starts, make sure the local CRDS cache has the current MOS
# operability reference (best-effort; offline-safe, ~once per process). Then
# load it — so the failed-/stuck-shutter map matches APT/MPT.
ensure_current_operability()
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

def _new_config(name: str) -> dict:
    """Create one empty MPT configuration.

    Each config owns its picks (`open_shutters`), visual highlights, undo
    history, AND its own saved pointing (`ra_deg`/`dec_deg`/`pa_v3`). A
    pointing of None means "inherit the live widgets when this config is
    first activated" (used for the freshly-created config). While a config
    is the *active* one, its live pointing is the global `state['ra_deg']`
    etc. + the widgets; the saved copy is refreshed when you switch away.
    """
    return {
        "name": name,
        "open_shutters": {},   # (q,s,d) -> OpenShutter
        "highlighted": set(),  # set of (q,s,d) tuples
        "history": [],         # undo stack of open_shutters snapshots
        "ra_deg": None,
        "dec_deg": None,
        "pa_v3": None,
    }


# The first config exists from the start; `state['open_shutters']`,
# `['highlighted']` and `['history']` are LIVE ALIASES of the active
# config's dicts so the ~40 existing in-place readers/mutators keep working
# unchanged. Only WHOLESALE reassignments (`state['open_shutters'] = {...}`)
# must go through `_set_open_shutters()` so the alias and the config slot
# stay in sync.
_config0 = _new_config("Config 1")

state: dict = {
    "image": None,
    "catalog": None,         # merged-active cache; rebuilt from "catalogs"
    "catalogs": [],          # list of {"name", "catalog", "enabled"} entries
                             # — the source of truth. `state["catalog"]` is
                             # recomputed on every add/remove/toggle.
    "tmp_sidecar_path": None,
    # ─── Multi-config (v1.4.0) ───────────────────────────────────────────
    "configs": [_config0],   # list of config dicts (see _new_config)
    "active_config": 0,      # index into configs the user is working on
    "n_configs": 1,          # how many configs are "live" in the plan (1..2)
    # Legacy keys — LIVE ALIASES of configs[active_config][...]:
    "open_shutters": _config0["open_shutters"],  # (q,s,d) -> OpenShutter
    "highlighted": _config0["highlighted"],      # set of (q,s,d) tuples
    "history": _config0["history"],              # undo snapshots (per-config)
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


def _active_config() -> dict:
    """The config dict the user is currently working on."""
    return state["configs"][state["active_config"]]


def _set_open_shutters(new: dict) -> None:
    """Replace the active config's open_shutters wholesale, keeping the
    legacy `state['open_shutters']` alias pointed at the same object."""
    _active_config()["open_shutters"] = new
    state["open_shutters"] = new


def _set_history(new: list) -> None:
    """Replace the active config's undo history, keeping the alias."""
    _active_config()["history"] = new
    state["history"] = new


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
SIDEBAR_W = 340      # left tab panel (Input / Pointing / Setting / MPT)
HELPPANEL_W = 340    # right help panel (Quick guide + rotating tip)

# Per-config accent colours — shared by the top-bar CONFIG chip, the MSA
# quadrant outline (active solid + idle dashed), and the MPT-viewer Cfg
# column so the "which config" colour language is identical everywhere.
# Five distinct, high-contrast accent hues — one per simultaneous
# config. Order is the cycling order of the CONFIG chip:
#   C1 blue · C2 magenta · C3 green · C4 orange · C5 violet.
_CONFIG_COLORS = ["#1f6fc0", "#b5179e", "#2a9d3a", "#e8590c", "#7048e8"]

# The maximum number of simultaneous configs equals the number of
# distinct accent colours, so every live config is always
# colour-distinguishable on the canvas, chip, and MPT viewer.
_MAX_CONFIGS = len(_CONFIG_COLORS)


def _config_color(idx: int) -> str:
    """Accent colour for config ``idx`` (0-based). Rotates through the
    palette (mod) so any index maps to a colour — though in practice
    ``idx`` never exceeds ``_MAX_CONFIGS - 1``."""
    idx = max(0, int(idx))
    return _CONFIG_COLORS[idx % len(_CONFIG_COLORS)]


# Catalogs with more sources than this trigger the full-page loading
# spinner when toggled on/off (re-rendering all their circles takes a beat).
_LARGE_CATALOG_N = 1000

# Inline styles for the header bar at the top of every modal dialog.
# Applied via `styles=` on the header `row()` so they survive Bokeh's
# per-model shadow root (document-level CSS doesn't reach inside
# shadow DOMs — that's why GlobalInlineStyleSheet doesn't work here).
# The `vmpt-modal-header` css_class is still needed for the drag JS
# (which queries by class), but the visual style is set here.
_MODAL_HEADER_STYLES = {
    "cursor": "move",
    "user-select": "none",
    "background": "linear-gradient(180deg, #eef4fc 0%, #d9e6f7 100%)",
    "border-bottom": "1px solid #c2d2e6",
    "border-radius": "6px 6px 0 0",
    "padding": "9px 14px 9px 16px",
    # Extend the header to the modal card's edges by negating its
    # internal padding (cards use 16px 18px).
    "margin": "-16px -18px 12px -18px",
    "display": "flex",
    "align-items": "center",
    "justify-content": "space-between",
    "min-height": "36px",
}

# Color palette for catalog markers when multiple catalogs are loaded.
# Cycled by load order (entry 0 → palette[0], entry 1 → palette[1], …).
# Picked to (a) read clearly on a dark astronomical image, (b) avoid
# the colors already in use by other overlay layers — red (open
# shutters), green (matched targets), cyan (highlighted), gold
# (fixed slits), orange (spec-overlap), lime (pointing handle).
CATALOG_COLOR_PALETTE = (
    "#ffd200",   # 1 — bright yellow (default for single-catalog use)
    "#ff66cc",   # 2 — magenta
    "#9aff8b",   # 3 — pale spring-green (well away from match green)
    "#ff9b3d",   # 4 — coral-orange
    "#b39bff",   # 5 — lavender
    "#7ad9ff",   # 6 — pale sky-blue (lighter than highlighted cyan)
    "#ffffff",   # 7 — white
    "#ff5e5e",   # 8 — salmon (pinker than open-shutter red)
)


# ---------------------------------------------------------------------------
# Bokeh widgets / glyphs
# ---------------------------------------------------------------------------

status = Div(
    text='<div style="color:#888">Ready.</div>',
    # The status bar lives OUTSIDE the scrollable sidebar column —
    # position:fixed pins it to the bottom-left of the viewport so a
    # long message never bleeds into / renders on top of whatever
    # tab content the user has scrolled to. The sidebar gets
    # padding-bottom = status height so its bottom-most widget is
    # never covered. Width matches SIDEBAR_W exactly so the bar
    # spans only under the sidebar.
    width=SIDEBAR_W,
    height=42,
    styles={
        "position": "fixed",
        "bottom": "0",
        "left": "0",
        "width": f"{SIDEBAR_W}px",
        "z-index": "100",
        "padding": "4px 8px",
        "font-size": "11.5px",
        "line-height": "1.35",
        "border-top": "1px solid #e0e6f0",
        "background": "#f7f9fc",
        "box-sizing": "border-box",
        "overflow": "hidden",
    },
)

# Path-based inputs are the primary way to load — no WebSocket size limit, no temp files.
fits_path_input = TextInput(title="FITS path (local)", value="", placeholder="/path/to/image.fits")
jpg_path_input = TextInput(title="JPG path (local)", value="", placeholder="/path/to/image.jpg")
sidecar_path_input = TextInput(title="Sidecar FITS path (WCS for JPG)", value="", placeholder="/path/to/wcs.fits")
catalog_path_input = TextInput(title="Catalog path (local)", value="", placeholder="/path/to/catalog.csv")
catalog_add_btn = Button(label="Add", button_type="primary", width=70)
catalog_edit_btn = Button(label="Edit catalog…", button_type="default",
                          width=SIDEBAR_W - 20)
# Dynamic list of loaded catalogs — populated by `_render_catalog_list()`
# whenever state["catalogs"] changes. Each row in this column is a
# row(CheckboxGroup, Button) pair: checkbox toggles the catalog
# on/off (recomputes the merged target overlay), button × removes it.
catalog_list_column = column(width=SIDEBAR_W - 20)

# Upload widgets work for small files but Bokeh's default WebSocket limit (~20 MB) will
# silently truncate larger ones. Start the server with --websocket-max-message-size if
# you want to use these for big files.
# "Browse…" buttons paired with each path TextInput. Clicking one opens
# a native file picker (via a tkinter subprocess) and writes the chosen
# path into the text input — which triggers the existing on_<path>_path
# callback. No upload, no WebSocket size limit.
fits_browse_btn = Button(label="Browse…", button_type="primary", width=120)
jpg_browse_btn = Button(label="Browse…", button_type="primary", width=120)
sidecar_browse_btn = Button(label="Browse…", button_type="primary", width=120)
catalog_browse_btn = Button(label="Browse…", button_type="primary", width=120)
mpt_json_browse_btn = Button(label="Browse…", button_type="primary", width=120)
mpt_csv_browse_btn = Button(label="Browse…", button_type="primary", width=120)
apt_path_browse_btn = Button(label="Browse…", button_type="primary", width=120)
session_save_browse_btn = Button(label="Browse…", button_type="primary", width=120)
session_load_browse_btn = Button(label="Browse…", button_type="primary", width=120)
export_dir_browse_btn = Button(label="Browse…", button_type="primary", width=120)


# Catalog filters — hide-able. Numeric thresholds; leave blank/empty to skip.
catalog_priority_input = TextInput(
    title="Show priority class ≤ (blank = all)", value="", placeholder="e.g. 3",
)
catalog_mag_input = TextInput(
    title="Show mag ≤ (blank = all)", value="", placeholder="e.g. 28",
)

# Half-width so RA + Dec sit side-by-side in the Pointing tab.
_HALF_W = (SIDEBAR_W - 30) // 2
ra_input = TextInput(title="Pointing RA (deg)", value="", width=_HALF_W)
dec_input = TextInput(title="Pointing Dec (deg)", value="", width=_HALF_W)

# V3 PA = position angle of the JWST V3 axis on sky. This is what drives the
# V2/V3 -> RA/Dec math. APT/MPT's "NIRSpec PA" is the *aperture* PA (APA),
# which differs by the V3IdlYAngle of NRS_FULL_MSA (~138.57 deg). We show
# both, synchronized.
v3pa_slider = Slider(title="V3 PA (deg)", start=0.0, end=360.0, step=0.1, value=0.0)
# V3 PA exact + APA share a row in the Pointing tab.
v3pa_input = TextInput(title="V3 PA (deg, exact)", value="0.0", width=_HALF_W)
apa_input = TextInput(
    title=f"NIRSpec APA — V3PA + {V3_IDL_Y_ANGLE:.2f}°",
    value=f"{V3_IDL_Y_ANGLE % 360.0:.2f}",
    width=_HALF_W,
)
pa_help_div = Div(text=(
    f"<small style='color:#7a8699'>APA = V3 PA + {V3_IDL_Y_ANGLE:.2f}° "
    "(NRS_FULL_MSA) · "
    "<a href='https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/"
    "jwst-position-angles-ranges-and-offsets' target='_blank'>JDox</a></small>"
), width=SIDEBAR_W - 20)

# Visibility window query (jwst_gtvt) — date input + button on one row.
visibility_date_input = TextInput(
    title="Visibility date (YYYY-MM-DD)", value="",
    placeholder="blank = today", width=_HALF_W,
)
visibility_btn = Button(label="Compute allowed V3 PA", button_type="primary",
                       width=_HALF_W, height=42)
visibility_div = Div(text="<small>Allowed V3 PA windows appear here.</small>",
                     width=SIDEBAR_W - 20)

# ── Pointing-optimizer widgets ───────────────────────────────────────────
# Derived from hMPT (Eisenstein, McCarty, Wu; CfA/Harvard) — see
# `app/optimizer.py` for the algorithm. UI exposes a small set of
# essentials inline + an "Advanced" foldout for grid resolution, the
# objective choice, and the PSF σ. After the user clicks Run, a column
# of result rows appears — clicking any row applies that pointing
# (RA, Dec, V3 PA) but does not auto-place picks (per user pref).

opt_dra_input = TextInput(title="ΔRA (arcsec)", value="30", width=_HALF_W)
opt_ddec_input = TextInput(title="ΔDec (arcsec)", value="30", width=_HALF_W)
opt_dpa_input = TextInput(title="ΔPA (deg)", value="30", width=_HALF_W)
opt_n_top_input = TextInput(title="Refine top N", value="5", width=_HALF_W)
opt_method_select = Select(
    title="Method",
    # (value, label) tuples — keep the internal value short (Python
    # code checks for "Democracy" / "Meritocracy" / "Hierarchy")
    # while the user sees a one-line clarifier.
    options=[
        ("Democracy",   "Democracy — most targets"),
        ("Meritocracy", "Meritocracy — highest sum of weights"),
        ("Hierarchy",   "Hierarchy — most top-priority targets (eMPT-style)"),
    ],
    value="Democracy",
    width=SIDEBAR_W - 20,
)
opt_method_help_div = Div(
    # The full method-comparison blurb. Hidden by default — the
    # dropdown's option labels ("Democracy — most targets" etc.) are
    # self-describing, and the inline blurb made the Pointing tab
    # taller than most laptops' viewport. Toggled by clicking the
    # ⓘ helper next to the dropdown.
    text=(
        "<small style='color:#5a6b85; line-height:1.4'>"
        "<b>Democracy</b>: maximises raw count, ignores priority/weight.<br>"
        "<b>Meritocracy</b>: maximises Σ weight of placed sources "
        "(requires <code>weight</code> column).<br>"
        "<b>Hierarchy</b>: strict tier order — best for top priority "
        "first, ties broken by next tier "
        "(requires <code>priority</code> column).</small>"
    ),
    width=SIDEBAR_W - 20,
    visible=False,
)
# Inline ⓘ that toggles the blurb on demand. Bokeh has no native
# tooltip widget, so we fake it with a click-to-expand Div + small
# button. Keeps the default Pointing-tab layout compact while still
# letting curious users reveal the method comparison without leaving
# the page.
opt_method_help_toggle = Button(
    label="ⓘ What do these mean?",
    button_type="default",
    width=210,
    height=24,
    css_classes=["vmpt-help-toggle"],
    # The document-level <style> can't reach a Button's <button> (it lives
    # in the widget's shadow root), so inject a stylesheet into that root
    # to render this as a quiet inline link rather than a boxed button.
    stylesheets=[InlineStyleSheet(css="""
      .bk-btn, button {
        background: transparent !important; border: 0 !important;
        box-shadow: none !important; color: #4a7ab8 !important;
        text-align: left !important; padding: 2px 0 !important;
        font-size: 12px !important; font-weight: 400 !important;
        cursor: pointer;
      }
      .bk-btn:hover, button:hover {
        color: #1a3b66 !important; text-decoration: underline;
      }
    """)],
)


def _toggle_method_help() -> None:
    opt_method_help_div.visible = not opt_method_help_div.visible
    opt_method_help_toggle.label = (
        "ⓘ Hide method help" if opt_method_help_div.visible
        else "ⓘ What do these mean?"
    )


opt_method_help_toggle.on_click(_toggle_method_help)
opt_centration_select = Select(
    title="Source centering",
    options=["UNCONSTRAINED", "ENTIRE_OPEN", "MIDPOINT",
             "CONSTRAINED", "TIGHTLY_CONSTRAINED"],
    value="UNCONSTRAINED",
    width=SIDEBAR_W - 20,
)
# Per-target centration override hint (v1.3.1+). Lives directly under
# `opt_centration_select` in the optimizer modal; populated by
# `_refresh_centration_override_hint()` whenever the catalog changes
# (catalog load, editor Apply, session reload). Empty unless ≥1 row
# carries a non-blank `centration` field.
opt_centration_override_hint = Div(
    text="",
    width=SIDEBAR_W - 20,
)
# Both inputs live in the Advanced-settings modal (rarely changed), so
# their width matches the other Advanced inputs (_ADV_INPUT_W).
_ADV_INPUT_W = 240
opt_priority_input = TextInput(
    title="Priority cutoff ≤ (blank = all)", value="", placeholder="e.g. 1",
    width=_ADV_INPUT_W,
)
# Global multi-config cap (v1.4.0). Default 1 → each source is observed in
# at most ONE config (disjoint configs, no duplicate pointing). Blank =
# unlimited. A per-source override (Constraints… popover / max_configs
# column / rule) wins over this default. Only consulted when n_configs > 1.
opt_global_max_configs_input = TextInput(
    title="Max configs per source (blank = unlimited)",
    value="1", placeholder="e.g. 1",
    width=_ADV_INPUT_W,
)

# Collision-protection group — opt-in. When enabled, sources matching
# the priority/weight rule are marked as "protected" and the optimizer
# drops other sources whose spectra would collide with theirs on the
# detector under the current Disperser / Filter (live-canvas orange
# rule, applied per pointing).
opt_protect_section_div = Div(
    text=("<div style='margin:6px 0 -2px 0; font-weight:600; "
          "color:#1a3b66'>Protect spectra from collision</div>"
          "<div style='font-size:11px; color:#5a6b85; line-height:1.35'>"
          "Drop targets whose spectra would overlap the protected set "
          "on the detector under the current Disperser / Filter.</div>"),
    width=SIDEBAR_W - 20,
)
opt_protect_enable_cb = CheckboxGroup(
    labels=["Enable collision protection"], active=[],
    width=SIDEBAR_W - 20,
)
opt_protect_mode_radio = RadioGroup(
    labels=["By priority ≤", "By weight ≥"], active=0,
    inline=True, width=SIDEBAR_W - 20,
)
opt_protect_threshold_input = TextInput(
    title="Threshold", value="1", placeholder="e.g. 1",
    width=SIDEBAR_W - 20,
)
opt_protect_status_div = Div(
    text="<small style='color:#5a6b85'>—</small>",
    width=SIDEBAR_W - 20,
)

# Advanced settings — surfaced via a pop-up modal so they don't bloat
# the Pointing tab. The widgets retain their values regardless of the
# modal's visibility; the optimizer reads them when Run is clicked.
opt_advanced_btn = Button(label="Advanced settings…",
                          button_type="default", width=SIDEBAR_W - 20)

opt_grid_n_ra_input = TextInput(title="Grid n_RA", value="20", width=_ADV_INPUT_W)
opt_grid_n_dec_input = TextInput(title="Grid n_Dec", value="20", width=_ADV_INPUT_W)
opt_grid_n_pa_input = TextInput(title="Grid n_PA", value="20", width=_ADV_INPUT_W)
opt_de_maxiter_input = TextInput(title="DE max iter", value="200", width=_ADV_INPUT_W)
opt_objective_select = Select(
    title="Objective",
    options=["number", "flux"],
    value="number",
    width=_ADV_INPUT_W,
)
opt_sigma_input = TextInput(
    title="Source σ (arcsec)", value="0.06", width=_ADV_INPUT_W,
)
opt_theta_input = TextInput(
    title="APT θ (DVA, deg) — 90 = none", value="90", width=_ADV_INPUT_W,
)
opt_advanced_modal_close_btn = Button(label="Done", button_type="primary",
                                      width=80)
opt_advanced_modal_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)

opt_advanced_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        # Above the optimizer config modal (z-index 1000) — Advanced
        # settings is opened FROM inside the config modal, so it must
        # stack on top of it, otherwise the config card covers the
        # Advanced card. Same reasoning for the card below.
        "z-index": "1001",
    },
)
opt_advanced_modal_card = column(
    row(
        Div(text="<h3>Advanced optimizer settings</h3>",
            sizing_mode="stretch_width"),
        opt_advanced_modal_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    Div(text="<div style='font-size:12px; color:#5a6b85'>"
             "Tune only if the defaults don't fit. Values stick after Done. "
             "<i>Drag the header to reposition.</i></div>",
        width=520),
    row(opt_grid_n_ra_input, opt_grid_n_dec_input, spacing=12),
    row(opt_grid_n_pa_input, opt_de_maxiter_input, spacing=12),
    opt_objective_select,
    opt_sigma_input,
    opt_theta_input,
    Div(text="<div style='font-size:12px; color:#1f4e87; font-weight:600; "
             "margin-top:4px'>Source selection</div>"
             "<div style='font-size:11px; color:#5a6b85'>Priority cutoff "
             "drops sources above the cutoff. Max configs per source caps "
             "how many configs may observe each source (1 = disjoint "
             "configs).</div>", width=520),
    row(opt_priority_input, opt_global_max_configs_input, spacing=12),
    opt_advanced_modal_close_btn,
    spacing=10,
    width=540,
    visible=False,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 18px",
        # Sits above the config modal card (z-index 1000) so it
        # actually overlays it when opened from inside.
        "z-index": "1002",
        "max-height": "85vh",
        "overflow-y": "auto",
    },
)

# ── Catalog editor modal ────────────────────────────────────────────────
# Sortable, editable spreadsheet view of one loaded catalog.
# `_cat_edit_source` is the working copy; on Apply we write it back to
# the underlying Catalog dataclass in state["catalogs"]. Edits are
# captured on `source.on_change("data")` into an undo / redo stack so
# the user can revert mistakes.
cat_edit_select = Select(
    title="Catalog to edit", options=[], value="",
    width=420,
)
# `_idx` carries the per-row row index, used by the trash-icon column's
# HTMLTemplateFormatter so its onclick handler knows which row to
# delete. Maintained on every populate / delete / undo.
_cat_edit_source = ColumnDataSource(data=dict(
    id=[], ra=[], dec=[], priority=[], mag=[], z=[], label=[], _idx=[],
    # Per-target spectral constraints (v1.3.0+). All optional; the
    # initial state for any newly-loaded catalog is "no constraint
    # set". `_has_constraint` is a derived 0/1 flag used by the
    # Constraints-column HTMLTemplateFormatter to pick the button's
    # colour (gray = unset, blue = at least one field set).
    lam_req=[], no_gap=[], extend_blue=[], extend_red=[], protect=[],
    centration=[], max_configs=[],
    _has_constraint=[],
))
# Tiny sink the JS-side delete handler writes into when the user
# clicks a 🗑️ icon. Python's `_on_cat_edit_delete_signal` listens for
# changes and removes the matching row from the working copy.
_cat_edit_delete_signal = ColumnDataSource(data=dict(idx=[-1], stamp=[0]))
# Same pattern for the per-row Constraints… button. The JS-side
# onclick writes the row index here; the Python handler opens the
# constraints popover pre-filled with that row's current values.
_cat_edit_constraint_signal = ColumnDataSource(data=dict(idx=[-1], stamp=[0]))
# Same pattern for the always-visible top-bar CONFIG chip: a JS onclick
# stamps `data`, and `_on_config_chip_signal` advances the active config
# (1→2→…→1) so the user can switch configs without leaving the canvas.
_config_chip_signal = ColumnDataSource(data=dict(stamp=[0]))

_NUM_FMT = NumberFormatter(format="0.[000000]")
# Trash icon column — Underscore.js template renders a clickable span
# whose onclick fires the JS function installed by the DocumentReady
# CustomJS below. Single-click → row removed.
_TRASH_TEMPLATE = (
    "<span style='cursor:pointer; font-size:16px; user-select:none' "
    "      title='Delete this row' "
    "      onclick='window.__vmpt_delete_row(<%= _idx %>)'>"
    "🗑️</span>"
)
# Per-row Constraints… button. The button's background flips between
# gray (no constraint set) and the vMPT primary blue (≥1 constraint
# set) based on the `_has_constraint` flag carried per-row in
# `_cat_edit_source`. Clicking writes the row index into
# `_cat_edit_constraint_signal` via the global JS function installed
# below, and the Python listener opens the popover pre-filled.
_CONSTRAINT_TEMPLATE = (
    "<span style=\"cursor:pointer; padding:1px 8px; border-radius:4px; "
    "font-size:11px; user-select:none; "
    "background:<%= _has_constraint ? '#1f6fc0' : '#e8eaef' %>; "
    "color:<%= _has_constraint ? 'white' : '#5a6b85' %>; "
    "border:1px solid <%= _has_constraint ? '#155a9b' : '#c8d0de' %>;\" "
    "title='Edit per-target spectral constraints for this row' "
    "onclick='window.__vmpt_open_constraints(<%= _idx %>)'>"
    "Edit…</span>"
)
cat_edit_table = DataTable(
    source=_cat_edit_source,
    columns=[
        # Every editable column uses StringEditor. NumberEditor's built-in
        # validator rejects blank input, which means optional columns
        # like Priority / Mag / z that legitimately have no value would
        # trap the user inside the editor. StringEditor accepts anything
        # (including ""); we coerce strings → floats on Apply. Numeric
        # columns are stored as pre-formatted strings in source.data
        # (NaN → ""), so the rendered cell already looks right.
        TableColumn(field="id",       title="ID",       editor=StringEditor(), width=110),
        TableColumn(field="ra",       title="RA (deg)", editor=StringEditor(), width=110),
        TableColumn(field="dec",      title="Dec (deg)",editor=StringEditor(), width=110),
        TableColumn(field="priority", title="Priority", editor=StringEditor(), width=80),
        TableColumn(field="mag",      title="Mag",      editor=StringEditor(), width=80),
        # `z` is stored as float; the float HTML formatter is added by
        # `_cat_edit_rebuild_columns()` (same pattern priority/weight
        # follow with their int formatter). Static init is plain
        # StringEditor with no formatter — only ever visible before the
        # first rebuild, which runs as part of the catalog load.
        TableColumn(field="z",        title="z",        editor=StringEditor(), width=80),
        TableColumn(field="label",    title="Label",    editor=StringEditor(), width=150),
        # Per-row Constraints button. Always visible — the column
        # picker doesn't hide this one because users can't see the
        # constraints any other way. Same `field=_idx` trick as the
        # trash icon: HTMLTemplateFormatter reads `_has_constraint`
        # for the colour and `_idx` for the onclick row index.
        TableColumn(field="_idx",     title="Constraints", width=82,
                    formatter=HTMLTemplateFormatter(template=_CONSTRAINT_TEMPLATE),
                    sortable=False),
        TableColumn(field="_idx",     title="🗑", width=34,
                    formatter=HTMLTemplateFormatter(template=_TRASH_TEMPLATE),
                    sortable=False),
    ],
    editable=True,
    # auto_edit=True maps to SlickGrid's `autoEdit` option, which
    # opens the editor as soon as a cell receives focus. Single-click
    # → click promotes the cell to focus → editor opens. Needs
    # `selectable=True` because SlickGrid's cell-focus model is gated
    # on the same machinery as row selection; setting selectable=False
    # silently disables cell editing too.
    auto_edit=True,
    sortable=True,
    selectable=True,
    reorderable=False,
    width=820, height=380,
    index_position=None,
)
# Column-visibility picker. Pre-populated when a catalog is opened
# with every column the loader saw (the 7 canonical columns plus any
# extras). The user ticks which to show in the table.
cat_edit_columns_choice = MultiChoice(
    title="Show columns",
    value=[], options=[],
    width=600,
)
cat_edit_new_col_input = TextInput(
    title="Add a custom column (e.g. reference, notes)",
    placeholder="column name", width=300,
)
cat_edit_new_col_btn = Button(label="Add column", button_type="default",
                              width=120)
cat_edit_compute_w_btn = Button(
    label="Compute w from p", button_type="default", width=170,
)
cat_edit_compute_p_btn = Button(
    label="Compute p from w", button_type="default", width=170,
)
cat_edit_compute_div = Div(
    text="<small style='color:#5a6b85'>"
         "<b>w ↔ p:</b> the optimizer's Meritocracy uses Weight; "
         "Hierarchy uses Priority. Use these to derive one from the "
         "other.</small>", width=460,
)

# ── Bulk max_configs rule (v1.4.1) ───────────────────────────────────────
# "Set max configs = N for sources where <condition>". The condition is a
# boolean expression over catalog columns (e.g. (mag_f444w > 27) & (z > 6))
# evaluated safely by `evaluate_catalog_condition`; matching rows get the
# chosen value, the rest keep theirs. Manual per-row edits stay available
# via the Constraints… popover. Use it for "faint sources → observe in 2
# configs" style cuts, or "(use global)" to clear an override.
cat_rule_value_select = Select(
    title="Set max configs =", value="2",
    options=["1", "2", "(use global)"], width=150,
)
cat_rule_condition_input = TextInput(
    title="for sources where",
    placeholder="e.g. (mag_f444w > 27) & (z > 6)", width=430,
)
cat_rule_apply_btn = Button(label="Apply rule", button_type="default",
                            width=110)
cat_rule_help_div = Div(
    text="<small style='color:#5a6b85'>Boolean expression over catalog "
         "columns — combine with <code>&amp;</code> / <code>|</code> / "
         "<code>~</code> and parenthesise each test (e.g. "
         "<code>(mag &gt; 27) &amp; (z &gt; 6)</code>). Functions: "
         "abs, log10, log, sqrt, exp, isin. Matching sources are updated; "
         "the rest keep their value. Syntax is checked before "
         "applying.</small>", width=620,
)
cat_rule_status_div = Div(text="", width=620)
# Inject CSS so SlickGrid cells are text-selectable (so the user can
# drag to highlight + Ctrl-C copy / Ctrl-V paste cell content while
# the editor input is active). Bokeh Div renders `<style>` content as
# real stylesheet rules when inserted into the DOM.
_cat_edit_css = Div(text="""
<style>
  .bk-data-table .slick-cell {
    user-select: text !important;
    -webkit-user-select: text !important;
    cursor: text;
  }
  .bk-data-table .slick-cell input {
    user-select: text !important;
    -webkit-user-select: text !important;
  }
  /* Top-right × dismiss buttons on pop-up modals. The Bokeh layout
     can't position the button absolutely, so we mark it with a CSS
     class and float it into the corner via CSS. */
  .vmpt-modal-x button {
    position: absolute;
    top: 6px; right: 8px;
    background: transparent;
    border: 1px solid transparent;
    color: #5a6b85;
    font-size: 20px;
    line-height: 18px;
    padding: 2px 8px;
    cursor: pointer;
    border-radius: 4px;
  }
  .vmpt-modal-x button:hover {
    color: #1a3b66;
    background: rgba(20, 30, 50, 0.06);
    border-color: #c8d0de;
  }
  /* The ⓘ help toggle next to the Method dropdown. Subtler than a
     normal button — looks like a link rather than a button. */
  .vmpt-help-toggle button, .vmpt-help-toggle .bk-btn {
    background: transparent !important;
    border: 1px solid transparent !important;
    box-shadow: none !important;
    color: #5a6b85 !important;
    font-weight: 400 !important;
    text-align: left;
    padding: 2px 0;
    font-size: 12px;
    cursor: pointer;
  }
  .vmpt-help-toggle button:hover, .vmpt-help-toggle .bk-btn:hover {
    color: #1a3b66 !important; text-decoration: underline;
  }
  /* ---- Shared polish for every pop-up dialog (.vmpt-modal-card) ---- */
  .vmpt-modal-header h3 {
    margin: 0; font-size: 15px; font-weight: 600;
    color: #21344f; letter-spacing: .01em;
  }
  /* Text inputs, selects, spinners inside dialogs: consistent, rounded,
     with a clear focus ring. Scoped to modals so the sidebar is untouched. */
  .vmpt-modal-card .bk-input {
    border: 1px solid #cbd5e3; border-radius: 6px;
    padding: 5px 9px; background: #fff;
    transition: border-color .12s, box-shadow .12s;
  }
  .vmpt-modal-card .bk-input:focus {
    border-color: #6f9fd8; outline: none;
    box-shadow: 0 0 0 3px rgba(95, 140, 210, 0.18);
  }
  /* Widget titles (labels above inputs): quieter, consistent. */
  .vmpt-modal-card .bk-input-group > label,
  .vmpt-modal-card label.bk-input-group-text {
    font-size: 11.5px; font-weight: 500; color: #51607a;
  }
  /* Default buttons in dialogs get the same rounded corners. */
  .vmpt-modal-card .bk-btn { border-radius: 6px; }
  /* A subtle section divider/label helper used inside the new dialogs. */
  .vmpt-modal-section {
    font-size: 10px; font-weight: 700; letter-spacing: .05em;
    text-transform: uppercase; color: #8a93a6; margin: 2px 0 0;
  }
</style>
""", width=0, height=0)
cat_edit_undo_btn = Button(label="↶ Undo", button_type="default", width=100)
cat_edit_redo_btn = Button(label="↷ Redo", button_type="default", width=100)
cat_edit_history_div = Div(
    text="<small style='color:#5a6b85'>0 edits</small>", width=160,
)
cat_edit_csv_path_input = TextInput(
    title="Save CSV to", placeholder="/path/to/edited.csv",
    width=380,
)
cat_edit_csv_browse_btn = Button(label="Browse…", button_type="default",
                                 width=80)
cat_edit_csv_save_btn = Button(label="Save as CSV",
                               button_type="default", width=120)
cat_edit_apply_btn = Button(label="Apply changes & close",
                            button_type="primary", width=200)
cat_edit_close_btn = Button(label="Cancel", button_type="default", width=80)
# Top-right × dismiss button. Functionally identical to "Cancel" —
# discards any unapplied edits and closes the modal. Mirrors the
# standard close affordance every user already knows from any
# dialog. Wired in `_close_cat_edit_modal` to share the existing
# close handler.
cat_edit_top_close_btn = Button(
    label="×",
    button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)

cat_edit_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        "z-index": "999",
    },
)
cat_edit_modal_card = column(
    row(
        Div(text="<h3>Edit catalog</h3>",
            sizing_mode="stretch_width"),
        cat_edit_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    Div(text="<div style='font-size:12px; color:#5a6b85'>"
             "Click a column header to sort. Double-click a cell to edit. "
             "Click 🗑️ in a row to delete it. <b>↶ Undo / ↷ Redo</b> "
             "revert / replay your edits. "
             "<b>Apply changes</b> commits to the live catalog so the "
             "eMPT bundle export reflects them; <b>Save as CSV</b> writes "
             "a standalone copy.</div>",
        width=820),
    cat_edit_select,
    cat_edit_columns_choice,
    row(cat_edit_new_col_input, cat_edit_new_col_btn, spacing=10),
    cat_edit_compute_div,
    row(cat_edit_compute_w_btn, cat_edit_compute_p_btn, spacing=10),
    cat_rule_help_div,
    row(cat_rule_value_select, cat_rule_condition_input, cat_rule_apply_btn,
        spacing=10),
    cat_rule_status_div,
    row(cat_edit_undo_btn, cat_edit_redo_btn, cat_edit_history_div,
        spacing=10),
    _cat_edit_css,
    cat_edit_table,
    Div(text="<b>Save as CSV</b>", width=820),
    row(cat_edit_csv_path_input, cat_edit_csv_browse_btn,
        cat_edit_csv_save_btn, spacing=10),
    row(cat_edit_apply_btn, cat_edit_close_btn, spacing=10),
    spacing=10,
    width=860,
    visible=False,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 20px",
        "z-index": "1000",
        "max-height": "92vh",
        "overflow-y": "auto",
    },
)

# ── Per-target Constraints… popover (v1.3.0+) ───────────────────────────
# Opens from inside the catalog editor when the user clicks the
# "Edit…" button in a row's Constraints column. Lets the user set
# `lam_req`, `no_gap`, `extend_blue`, `extend_red`, and `protect`
# for the row whose index is stored in the constraint-click signal.
# All edits stay scoped to the editor's working source until the
# user clicks Apply (writes back into the catalog), matching the
# rest of the editor's "stage then commit" workflow.

# Index of the row currently being edited. Updated by
# `_on_cat_edit_constraint_signal` when the user clicks Edit… on a
# row; reset to -1 when the popover closes.
_cat_constraints_row_idx: int = -1

cat_constraints_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)
cat_constraints_title_div = Div(
    text="<div style='font-size:12px; color:#5a6b85'>"
         "These constraints apply only when the optimizer evaluates "
         "this source. Empty / unchecked = no constraint.</div>",
    width=420,
)
cat_constraints_row_label = Div(
    text="<small style='color:#5a6b85'>Editing row —</small>",
    width=420,
)
cat_constraints_lam_input = TextInput(
    title="Required λ ranges (μm; semicolon-separated)",
    placeholder="e.g. 1.0-1.3; 1.5-1.8",
    value="",
    width=420,
)
cat_constraints_lam_warn = Div(
    text="", width=420,
)
cat_constraints_checks = CheckboxGroup(
    labels=[
        "Forbid detector gap inside spectrum (no_gap)",
        "Extend to bluest λ of disperser",
        "Extend to reddest λ of disperser",
        "Protect this source from spectral collision",
    ],
    active=[],
    width=420,
)
# Per-target source-centering override (v1.3.1+). "(use global)" → ""
# in storage; the rest are the same five labels as the optimizer's
# global Source-centering Select. The override is **unconditional** —
# whatever the user picks here wins, even when it's laxer than the
# global setting.
cat_constraints_centration_select = Select(
    title="Source centering override (blank = use optimizer global)",
    options=[
        "(use global)",
        "UNCONSTRAINED",
        "ENTIRE_OPEN",
        "MIDPOINT",
        "CONSTRAINED",
        "TIGHTLY_CONSTRAINED",
    ],
    value="(use global)",
    width=420,
)
cat_constraints_centration_hint = Div(
    text="<small style='color:#5a6b85'>"
         "Wins over the optimizer's global Source-centering "
         "setting for this row only.</small>",
    width=420,
)
# Per-target multi-config cap (v1.4.0+). "(use global)" → "" in storage;
# "1"/"2" cap how many MPT configs the optimizer may place this source in.
cat_constraints_max_configs_select = Select(
    title="Max MPT configs this source may be observed in",
    options=["(use global)", "1", "2"],
    value="(use global)",
    width=420,
)
cat_constraints_max_configs_hint = Div(
    text="<small style='color:#5a6b85'>"
         "Caps how many of the planned MPT configs may observe this "
         "source. Blank = use the optimizer's global default.</small>",
    width=420,
)
cat_constraints_apply_btn = Button(
    label="Apply", button_type="primary", width=80,
)
cat_constraints_cancel_btn = Button(
    label="Cancel", button_type="default", width=80,
)
cat_constraints_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        # Stacks above the catalog editor's own backdrop (z-index 999)
        # so the popover doesn't sit behind it.
        "z-index": "1010",
    },
)
cat_constraints_modal_card = column(
    row(
        Div(text="<h3>Per-target spectral constraints</h3>",
            sizing_mode="stretch_width"),
        cat_constraints_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    cat_constraints_title_div,
    cat_constraints_row_label,
    cat_constraints_lam_input,
    cat_constraints_lam_warn,
    cat_constraints_checks,
    cat_constraints_centration_select,
    cat_constraints_centration_hint,
    cat_constraints_max_configs_select,
    cat_constraints_max_configs_hint,
    row(cat_constraints_apply_btn, cat_constraints_cancel_btn, spacing=10),
    spacing=10,
    width=460,
    visible=False,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 18px",
        "z-index": "1011",
        "max-height": "85vh",
        "overflow-y": "auto",
    },
)


# Width leaves room for the 80 px Cancel button (+10 px spacing) inside
# the modal body (SIDEBAR_W + 30 = 370 px); 270 + 10 + 80 = 360 ≤ 370 so
# the Cancel button no longer overflows the card and gets clipped.
opt_run_btn = Button(label="Run optimization",
                    button_type="primary", width=SIDEBAR_W - 70)
# The Pointing-tab CTA that opens the optimizer's config modal. The
# config widgets used to live inline in the Pointing tab; on a 913 px
# laptop screen they pushed the actual Run button below the fold.
# Wrapping them in a modal keeps the Pointing tab compact.
opt_open_btn = Button(
    label="Open optimizer…",
    button_type="primary",
    width=SIDEBAR_W - 20,
)

# ── MPT configurations (v1.4.0) ─────────────────────────────────────────
# Multiple MPT configs mirror JWST APT/MPT: each config is an independent
# pointing + set of open shutters (a separate exposure). The user works on
# one config at a time ("Working on"); manual shutter opens land only on
# the active config. Default is a single config so v1.3.x behaviour is
# unchanged until the user opts in.
mpt_num_configs_spinner = Spinner(
    low=1, high=_MAX_CONFIGS, step=1, value=1, width=_HALF_W,
    title="Number of configs",
)
mpt_config_select = Select(
    title="Working on", value="Config 1", options=["Config 1"],
    width=_HALF_W,
)
# Prominent "you are here" banner — only shown in multi-config mode so the
# user is never confused about which config a manual open lands in. Updated
# by `_refresh_active_config_banner()` on every switch / count change.
mpt_active_config_div = Div(text="", width=SIDEBAR_W - 20, visible=False)
mpt_view_btn = Button(label="View MPT catalog…",
                      button_type="default", width=SIDEBAR_W - 20)

# Read-only "MPT catalog viewer" — one row per (config, selected source).
# Mirrors the input catalog editor's DataTable but is never editable; it
# stays empty until a shutter is opened by hand or by the optimizer.
# Numeric columns are STORED as floats (NaN = missing) so the DataTable
# sorts them numerically — string storage made "101" sort before "11".
# These NaN-safe Underscore templates render them (integer for
# Pri/Weight/Q/s/d, fixed decimals elsewhere) and blank a missing value.
def _mpt_num_template(decimals: int) -> str:
    render = ("Math.round(value)" if decimals == 0
              else f"value.toFixed({decimals})")
    return (
        "<%= (value === null || value === undefined || value === '' || "
        "(typeof value === 'number' && isNaN(value))) ? '' : "
        f"(typeof value === 'number' ? {render} : value) %>"
    )


_MPT_FMT_INT = HTMLTemplateFormatter(template=_mpt_num_template(0))
_MPT_FMT_COORD = HTMLTemplateFormatter(template=_mpt_num_template(6))
_MPT_FMT_Z = HTMLTemplateFormatter(template=_mpt_num_template(4))
_MPT_FMT_MAG = HTMLTemplateFormatter(template=_mpt_num_template(2))
_MPT_FMT_UM = HTMLTemplateFormatter(template=_mpt_num_template(3))  # μm, 3 dp
# field → formatter for the numeric columns (others render as plain text).
_MPT_VIEW_NUM_FMT = {
    "ra": _MPT_FMT_COORD, "dec": _MPT_FMT_COORD,
    "priority": _MPT_FMT_INT, "weight": _MPT_FMT_INT,
    "q": _MPT_FMT_INT, "s": _MPT_FMT_INT, "d": _MPT_FMT_INT,
    "z": _MPT_FMT_Z, "mag": _MPT_FMT_MAG,
    "lam_blue": _MPT_FMT_UM, "lam_red": _MPT_FMT_UM,
}
_mpt_view_source = ColumnDataSource(data=dict(
    config=[], id=[], ra=[], dec=[], priority=[], weight=[], mag=[], z=[],
    lam_blue=[], lam_red=[], gap=[], q=[], s=[], d=[], label=[],
))
# Toggleable columns (Cfg + Source ID are always shown). (field, title, width).
# "Role" (the internal target/sky nod-shutter attribute) is intentionally
# NOT shown — it isn't a per-source property and was misleading.
# λ_blue / λ_red / Gap are the spectrum's on-detector wavelength coverage
# and the NRS1/NRS2 detector-gap range for the source's shutter under the
# current Disperser / Filter (computed on open via wavelengths.cutoffs).
_MPT_VIEW_COLS = [
    ("ra", "RA (deg)", 92), ("dec", "Dec (deg)", 92),
    ("priority", "Pri", 44), ("weight", "Weight", 66), ("mag", "Mag", 56),
    ("z", "z", 60),
    ("lam_blue", "λ_blue", 64), ("lam_red", "λ_red", 64), ("gap", "Gap (μm)", 104),
    ("q", "Q", 34), ("s", "s", 42), ("d", "d", 42),
    ("label", "Label", 120),
]
_MPT_VIEW_TITLE_TO_FIELD = {t: f for f, t, _ in _MPT_VIEW_COLS}
# Colour the Cfg cell to match the active-config chip language
# (Config 1 → blue, Config 2 → magenta; any further config → green) so
# the per-config grouping reads at a glance in the multi-config viewer.
_MPT_VIEW_CFG_FORMATTER = HTMLTemplateFormatter(template=(
    "<span style='font-weight:600; color:"
    "<%= (value==1)?'#1f6fc0':(value==2)?'#b5179e':'#2a9d3a' %>'>"
    "<%= value %></span>"
))
mpt_view_columns_choice = MultiChoice(
    title="Show columns (drag chips to reorder; remove to hide)",
    value=[t for _, t, _ in _MPT_VIEW_COLS],
    options=[t for _, t, _ in _MPT_VIEW_COLS],
    width=1040,
)
mpt_view_table = DataTable(
    source=_mpt_view_source,
    columns=[
        TableColumn(field="config", title="Cfg", width=42,
                    formatter=_MPT_VIEW_CFG_FORMATTER),
        TableColumn(field="id", title="Source ID", width=104),
        *[TableColumn(field=f, title=t, width=w,
                      **({"formatter": _MPT_VIEW_NUM_FMT[f]}
                         if f in _MPT_VIEW_NUM_FMT else {}))
          for f, t, w in _MPT_VIEW_COLS],
    ],
    editable=False,
    sortable=True,
    selectable=True,
    reorderable=True,   # drag a row's handle to reorder
    index_position=0,   # drag-handle / row-number gutter
    width=1080,
    height=360,
    autosize_mode="none",
)


def _mpt_view_rebuild_columns() -> None:
    """Rebuild the viewer's columns from the picker. Cfg + Source ID are
    always present; the chosen optional columns follow in the picker's
    (drag-reorderable) order."""
    chosen = list(mpt_view_columns_choice.value or [])
    cols = [
        TableColumn(field="config", title="Cfg", width=42,
                    formatter=_MPT_VIEW_CFG_FORMATTER),
        TableColumn(field="id", title="Source ID", width=104),
    ]
    width_of = {t: w for _, t, w in _MPT_VIEW_COLS}
    for title in chosen:
        f = _MPT_VIEW_TITLE_TO_FIELD.get(title)
        if f:
            kw = ({"formatter": _MPT_VIEW_NUM_FMT[f]}
                  if f in _MPT_VIEW_NUM_FMT else {})
            cols.append(TableColumn(field=f, title=title,
                                    width=width_of.get(title, 80), **kw))
    mpt_view_table.columns = cols


mpt_view_columns_choice.on_change(
    "value", lambda a, o, n: _mpt_view_rebuild_columns())
mpt_view_summary_div = Div(text="", width=1040)
mpt_view_top_close_btn = Button(label="×", button_type="default",
                                width=32, height=28,
                                css_classes=["vmpt-modal-x"])
mpt_view_close_btn = Button(label="Close", button_type="default", width=90)
mpt_view_csv_path_input = TextInput(
    title="Save selected list as CSV", value="", width=620,
    placeholder="/path/to/mpt_selected.csv",
)
mpt_view_csv_save_btn = Button(label="Save as CSV",
                               button_type="default", width=110)
mpt_view_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        "z-index": "999",
    },
)
mpt_view_modal_card = column(
    row(
        Div(text="<h3>MPT catalog — selected sources</h3>",
            sizing_mode="stretch_width"),
        mpt_view_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    mpt_view_summary_div,
    mpt_view_columns_choice,
    mpt_view_table,
    row(mpt_view_csv_path_input, mpt_view_csv_save_btn, spacing=10),
    row(mpt_view_close_btn, spacing=10),
    spacing=10,
    width=1100,
    visible=False,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 20px",
        "z-index": "1000",
        "max-height": "92vh",
        "overflow-y": "auto",
    },
)
opt_status_div = Div(
    # Updated live by `_refresh_opt_status_div` from the catalog +
    # method state. Default is a neutral placeholder for when there's
    # no catalog at all.
    text=("<small style='color:#5a6b85'>"
          "Load a catalog (Input tab) before running.</small>"),
    width=SIDEBAR_W - 20,
)
# Config modal: wraps every optimizer input + the Run button so the
# Pointing tab can stay short. Built later (after all widgets are
# defined); see `opt_config_modal_card` and `_open_opt_config_modal`.
opt_config_close_btn = Button(
    label="Cancel",
    button_type="default", width=80,
)
opt_config_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)
opt_config_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        "z-index": "999",
    },
)

# The actual config modal body. All the optimizer-config widgets
# (Method, ΔRA/ΔDec/ΔPA, Refine top N, Source centering, Priority
# cutoff, Protect-spectra group, Advanced settings…, Run optimization,
# status) live here. The Pointing tab now just shows a single
# "Open optimizer…" button that flips this card's `visible`.
opt_config_modal_card = column(
    row(
        Div(text="<h3>MSA pointing optimizer</h3>",
            sizing_mode="stretch_width"),
        opt_config_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    Div(text="<div style='font-size:12px; color:#5a6b85'>"
             "Searches (RA, Dec, V3 PA) within the box below for the "
             "best placement. Adjust the inputs, then <b>Run</b>.</div>",
        width=SIDEBAR_W + 30),
    opt_method_select,
    opt_method_help_toggle,
    opt_method_help_div,
    row(opt_dra_input, opt_ddec_input, spacing=12),
    row(opt_dpa_input, opt_n_top_input, spacing=12),
    opt_centration_select,
    opt_centration_override_hint,
    opt_protect_section_div,
    opt_protect_enable_cb,
    opt_protect_mode_radio,
    opt_protect_threshold_input,
    opt_protect_status_div,
    opt_advanced_btn,
    row(opt_run_btn, opt_config_close_btn, spacing=10),
    opt_status_div,
    spacing=10,
    width=SIDEBAR_W + 70,
    visible=False,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 18px",
        "z-index": "1000",
        "max-height": "92vh",
        "overflow-y": "auto",
    },
)
# Results live in a column that's rebuilt after every run. Each row
# is a Button labelled with the rank + score + delta-pointing.
# (Kept for callers that prefer an inline list; the primary surface
# is now the modal dialog below.)
opt_results_column = column(width=SIDEBAR_W - 20)

# ── Optimizer pop-up dialog ──────────────────────────────────────────────
# Two-phase modal: progress (live progress while grid + DE run) and
# results (top-N candidates with Apply buttons). Realised as a
# position-fixed Bokeh column overlaying the page; a sibling Div
# renders the semi-transparent backdrop.

opt_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.45)",
        "z-index": "999",
    },
)

# Dialog header.
opt_modal_title = Div(
    text="<h3 style='margin:0 0 6px 0; color:#1a3b66'>"
         "MSA pointing optimization</h3>",
    width=560,
)

# Progress section (shown while running).
# Spinner Div — text is set ONCE at construction and never updated,
# so the CSS animation keeps spinning continuously. (Bokeh's Div.text
# update replaces innerHTML, which would restart any child animation.)
opt_modal_progress_spinner = Div(
    text=(
        "<style>"
        "@keyframes vmpt-spin { from { transform: rotate(0deg); }"
        "                       to   { transform: rotate(360deg); } }"
        "@keyframes vmpt-stripe { 0%   { background-position: 0 0; }"
        "                         100% { background-position: 32px 0; } }"
        "@keyframes vmpt-pulse { 0%, 100% { box-shadow: 0 0 8px rgba(50,115,220,0.35); }"
        "                        50%      { box-shadow: 0 0 18px rgba(50,115,220,0.7); } }"
        "</style>"
        "<div style='display:inline-block; width:18px; height:18px;"
        " border:3px solid #c9d4e8; border-top-color:#3273dc;"
        " border-radius:50%;"
        " animation: vmpt-spin 0.85s linear infinite;"
        " vertical-align:middle;'></div>"
    ),
    width=26, height=26,
)
opt_modal_progress_text = Div(
    text="<i>Starting…</i>", width=520,
    styles={"font-size": "13px", "padding": "4px 0", "line-height": "24px"},
)
# Static bar HTML — set ONCE. The inner fill uses `width: var(--vmpt-pct)`,
# so updating the wrapper's `styles` dict (which Bokeh applies without
# replacing innerHTML) is enough to change the fill width. This keeps
# the stripe animation running continuously without restarts.
_BAR_HTML = (
    "<div style='background:linear-gradient(180deg, #dde3ec 0%, #eaeff7 100%);"
    " border-radius:8px; height:16px; overflow:hidden;"
    " box-shadow:inset 0 1px 3px rgba(0,0,0,0.07);'>"
    "<div style='"
    " width: var(--vmpt-pct, 0%);"
    " height: 100%;"
    " background-color: #3273dc;"
    " background-image: linear-gradient(135deg,"
    "   rgba(255,255,255,0.30) 25%, transparent 25%,"
    "   transparent 50%, rgba(255,255,255,0.30) 50%,"
    "   rgba(255,255,255,0.30) 75%, transparent 75%);"
    " background-size: 32px 32px;"
    " animation: vmpt-stripe 0.8s linear infinite,"
    "            vmpt-pulse  2.2s ease-in-out infinite;"
    " transition: width 0.3s ease-out;"
    " border-radius: 8px;"
    "'></div></div>"
)
opt_modal_progress_bar = Div(
    text=_BAR_HTML,
    width=560,
    styles={"--vmpt-pct": "0%"},
)
opt_modal_progress_box = column(
    row(opt_modal_progress_spinner, opt_modal_progress_text, spacing=8),
    opt_modal_progress_bar,
    spacing=4,
    width=560,
)

# Results section (shown when done). Built as `column(header, row1,
# row2, …)` where each result row is itself `row(cells_div, apply_btn)`
# so the Apply button lines up natively with its cells — the previous
# parallel "HTML table + buttons column" pattern drifted out of
# alignment because Bokeh column spacing accumulated between buttons.
opt_modal_results_summary = Div(text="", width=820)
opt_modal_results_rows = column(spacing=0, width=820)
opt_modal_results_box = column(
    opt_modal_results_summary,
    opt_modal_results_rows,
    spacing=4,
    width=820,
    visible=False,
)

opt_modal_close_btn = Button(label="Close", button_type="default", width=80)
opt_modal_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)

opt_modal_card = column(
    row(
        Div(text="<h3>MSA pointing optimization</h3>",
            sizing_mode="stretch_width"),
        opt_modal_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    opt_modal_progress_box,
    opt_modal_results_box,
    opt_modal_close_btn,
    visible=False,
    spacing=10,
    # Sized to fit the widest table — the multi-config combined view
    # (rank + combined + per-config Cfg/Δ/Score + Apply ≈ 800 px) plus
    # the modal's inner padding (~36 px) — so neither the Score column
    # nor the Apply button needs horizontal scrolling.
    width=860,
    css_classes=["vmpt-modal-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 18px",
        "z-index": "1000",
        "min-width": "700px",
        "max-width": "94vw",
        # `max-height` + overflow-y keeps the card from running off
        # the bottom of the viewport when there are many result rows.
        "max-height": "85vh",
        "overflow-y": "auto",
    },
)

# In-flight optimizer state — only one run at a time. Reset on each
# new run via `on_optimize`.
_opt_run: dict = {}

# Hidden trigger used by the Apply-button JS confirm dialog. The JS
# callback writes "<ra>,<dec>,<pa>,<stamp>" here once the user OKs
# the confirm; Python's on_change handler parses + applies. The
# stamp guarantees a fresh `change` event even when the same row is
# re-applied (Bokeh dedupes identical values).
opt_apply_trigger = TextInput(value="", visible=False)
# v1.4.0 (N configs since v1.5.0): applies ALL configs of a multi-config
# plan in one click. The JS confirm sets a fresh stamp; the Python handler
# reads each config's best pointing out of `_opt_run["pass_results"]` and
# applies it to that config.
opt_apply_both_trigger = TextInput(value="", visible=False)

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
help_toggle_btn = Button(label="Show help", button_type="default", width=110)

# ── Rotating tip strip ───────────────────────────────────────────────────────
# The help panel shows a rotating one-liner tip on top, and the full
# reference below it. Tips fade in/out every 15 s so the help feels alive
# without being annoying. See `_TIPS` for content; rotation is wired by
# `_advance_tip` registered as a periodic callback further down.
_TIPS = [
    ("🎯", "Pick mode", "Click anywhere on the image — vMPT snaps to the nearest operable shutter and opens an <b>N-shutter slitlet</b> (set N=1/2/3/5 in the <b>Setting</b> tab)."),
    ("⌨️", "Pan with the keyboard", "Glide the view with <b>W&nbsp;A&nbsp;S&nbsp;D</b> or the <b>arrow keys</b> — like dragging, but hands stay on the keys. Hold <b>Shift</b> for bigger steps. (Pauses while you're typing in a field.)", True),
    ("🔲", "Space toggles one shutter", "Hover any operable shutter and tap <b>Space</b> to open or close <b>just that shutter</b> — independent of the N-shutter slitlet size. Handy for fine-tuning a mask cell by cell.", True),
    ("✋", "Move the pointing", "<b>Shift + click</b> anywhere on the image to recentre the pointing on that spot. The <span style='color:#2e9b3f;font-weight:600'>lime cross</span> marks the current pointing."),
    ("🔁", "Toggle a slitlet", "Click an already-open shutter to close it. Its slitlet siblings come down with it."),
    ("🔭", "Pick a roll", "In the <b>Pointing</b> tab, enter a visibility date and click <b>Compute allowed V3 PA</b>. jwst_gtvt reports the valid window for the date."),
    ("🌈", "Wavelength check", "Hover any open shutter to see its λ<sub>blue</sub> / λ<sub>red</sub> and the NRS1 / NRS2 detector-gap range for the current disperser."),
    ("⚠️", "Orange = collision", "Orange-tinted shutters share a dispersed-y row with an open or stuck-open shutter — opening them would put two spectra on the same detector pixels."),
    ("🪞", "Cross-quadrant", "Spec-overlap correctly pairs Q1↔Q3 (NRS1) and Q2↔Q4 (NRS2). A pick in Q1 will never light up Q2 or Q4."),
    ("💎", "Catalog match", "Open a shutter with a catalog source inside it — vMPT auto-tags the slitlet with that source's ID. Status bar names the match."),
    ("📤", "Export bundle", "<b>MPT</b> tab → <b>Export eMPT bundle</b> writes a folder with <code>MPT_plan.json</code>, an APT-importable <code>.cat</code> target list, and the eMPT pipeline's three files."),
    ("⏪", "Undo picks", "<b>Setting</b> tab → <b>Undo last</b> reverts the most recent slitlet open/close action. History is 50 deep."),
    ("📐", "Slitlet sizes", "N=2 means clicked-shutter + one row of lower-y on the detector. N=3/5 are centred on the click. Switch any time in <b>Setting</b>."),
    ("🛰️", "Two ways to load APT", "<b>MPT</b> tab → either point at a local <code>.aptx</code>, or just type a JWST program ID (e.g. <code>1208</code>) and vMPT pulls it from STScI."),
    ("🚀", "Pre-load via run.sh", "Start the app with files ready: <code>./run.sh --fits img.fits --catalog tgts.csv</code>. Use <code>--jpg + --wcs</code> for JPG/sidecar. <code>--port 5010</code> picks a different port."),
    ("🧮", "Optimize MSA pointing", "<b>Pointing</b> tab → bottom panel. Set ΔRA/ΔDec/ΔPA (zero on any axis = freeze it) and click Run. Get the top 10 (RA, Dec, V3 PA) ranked by sources placed."),
    ("🔀", "Plan up to 5 configs", "Need more than one pointing? Set <b>Number of configs</b> (<b>Pointing</b> → MPT configurations) up to 5. Each gets its own colour, and the optimizer fills them in turn so they cover different sources, not duplicates.", True),
    ("📝", "Edit catalog inline", "<b>Input</b> tab → <b>Edit catalog…</b>. Double-click any cell to edit. Click 🗑️ on a row to delete it. ↶ Undo / ↷ Redo revert mistakes. Save as CSV or commit back to the live catalog."),
    ("🎴", "Layer multiple catalogs", "Click <b>Add</b> in the Input tab to layer several catalogs. Each gets its own colour. ▲ / ▼ reorder the stack; ✕ removes one; checkbox toggles visibility."),
    ("🖼️", "Pixel-perfect canvas", "Image pixels are always rendered 1:1 — resizing the window letterboxes around the canvas instead of stretching the image. Aspect lock is set automatically when you load."),
    ("🆔", "Big IDs auto-shrink", "Catalog IDs ≥ 10⁷ are taken mod 10⁷ on load — APT MPT wants compact integers. The original string token survives in the Label column for traceability."),
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
    the fade-in animation. A tip may carry an optional 4th element — a truthy
    `is_new` flag — which adds a red NEW pill to its header."""
    tip = _TIPS[idx % len(_TIPS)]
    emoji, header, body = tip[0], tip[1], tip[2]
    is_new = len(tip) > 3 and tip[3]
    badge = (
        '<span style="background:#e8453c; color:#fff; font-size:8.5px; '
        'font-weight:800; letter-spacing:0.5px; padding:1px 5px; '
        'border-radius:8px; margin-left:7px; vertical-align:middle;">NEW</span>'
        if is_new else ''
    )
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
        f'      Tip · {header}{badge}'
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
  /* Quick guide is a stack of collapsible <details> folds. Scoped to
     .vmpt-help; `max-width:100%` + break rules keep every element inside
     the 340-px panel so nothing spills past the box edge. */
  .vmpt-help { overflow-wrap: anywhere; word-break: break-word; }
  .vmpt-help * { max-width: 100%; box-sizing: border-box; }
  .vmpt-help h3 { margin: 0 0 8px 0; font-size: 14px; }
  .vmpt-help b  { color: #1a3b66; }
  .vmpt-help details { margin: 0 0 4px 0; border: 1px solid #e3e6ea;
                       border-radius: 5px; background: #fff; }
  .vmpt-help summary { cursor: pointer; padding: 5px 8px; font-weight: 700;
                       color: #1a3b66; border-radius: 5px;
                       -webkit-user-select: none; user-select: none; }
  .vmpt-help summary:hover { background: #eef2f7; }
  .vmpt-help details[open] > summary { border-bottom: 1px solid #e3e6ea;
                       border-radius: 5px 5px 0 0; }
  .vmpt-help .fold-body { padding: 2px 10px 6px 10px; }
  .vmpt-help ul { margin: 2px 0 4px 0; padding-left: 15px; }
  .vmpt-help ul ul { margin: 1px 0; padding-left: 12px; }
  .vmpt-help li { margin: 2px 0; }
  .vmpt-help code { font-size: 11px; padding: 0 2px;
                    background: #ececec; border-radius: 2px;
                    word-break: break-all; }
</style>
<div class="vmpt-help">
<h3>Quick guide</h3>

<details open>
  <summary>1 · Load an image</summary>
  <div class="fold-body"><ul>
    <li>One-click <b>Load Abell 370 example</b> or <b>Load RXCJ0600 example</b> from the <b>Input</b> tab — fastest.</li>
    <li>Or a local <b>FITS</b> path (with WCS), or a <b>JPG + sidecar FITS</b> pair.</li>
  </ul></div>
</details>

<details>
  <summary>2 · (Optional) load target catalog</summary>
  <div class="fold-body"><ul>
    <li>CSV / ASCII / FITS with at least <code>ID, RA, DEC</code>.</li>
    <li>Targets render as yellow circles. A shutter containing a catalog source auto-tags the slitlet on click.</li>
  </ul></div>
</details>

<details>
  <summary>3 · Aim the MSA</summary>
  <div class="fold-body"><ul>
    <li><b>V3 PA</b> drives the math; <b>NIRSpec APA</b> = V3 PA + 138.575° (mod 360).</li>
    <li><b>Shift + click</b> to move pointing. The <span style='color:#2e9b3f;font-weight:600'>lime cross</span> marks it.</li>
    <li>Type a date in <b>Visibility</b> → <b>Compute allowed V3 PA</b> to query jwst_gtvt.</li>
  </ul></div>
</details>

<details>
  <summary>4 · Optimize MSA</summary>
  <div class="fold-body"><ul>
    <li><b>Pointing</b> tab → <b>Open optimizer…</b>. Pick a method: <b>Democracy</b> (count), <b>Meritocracy</b> (Σ weight), <b>Hierarchy</b> (priority tiers).</li>
    <li>Set the search radius + optional <b>collision protection</b>, then <b>Run</b>. Top solutions appear in a table; <b>Apply #N</b> sets the pointing and auto-opens slitlets.</li>
    <li>Multiple exposures? Bump <b>Number of configs</b> (≤ 5); <b>Max configs per source = 1</b> keeps later configs off earlier targets.</li>
  </ul></div>
</details>

<details>
  <summary>5 · Hand-pick shutters</summary>
  <div class="fold-body"><ul>
    <li>Pick the <b>N-shutter slitlet</b> size (1/2/3/5) in <b>Settings</b>.</li>
    <li><b>Click</b> → opens an N-shutter slitlet at the nearest operable shutter. Click an open shutter to close the slitlet.</li>
    <li><b>Hover + Space</b> → toggle <b>just the one</b> hovered shutter open/close (ignores the slitlet size; operable shutters only).</li>
  </ul></div>
</details>

<details>
  <summary>Layers &amp; colours</summary>
  <div class="fold-body">
  <ul>
    <li><span style='background:silver;padding:0 4px'>silver</span> = operable</li>
    <li><span style='color:#d63d3d;font-weight:700'>red fill</span> = your picks</li>
    <li><span style='color:#b30000;font-weight:700'>dark red</span> = stuck-open</li>
    <li><span style='color:gold;font-weight:700'>gold</span> = fixed slits</li>
    <li><span style='color:#ddd200;font-weight:700'>yellow ○</span> = catalog target · <span style='color:#2e9b3f;font-weight:700'>green ○</span> = matched</li>
  </ul>
  <b>Spec-overlap</b> (MPT colours; alpha stacks with #sources):
  <ul>
    <li><span style='color:#d96272;font-weight:700'>pink</span> = <b>Mask Stuck</b> — operable shutter where only stuck-open spectra land (no collision).</li>
    <li><span style='color:#e26a00;font-weight:700'>orange</span> = <b>Masked</b> — operable shutter where at least one user-open's spectrum lands (no collision).</li>
    <li><span style='color:#a050b8;font-weight:700'>purple</span> = <b>Mask Conflict</b> — two open slitlets crowd with no operable buffer row between them. Only the rows where they crowd go purple (the ±2-row window between the two slitlets), not the whole spectrum; the band beyond reverts to orange.</li>
  </ul>
  </div>
</details>

<details>
  <summary>Save / share / export</summary>
  <div class="fold-body"><ul>
    <li><b>MPT</b> tab → <b>Save session</b> writes a bundle.</li>
    <li><b>Load session</b> — point at <code>MPT_plan.json</code> or <code>vMPT_workspace.json</code>; the sibling auto-loads.</li>
    <li><b>Export eMPT bundle</b> writes a timestamped folder:
      <ul>
        <li><code>MPT_plan.json</code> + <code>&lt;catalog&gt;.cat</code> → APT MPT</li>
        <li><code>vMPT_workspace.json</code> → vMPT round-trip</li>
        <li><code>eMPT_*</code> three files → eMPT pipeline</li>
      </ul>
    </li>
  </ul></div>
</details>

<details>
  <summary>Pan / zoom / shortcuts</summary>
  <div class="fold-body"><ul>
    <li><b>Wheel</b>: zoom · <b>Drag</b> or <b>W A S D</b> / <b>arrows</b>: pan (<b>Shift</b> = bigger steps) · <b>Box zoom</b>: toolbar → drag</li>
    <li><b>Hover + Space</b>: toggle the single hovered shutter · <b>Reset</b>: toolbar · <b>Undo</b>: Settings → <b>Undo last</b></li>
  </ul></div>
</details>

</div>
<p style='margin:6px 2px 0 2px'>📖 Full documentation at <a href='https://vmpt.readthedocs.io/' target='_blank' rel='noopener' style='color:#1a3b66; font-weight:600;'>vmpt.readthedocs.io</a></p>
""",
)


def on_help_toggle():
    showing = not help_div.visible
    help_div.visible = showing
    tip_div.visible = showing
    help_toggle_btn.label = "Hide help" if showing else "Show help"
    # Resize the column itself so the help panel actually gives back
    # horizontal real-estate to the figure column when collapsed. The
    # figure itself has fixed `frame_width`/`frame_height`, so the
    # canvas pixel aspect doesn't change — only the empty space around
    # it grows / shrinks.
    help_panel.width = HELPPANEL_W if showing else 130


# Expanded by default (v1.3.3+) — first-run users get the Quick
# guide + rotating tip without having to discover the "Show help"
# button. Returning users can collapse the panel with one click;
# the button label flips to "Hide help" when expanded.
help_div.visible = True
tip_div.visible = True
help_toggle_btn.label = "Hide help"
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
    # These must match the keys in `_OVERLAY_LAYER_CONFIG` exactly, or the
    # alpha/stroke sliders find no config and silently do nothing. (The
    # v1.3.1 split of "Overlapping shutters" into the three MPT-style
    # spec-overlap layers had left a single dead "Overlapping shutters"
    # entry here, so those sliders were inert.)
    options=[
        "Operable shutters",
        "Mask Stuck (pink)",
        "Masked (overlapping warning)",
        "Mask Conflict (purple)",
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
    start=0.0, end=3.0, step=0.05, value=1.0,
    title="Stroke (px)",
    width=SIDEBAR_W - 40,
)

# Canvas size — two independent sliders for the X (width) and Y
# (height) of the figure's drawing frame in pixels.
# `refresh_image_glyph()` reads state["frame_x"] / state["frame_y"]
# and sets fig.frame_width / fig.frame_height directly.
#
# `match_aspect=True` on the figure enforces the **per-pixel
# square** constraint: 1 data unit in x renders at the same screen
# size as 1 data unit in y, no matter what the canvas X/Y are.
# That's what keeps the science image's pixels visually square AND
# the NIRSpec shutters at their correct geometric ratios. Bokeh
# achieves this by EXPANDING whichever data range is "too short"
# given the canvas dimensions — when the user sets a non-image-
# aspect canvas, the image stays at its native pixel shape and the
# extra canvas area shows empty space (the user can pan into it).
# We deliberately do NOT couple X/Y to the image's W:H — letting
# them roam free is the point of having two sliders. Default 800x600.
# Compact `W: [n] x H: [n]` pair — Spinner instead of Slider. The
# Spinner has up/down arrows, accepts free typing, and only commits
# on blur / Enter / arrow click so it's naturally throttled (no
# need for `value_throttled`).  Layout-wise the two spinners go
# in a `row()` next to a short label so the Settings tab stays
# compact.
_CANVAS_SIZE_MIN, _CANVAS_SIZE_MAX, _CANVAS_SIZE_STEP = 400, 1600, 50
canvas_x_spinner = Spinner(
    low=_CANVAS_SIZE_MIN, high=_CANVAS_SIZE_MAX, step=_CANVAS_SIZE_STEP,
    value=800, width=88, title="Width (X)",
)
canvas_y_spinner = Spinner(
    low=_CANVAS_SIZE_MIN, high=_CANVAS_SIZE_MAX, step=_CANVAS_SIZE_STEP,
    value=600, width=88, title="Height (Y)",
)

undo_btn = Button(label="Undo last", button_type="default")
clear_btn = Button(label="Clear open", button_type="warning")
# Reset-to-defaults button — handler wired further down in the
# preferences-init section (forward-declared here so the Settings
# tab layout can include it).
reset_prefs_btn = Button(
    label="Reset display to defaults",
    button_type="default", width=SIDEBAR_W - 20,
)

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


# ─────────────────────────────────────────────────────────────────────
# MPT-tab pop-up dialogs — Import / Save / Export
# ─────────────────────────────────────────────────────────────────────
# The MPT tab is just three buttons that open these dialogs, so the
# sidebar stays short and the three workflows (bring something IN / save
# a vMPT session / write the APT deliverables) read as distinct actions.
# All three reuse the shared `.vmpt-modal-*` styling. Browse-only: each
# Browse button fills the path field shown beneath it (no Edit toggle).
_MPT_DLG_W = 560

# Give the reused import/save/export widgets one consistent dialog width…
for _w in (mpt_json_path_input, mpt_csv_path_input, apt_path_input,
           apt_program_input, session_load_path_input,
           session_save_path_input, export_dir_input,
           mpt_plan_select, apt_plan_select,
           mpt_load_btn, mpt_csv_load_btn, apt_fetch_btn, apt_load_btn,
           session_load_btn, session_save_btn, export_btn):
    _w.width = _MPT_DLG_W
# …and show the path fields directly (the chosen path reads as text).
for _w in (mpt_json_path_input, mpt_csv_path_input, apt_path_input,
           session_load_path_input, session_save_path_input,
           export_dir_input):
    _w.visible = True


def _mpt_dlg_card_styles():
    return {
        "position": "fixed", "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)", "background": "white",
        "border": "1px solid #c0c8d6", "border-radius": "8px",
        "box-shadow": "0 12px 36px rgba(0, 30, 80, 0.32)",
        "padding": "16px 18px", "z-index": "1000",
        "max-height": "90vh", "overflow-y": "auto",
    }


def _mpt_dlg_backdrop():
    return Div(text="", width=0, height=0, visible=False, styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)", "z-index": "999",
    })


def _mpt_dlg_header(title, close_x):
    return row(
        Div(text=f"<h3>{title}</h3>", sizing_mode="stretch_width"),
        close_x,
        css_classes=["vmpt-modal-header"], styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    )


def _mpt_dlg_caption(text):
    return Div(text=f"<div style='font-size:12px; color:#5a6b85; "
                    f"margin:-2px 0 2px'>{text}</div>", width=_MPT_DLG_W)


def _mpt_dlg_section(text):
    return Div(text=f"<div class='vmpt-modal-section'>{text}</div>",
               width=_MPT_DLG_W)


# ── Import dialog (dropdown-picker chooses the source) ───────────────
import_modal_close_x = Button(label="×", button_type="default", width=32,
                              height=28, css_classes=["vmpt-modal-x"])
import_modal_close_btn = Button(label="Close", button_type="default", width=90)
import_source_select = Select(
    title="Import from", value="APT / MPT plan (JSON)",
    options=["APT / MPT plan (JSON)", "Shutter mask (CSV)",
             "APT program (.aptx or ID)", "vMPT session"],
    width=_MPT_DLG_W,
)
_imp_grp_json = column(
    _mpt_dlg_section("Plan JSON exported by APT / MPT"),
    row(mpt_json_browse_btn), mpt_json_path_input, mpt_plan_select, mpt_load_btn,
    spacing=6, width=_MPT_DLG_W,
)
_imp_grp_csv = column(
    _mpt_dlg_section("Shutter-mask CSV (open mask only)"),
    row(mpt_csv_browse_btn), mpt_csv_path_input, mpt_csv_load_btn,
    spacing=6, width=_MPT_DLG_W, visible=False,
)
_imp_grp_apt = column(
    _mpt_dlg_section("APT .aptx file on disk, or a JWST program ID from STScI"),
    row(apt_path_browse_btn), apt_path_input, apt_program_input,
    apt_fetch_btn, apt_plan_select, apt_load_btn,
    spacing=6, width=_MPT_DLG_W, visible=False,
)
_imp_grp_session = column(
    _mpt_dlg_section("Restore a saved vMPT session / workspace"),
    row(session_load_browse_btn), session_load_path_input, session_load_btn,
    spacing=6, width=_MPT_DLG_W, visible=False,
)
_IMPORT_GROUPS = {
    "APT / MPT plan (JSON)": _imp_grp_json,
    "Shutter mask (CSV)": _imp_grp_csv,
    "APT program (.aptx or ID)": _imp_grp_apt,
    "vMPT session": _imp_grp_session,
}


def _on_import_source(attr, old, new):
    for _k, _g in _IMPORT_GROUPS.items():
        _g.visible = (_k == new)


import_source_select.on_change("value", _on_import_source)
import_modal_backdrop = _mpt_dlg_backdrop()
import_modal_card = column(
    _mpt_dlg_header("Import", import_modal_close_x),
    _mpt_dlg_caption("Bring an existing plan, shutter mask, APT program, or "
                     "saved vMPT session into the canvas."),
    import_source_select,
    _imp_grp_json, _imp_grp_csv, _imp_grp_apt, _imp_grp_session,
    row(import_modal_close_btn),
    spacing=10, width=_MPT_DLG_W + 36, visible=False,
    css_classes=["vmpt-modal-card"], styles=_mpt_dlg_card_styles(),
)

# ── Save dialog ──────────────────────────────────────────────────────
save_modal_close_x = Button(label="×", button_type="default", width=32,
                            height=28, css_classes=["vmpt-modal-x"])
save_modal_close_btn = Button(label="Close", button_type="default", width=90)
save_modal_backdrop = _mpt_dlg_backdrop()
save_modal_card = column(
    _mpt_dlg_header("Save session", save_modal_close_x),
    _mpt_dlg_caption("Write a vMPT session bundle you can re-open later or "
                     "share with a collaborator."),
    _mpt_dlg_section("Save to"),
    row(session_save_browse_btn), session_save_path_input, session_save_btn,
    row(save_modal_close_btn),
    spacing=10, width=_MPT_DLG_W + 36, visible=False,
    css_classes=["vmpt-modal-card"], styles=_mpt_dlg_card_styles(),
)

# ── Export dialog ────────────────────────────────────────────────────
export_modal_close_x = Button(label="×", button_type="default", width=32,
                              height=28, css_classes=["vmpt-modal-x"])
export_modal_close_btn = Button(label="Close", button_type="default", width=90)
export_modal_backdrop = _mpt_dlg_backdrop()
export_modal_card = column(
    _mpt_dlg_header("Export to APT", export_modal_close_x),
    _mpt_dlg_caption(f"Write an eMPT bundle + {MPT_PLAN_FILENAME} + an "
                     f"APT-importable .cat for this plan."),
    _mpt_dlg_section("Output directory"),
    row(export_dir_browse_btn), export_dir_input, export_btn,
    row(export_modal_close_btn),
    spacing=10, width=_MPT_DLG_W + 36, visible=False,
    css_classes=["vmpt-modal-card"], styles=_mpt_dlg_card_styles(),
)


def _open_import_modal():
    import_modal_backdrop.visible = True
    import_modal_card.visible = True


def _close_import_modal():
    import_modal_backdrop.visible = False
    import_modal_card.visible = False


def _open_save_modal():
    save_modal_backdrop.visible = True
    save_modal_card.visible = True


def _close_save_modal():
    save_modal_backdrop.visible = False
    save_modal_card.visible = False


def _open_export_modal():
    export_modal_backdrop.visible = True
    export_modal_card.visible = True


def _close_export_modal():
    export_modal_backdrop.visible = False
    export_modal_card.visible = False


import_modal_close_x.on_click(_close_import_modal)
import_modal_close_btn.on_click(_close_import_modal)
save_modal_close_x.on_click(_close_save_modal)
save_modal_close_btn.on_click(_close_save_modal)
export_modal_close_x.on_click(_close_export_modal)
export_modal_close_btn.on_click(_close_export_modal)

# MPT-tab launcher buttons (the tab itself is just these three).
mpt_open_import_btn = Button(label="📥  Import…", button_type="primary",
                             width=SIDEBAR_W - 20, height=40)
mpt_open_save_btn = Button(label="💾  Save session…", button_type="primary",
                           width=SIDEBAR_W - 20, height=40)
mpt_open_export_btn = Button(label="📤  Export to APT…", button_type="success",
                             width=SIDEBAR_W - 20, height=40)
mpt_open_import_btn.on_click(_open_import_modal)
mpt_open_save_btn.on_click(_open_save_modal)
mpt_open_export_btn.on_click(_open_export_modal)


# ─────────────────────────────────────────────────────────────────────
# Input-tab pop-up dialogs — Load image / Load catalog
# ─────────────────────────────────────────────────────────────────────
# Same pattern as the MPT tab: the Input tab is a couple of launcher
# buttons (+ the live catalog list), and the file-picking lives in
# dialogs. Reuses the shared `_mpt_dlg_*` helpers + modal CSS.


def _tab_caption(text):
    """Small grey one-liner under a launcher button on a tab."""
    return Div(text=f"<div style='font-size:11px; color:#7a8699; "
                    f"margin:0 0 12px 2px'>{text}</div>", width=SIDEBAR_W - 20)


def _section_header(label, tip=""):
    """Uppercase, divider-topped group heading shared by the tabs. When
    `tip` is given the title shows a native hover tooltip (and a help
    cursor) so users get inline documentation for that group."""
    attr = f' title="{tip}"' if tip else ""
    cursor = "help" if tip else "default"
    return Div(
        text=f"<div{attr} style='font-size:11px; font-weight:700; "
             f"letter-spacing:.04em; text-transform:uppercase; color:#3f4d66; "
             f"cursor:{cursor}; margin:5px 0 0; border-top:1px solid #d3dae6; "
             f"padding-top:4px'>{label}</div>",
        width=SIDEBAR_W - 20, margin=(0, 0, 0, 0),
    )


# Dialog widths for the relocated pickers (Browse-only; path shown beneath).
for _w in (fits_path_input, jpg_path_input, sidecar_path_input,
           catalog_path_input):
    _w.width = _MPT_DLG_W
    _w.visible = True
for _w in (example_a370_btn, example_r0600_btn):
    _w.width = (_MPT_DLG_W - 8) // 2
catalog_add_btn.width = 140

# ── Load image dialog (dropdown picker: example / FITS / JPG+WCS) ─────
load_image_close_x = Button(label="×", button_type="default", width=32,
                            height=28, css_classes=["vmpt-modal-x"])
load_image_close_btn = Button(label="Close", button_type="default", width=90)
load_image_source_select = Select(
    title="Image source", value="Example field",
    options=["Example field", "FITS image", "JPG / PNG + WCS sidecar"],
    width=_MPT_DLG_W,
)
_img_grp_example = column(
    _mpt_dlg_section("Quick-start example fields"),
    row(example_a370_btn, example_r0600_btn, spacing=8),
    spacing=6, width=_MPT_DLG_W,
)
_img_grp_fits = column(
    _mpt_dlg_section("Local FITS image (WCS read from its header)"),
    row(fits_browse_btn), fits_path_input,
    spacing=6, width=_MPT_DLG_W, visible=False,
)
_img_grp_jpg = column(
    _mpt_dlg_section("JPG / PNG + a WCS sidecar FITS"),
    _mpt_dlg_caption("Pick the WCS sidecar first, then the image."),
    row(sidecar_browse_btn), sidecar_path_input,
    row(jpg_browse_btn), jpg_path_input,
    spacing=6, width=_MPT_DLG_W, visible=False,
)
_IMAGE_GROUPS = {
    "Example field": _img_grp_example,
    "FITS image": _img_grp_fits,
    "JPG / PNG + WCS sidecar": _img_grp_jpg,
}


def _on_image_source(attr, old, new):
    for _k, _g in _IMAGE_GROUPS.items():
        _g.visible = (_k == new)


load_image_source_select.on_change("value", _on_image_source)
load_image_modal_backdrop = _mpt_dlg_backdrop()
load_image_modal_card = column(
    _mpt_dlg_header("Load image", load_image_close_x),
    _mpt_dlg_caption("Load a background image — an example field, a FITS file, "
                     "or a JPG/PNG with a WCS sidecar."),
    load_image_source_select,
    _img_grp_example, _img_grp_fits, _img_grp_jpg,
    row(load_image_close_btn),
    spacing=10, width=_MPT_DLG_W + 36, visible=False,
    css_classes=["vmpt-modal-card"], styles=_mpt_dlg_card_styles(),
)

# ── Load catalog dialog ──────────────────────────────────────────────
load_catalog_close_x = Button(label="×", button_type="default", width=32,
                              height=28, css_classes=["vmpt-modal-x"])
load_catalog_close_btn = Button(label="Close", button_type="default", width=90)
load_catalog_modal_backdrop = _mpt_dlg_backdrop()
load_catalog_modal_card = column(
    _mpt_dlg_header("Load catalog", load_catalog_close_x),
    _mpt_dlg_caption("Add a target catalog (CSV / ASCII / FITS with ID, RA, "
                     "Dec). Add several to layer them; toggle / reorder / "
                     "remove them from the Input tab."),
    _mpt_dlg_section("Catalog file"),
    row(catalog_browse_btn), catalog_path_input, catalog_add_btn,
    row(load_catalog_close_btn),
    spacing=10, width=_MPT_DLG_W + 36, visible=False,
    css_classes=["vmpt-modal-card"], styles=_mpt_dlg_card_styles(),
)


def _open_load_image_modal():
    load_image_modal_backdrop.visible = True
    load_image_modal_card.visible = True


def _close_load_image_modal():
    load_image_modal_backdrop.visible = False
    load_image_modal_card.visible = False


def _open_load_catalog_modal():
    load_catalog_modal_backdrop.visible = True
    load_catalog_modal_card.visible = True


def _close_load_catalog_modal():
    load_catalog_modal_backdrop.visible = False
    load_catalog_modal_card.visible = False


load_image_close_x.on_click(_close_load_image_modal)
load_image_close_btn.on_click(_close_load_image_modal)
load_catalog_close_x.on_click(_close_load_catalog_modal)
load_catalog_close_btn.on_click(_close_load_catalog_modal)

load_image_open_btn = Button(label="📷  Load image…", button_type="primary",
                             width=SIDEBAR_W - 20, height=40)
load_catalog_open_btn = Button(label="🎯  Load catalog…", button_type="primary",
                               width=SIDEBAR_W - 20, height=40)
load_image_open_btn.on_click(_open_load_image_modal)
load_catalog_open_btn.on_click(_open_load_catalog_modal)


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
    text="<i>Loading vMPT… pick an example from the Input tab to begin.</i>",
)

# ── Top stats bar — order + visibility (v1.3.0+) ────────────────────
# The bar above the figure is composed of six "cells": image name,
# RA/Dec, V3 PA + APA, disperser/filter, open-shutter count, and
# spec-conflict count. Each is keyed by a short identifier; the
# corresponding HTML span is built inside `refresh_overlays_light`.
# Settings → Top stats bar exposes a MultiChoice that lets users pick
# which cells to show AND the display order — the order of the
# `value` list IS the on-screen order.
STATS_BAR_CELL_LABELS: dict[str, str] = {
    "image":     "Image filename",
    "radec":     "RA · Dec",
    "pa":        "V3 PA + APA",
    "disperser": "Disperser / Filter",
    "open":      "Open shutters",
    "conflicts": "Conflict shutters",
}
# Default order = the canonical order the bar shipped with in v1.0.
STATS_BAR_DEFAULT_ORDER: tuple[str, ...] = (
    "image", "radec", "pa", "disperser", "open", "conflicts",
)
# Reverse lookup: label → key (used when reading the picker value
# back into state). Built once.
_STATS_BAR_LABEL_TO_KEY: dict[str, str] = {
    v: k for k, v in STATS_BAR_CELL_LABELS.items()
}
# The picker widget — defined now so it's stable for layout, wired
# below in the Settings tab.
stats_bar_choice = MultiChoice(
    title="Top stats bar — pick which cells to show (drag chips or "
          "click in pick order to reorder)",
    options=[STATS_BAR_CELL_LABELS[k] for k in STATS_BAR_DEFAULT_ORDER],
    value=[STATS_BAR_CELL_LABELS[k] for k in STATS_BAR_DEFAULT_ORDER],
    width=SIDEBAR_W - 20,
)

# ── Catalog hover tooltip — order + visibility (v1.3.0+) ────────────
# The HoverTool's tooltip HTML for the catalog-target glyph is built
# from a per-field list, picked via Settings → Catalog hover. Each
# entry below maps a stable key → (display label for the picker,
# tooltip-fragment template). The template references the columns
# carried in `src_targets.data` (which `refresh_overlays` populates).
# Adding a new hover field = add an entry here + populate the column
# in src_targets.data.
CATALOG_HOVER_FIELDS: dict[str, tuple[str, str]] = {
    "id": ("ID", '<b>@id</b>'),
    "radec": (
        "RA, Dec",
        '@ra{0.0000}, @dec{0.0000}',
    ),
    "priority": ("Priority", '<span style="color:#888">Pr </span>@pr'),
    "weight": ("Weight", '<span style="color:#888">W </span>@wt'),
    "mag": ("Magnitude", '<span style="color:#888">Mag </span>@mag'),
    "z": ("Redshift", '<span style="color:#888">z </span>@z'),
    "label": ("Label", '<span style="color:#888">L </span>@label'),
    "constraints": (
        "Constraints (λ·G·B·R·🛡)",
        '<span style="color:#888">C </span>@constr',
    ),
}
CATALOG_HOVER_DEFAULT_ORDER: tuple[str, ...] = (
    "id", "radec", "priority",
)
_CATALOG_HOVER_LABEL_TO_KEY: dict[str, str] = {
    v[0]: k for k, v in CATALOG_HOVER_FIELDS.items()
}
catalog_hover_choice = MultiChoice(
    title="Catalog hover — pick which fields show in the target "
          "tooltip (drag chips or click in pick order to reorder)",
    options=[CATALOG_HOVER_FIELDS[k][0]
             for k in CATALOG_HOVER_FIELDS],
    value=[CATALOG_HOVER_FIELDS[k][0]
           for k in CATALOG_HOVER_DEFAULT_ORDER],
    width=SIDEBAR_W + 70,
)

# ── Pop-up "Customise…" dialogs (v1.3.0+) ──────────────────────────
# Both pickers used to live inline in the Settings tab. They worked,
# but the chip widget needed ~250 vertical px each which pushed the
# Actions group (Undo / Clear) off the bottom of the viewport on a
# 913 px laptop screen. Wrap each picker in its own modal dialog
# (mirrors the optimizer-config and catalog-editor patterns) so the
# Settings tab stays compact and the chip area gets the whole modal
# width for chip reflow.

# Resize the stats-bar picker to match the modal width too.
stats_bar_choice.width = SIDEBAR_W + 70

# Open / close buttons for both modals — the open buttons live in the
# Settings tab; close buttons live inside the modal cards.
stats_bar_open_btn = Button(
    label="Customise stats bar…",
    button_type="default", width=SIDEBAR_W - 20,
)
stats_bar_modal_close_btn = Button(
    label="Done", button_type="primary", width=80,
)
stats_bar_modal_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)
stats_bar_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        "z-index": "999",
    },
)

catalog_hover_open_btn = Button(
    label="Customise catalog hover…",
    button_type="default", width=SIDEBAR_W - 20,
)
catalog_hover_modal_close_btn = Button(
    label="Done", button_type="primary", width=80,
)
catalog_hover_modal_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)
catalog_hover_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.40)",
        "z-index": "999",
    },
)

# Glyph data sources
src_image = ColumnDataSource(data=dict(image=[], x=[], y=[], dw=[], dh=[]))
src_msa_outline = ColumnDataSource(data=dict(xs=[], ys=[]))
# Quadrant outlines of the OTHER (idle) MPT configs — drawn faint so the
# user sees where the other exposure sits while editing the active one.
src_idle_msa_outline = ColumnDataSource(data=dict(xs=[], ys=[], cfg=[], color=[]))
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
# Spec-overlap shutters, split into three categories that follow
# APT MPT's color encoding (v1.3.1+ change — was a single src_spec_
# overlap orange layer in v1.0–v1.3.0):
#   * stuck  → pink   (shutter's spectrum overlaps with a stuck-open
#                      shutter's dispersion only)
#   * user   → orange (overlaps with a user-open shutter's dispersion
#                      only — the "Masked" category in APT MPT)
#   * both   → purple (overlaps with BOTH stuck-open AND user-open;
#                      this is "Mask Conflict" in APT MPT)
# Per-polygon fill_alpha column lets the alpha stack with the number
# of overlapping sources: base × n, capped at 1.0.
src_spec_overlap_stuck = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], fill_alpha=[])
)
src_spec_overlap_user = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], fill_alpha=[])
)
src_spec_overlap_both = ColumnDataSource(
    data=dict(xs=[], ys=[], q=[], s=[], d=[], fill_alpha=[])
)
# Backwards-compatibility alias: the old `src_spec_overlap` symbol
# pointed at the user-overlap source. Anything still importing it
# (tests / external scripts) keeps working.
src_spec_overlap = src_spec_overlap_user
src_fixed_slits = ColumnDataSource(data=dict(xs=[], ys=[], name=[]))
src_targets = ColumnDataSource(data=dict(
    x=[], y=[], id=[], ra=[], dec=[], pr=[],
    # Extra per-target columns for the customisable hover tooltip
    # (Settings → Catalog hover, v1.3.0+). Defaults are empty so any
    # row without a value renders blank in the tooltip rather than
    # throwing.
    wt=[], mag=[], z=[], label=[], constr=[],
    # Per-target rendering. `line_color`: matched sources (id == an
    # open-shutter's target_id) flip to green; unmatched take their
    # per-catalog palette colour. `line_alpha`: decays with the
    # catalog's z-depth so earlier-loaded catalogs read more strongly
    # than later-loaded reference catalogs that sit beneath them.
    line_color=[], line_width=[], line_alpha=[],
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
    # Canvas SIZE in pixels is fixed via frame_width / frame_height —
    # these are the dimensions of the actual data-drawing area. The
    # outer figure (axes + toolbar around the frame) grows by ~70 px
    # in each direction. We update frame_width / frame_height every
    # image load (in `refresh_image_glyph`) to maintain the image's
    # pixel aspect EXACTLY — so the image renders at its native ratio
    # regardless of window size. `sizing_mode` is intentionally not
    # set; the figure stays a fixed pixel block in the layout and the
    # surrounding column letterboxes around it.
    frame_width=800,
    frame_height=800,
    match_aspect=True,
    # IMPORTANT: leave x_range / y_range at the default DataRange1d.
    # Per Bokeh docs match_aspect=True only works with DataRange1d —
    # switching to explicit Range1d silently breaks the aspect lock and
    # rectangular images (e.g. 2200×2500 FITS) get stretched horizontally
    # by ~factor canvas_aspect/data_aspect.
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
# Pin DataRange1d to track only the image extent. Otherwise overlay
# renderers whose data lies outside [0, W] × [0, H] (e.g. MSA outline
# at corner shutters when the image is smaller than the MSA, or
# catalog markers at the image edges) extend the auto-range and
# break the aspect lock — the symptom is the image stretching
# horizontally when both image and catalog are loaded together
# (e.g. via run.sh --jpg --wcs --catalog).
fig.x_range.renderers = [img_glyph]
fig.y_range.renderers = [img_glyph]

# Bottom-to-top render order: image, operable shutters, stuck-open shutters,
# open shutters (user picks), MSA outline, fixed slits,
# targets, pointing handle.
bg_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_bg_shutters,
    line_color="silver", line_alpha=0.20, line_width=1.0,
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
# User-opens render BELOW spec-overlap so contamination on a
# user-opened shutter is visible — without this z-order the red
# user-open fill would hide the pink/orange/purple stripe that
# flags "your own pick is contaminated by another open shutter".
open_shutters_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_open_shutters,
    line_color="#ff3333", line_width=1.5,
    fill_color="#ff8888", fill_alpha=0.35,
)
# Spec-overlap shutters in three MPT-faithful colors (v1.3.1+).
# Each glyph has `fill_alpha="fill_alpha"` so the source's per-
# polygon alpha column drives transparency — alpha stacks with the
# number of overlapping dispersion sources (base × n, capped at 1).
# Edge defaults to a thin 0.5 px outline (line_alpha 0.6) so the
# masked shutters read clearly even at low fill alpha; the
# appearance picker's stroke slider tunes it (0 hides the edge).
#
# Rendered AFTER `open_shutters_glyph` so contamination on a
# user-opened shutter (open-vs-open and open-vs-stuck conflicts)
# overlays the red fill — making mutual conflicts among picked
# shutters visible. The user-open's red still shows through under
# the spec-overlap's alpha.
spec_overlap_stuck_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap_stuck,
    line_color="#d96272", line_alpha=0.6, line_width=0.5,
    # Pink (APT MPT "Mask Stuck" colour): a shutter whose
    # spectrum overlaps with a stuck-open's dispersion.
    fill_color="#ff9aa0", fill_alpha="fill_alpha",
)
spec_overlap_user_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap_user,
    line_color="#d97a00", line_alpha=0.6, line_width=0.5,
    # Orange (APT MPT "Masked" colour): a shutter whose spectrum
    # overlaps with a user-opened shutter's dispersion.
    fill_color="orange", fill_alpha="fill_alpha",
)
spec_overlap_both_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap_both,
    line_color="#6f3a8a", line_alpha=0.6, line_width=0.5,
    # Purple (APT MPT "Mask Conflict" colour): a shutter whose
    # spectrum overlaps with BOTH stuck-open AND user-open
    # dispersion sources.
    fill_color="#a050b8", fill_alpha="fill_alpha",
)
# Backwards-compatibility alias: the old `spec_overlap_glyph` symbol
# pointed at the orange (now: user-overlap) glyph. The overlay-
# appearance config below still references it under the old key.
spec_overlap_glyph = spec_overlap_user_glyph
# Idle configs' outlines render UNDER the active outline (drawn first):
# faint dashed blue with a whisper of fill so they read as "the other
# config's footprint" without competing with the active solid-blue grid.
idle_msa_outline_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_idle_msa_outline,
    line_color="color", line_width=1.0, line_dash="dashed",
    line_alpha=0.7, fill_alpha=0.0,  # boundary only — no fill in the box
)
# Active config's quadrant outline. `line_color` is repainted per active
# config in `refresh_overlays` (Config 1 blue, Config 2 magenta, …) to
# match the top-bar CONFIG chip; the initial value is the Config-1 colour.
msa_outline_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_msa_outline,
    line_color=_config_color(0), line_width=1.5, fill_alpha=0.0,
)
fixed_slits_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_fixed_slits,
    line_color="gold", line_width=2.0, fill_color="gold", fill_alpha=0.18,
)
target_glyph = fig.scatter(
    x="x", y="y", source=src_targets,
    size=10, marker="circle",
    line_color="line_color",   # field-driven: per-catalog colour, green when matched
    line_width="line_width",
    line_alpha="line_alpha",   # field-driven: decays with catalog z-depth
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
# Catalog hover tool — tooltip content is rebuilt at runtime by
# `_refresh_catalog_hover_tooltip` whenever the user reorders / toggles
# the Settings → Catalog hover picker. Field list is small enough that
# we just store the HoverTool on the module so the callback can mutate
# its `tooltips` attribute in place.
catalog_hover = HoverTool(
    renderers=[target_glyph],
    tooltips=(
        f'<div style="{_TIP_BASE_STYLE} color:#1a3b66;">'
        f'  <b>@id</b>'
        f'  <span style="color:#888;"> · </span>@ra{{0.0000}}, @dec{{0.0000}}'
        f'  <span style="color:#888;"> · Pr </span>@pr'
        f'</div>'
    ),
)
fig.add_tools(catalog_hover)

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
    """Trace each quadrant outline along the OUTER corners of its 4 corner
    shutters, so the blue outline hugs the rendered shutter area.

    Previously the outline ran through the corner shutter *centres*, which
    left the outermost half-shutter (~0.23″, ~3 screen px) sticking out
    past the blue line. For each quadrant-corner shutter we pick the one of
    its 4 footprint corners that is farthest from the quadrant's centre, so
    the polygon encloses every shutter regardless of the MSA's V2/V3
    rotation/shear.
    """
    xs_all, ys_all = [], []
    rows = [0, 0, 170, 170]
    cols = [0, 364, 364, 0]
    tdv2 = _SHUTTER_CORNER_TEMPLATE[:, 0]
    tdv3 = _SHUTTER_CORNER_TEMPLATE[:, 1]
    for q in range(4):
        cen_v2 = float(V2_MSA[q].mean())
        cen_v3 = float(V3_MSA[q].mean())
        v2_out, v3_out = [], []
        for r, c in zip(rows, cols):
            v2c = float(V2_MSA[q, r, c])
            v3c = float(V3_MSA[q, r, c])
            cand_v2 = v2c + tdv2
            cand_v3 = v3c + tdv3
            k = int(np.argmax((cand_v2 - cen_v2) ** 2
                              + (cand_v3 - cen_v3) ** 2))
            v2_out.append(cand_v2[k])
            v3_out.append(cand_v3[k])
        x, y = _project_v2v3_offsets_to_pixel(
            np.array(v2_out) - MSA_V2_REF, np.array(v3_out) - MSA_V3_REF,
            pa_v3, fid_pix, jinv,
        )
        xs_all.append([[x.tolist()]])
        ys_all.append([[y.tolist()]])
    return dict(xs=xs_all, ys=ys_all)


def _idle_config_outlines_data(wcs) -> dict:
    """Quadrant outlines for every non-active live config that sits at its
    own distinct pointing, each projected at THAT config's pointing.

    Returns empty when single-config, or for configs that have never been
    positioned (their pointing is None → they share the active pointing, so
    drawing them would just overlap the active grid)."""
    n = int(state.get("n_configs", 1))
    if n <= 1:
        return dict(xs=[], ys=[], cfg=[], color=[])
    active = int(state.get("active_config", 0))
    a_ra = state.get("ra_deg")
    a_dec = state.get("dec_deg")
    a_pa = state.get("pa_v3")
    xs_all, ys_all, cfg_all, color_all = [], [], [], []
    for ci in range(min(n, len(state["configs"]))):
        if ci == active:
            continue
        cfg = state["configs"][ci]
        ra_c, dec_c, pa_c = cfg.get("ra_deg"), cfg.get("dec_deg"), cfg.get("pa_v3")
        if ra_c is None or dec_c is None or pa_c is None:
            continue
        try:
            # Skip if it coincides with the active pointing (would overlap).
            if (a_ra is not None and abs(float(ra_c) - float(a_ra)) < 1e-9
                    and abs(float(dec_c) - float(a_dec)) < 1e-9
                    and abs(float(pa_c) - float(a_pa)) < 1e-9):
                continue
            sc = SkyCoord(float(ra_c), float(dec_c), unit=u.deg, frame="icrs")
            fx, fy = _world_to_pixel(sc, wcs)
            fp = (float(fx), float(fy))
            jv = _compute_wcs_jacobian(wcs, fp[0], fp[1])
            d = _msa_outline_polygons(float(pa_c), fp, jv)
        except Exception:  # noqa: BLE001
            continue
        for poly_xs, poly_ys in zip(d["xs"], d["ys"]):
            xs_all.append(poly_xs)
            ys_all.append(poly_ys)
            cfg_all.append(ci + 1)
            color_all.append(_config_color(ci))
    return dict(xs=xs_all, ys=ys_all, cfg=cfg_all, color=color_all)


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
        cut = cutoffs(
            v2c, v3c, state["disperser"], state["filter"],
            q=q, s=s, d=d,
        )
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


def _grating_adjacency_relaxed(lo1, hi1, d1, lo2, hi2, d2, min_colsep) -> bool:
    """Whether a SAME-quadrant conflict between two slitlets should be demoted
    from purple (Mask Conflict) to orange (Masked).

    True only for a pure no-buffer ADJACENCY — the two slitlets' row ranges
    are disjoint with EXACTLY 0 rows between them (`b_lo = a_hi + 1`) — under a
    grating (``min_colsep`` is the disperser's column threshold, not None) with
    the two slitlets at least ``min_colsep`` columns apart. Long grating
    spectra run parallel, so such a column-offset diagonal step isn't a real
    collision. Real row overlap (rows intersect) and PRISM (``min_colsep`` is
    None) always return False → they stay purple.
    """
    if min_colsep is None:
        return False
    if hi1 < lo2:
        gap = lo2 - hi1 - 1          # slitlet 1 is the lower one
    elif hi2 < lo1:
        gap = lo1 - hi2 - 1          # slitlet 2 is the lower one
    else:
        return False                 # rows overlap → real conflict, stay purple
    return gap == 0 and abs(d1 - d2) >= min_colsep


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

    # MSA outline (active config). Repaint its line colour to the active
    # config's accent so the solid boundary matches the top-bar chip
    # (Config 1 blue, Config 2 magenta, …).
    show_outline = 0 in layers_box.active
    msa_outline_glyph.glyph.line_color = _config_color(
        int(state.get("active_config", 0)))
    if show_outline:
        src_msa_outline.data = _msa_outline_polygons(pa_v3, fid_pix, jinv)
    else:
        src_msa_outline.data = dict(xs=[], ys=[])
    # Faint outlines of the OTHER configs (multi-config mode), each at its
    # own pointing and in its own accent colour, so the user sees where the
    # idle exposure(s) sit and which is which.
    if show_outline and int(state.get("n_configs", 1)) > 1:
        src_idle_msa_outline.data = _idle_config_outlines_data(wcs)
    else:
        src_idle_msa_outline.data = dict(xs=[], ys=[], cfg=[], color=[])

    # Cull to the current visible figure bounds (post-zoom). Use the ACTUAL
    # range — NOT clamped to the image extent — so shutters that project
    # outside the image cutout still render whenever that region is on
    # screen (the MSA can extend past a small image, and the image is often
    # letterboxed inside the canvas, leaving operable area beyond [0,W]×
    # [0,H]). `on_ranges_update` re-runs refresh_overlays on every pan/zoom,
    # so off-image shutters appear as soon as the user scrolls to them. All
    # three shutter-flavour layers (operable, stuck-open, spectral-overlap)
    # reuse the one mask.
    try:
        vx0, vx1 = float(fig.x_range.start), float(fig.x_range.end)
        vy0, vy1 = float(fig.y_range.start), float(fig.y_range.end)
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
    # Group user-opens into slitlets. An N-shutter slitlet opened on
    # one source (same target_id, same (q, d), N consecutive s rows)
    # disperses one continuous spectrum onto the detector — not N
    # independent ones. Without this grouping a 3-shutter slitlet
    # would increment its candidates' conflict count by 3 instead of
    # 1, tripling the alpha and producing a darker stripe than MPT
    # shows.
    #
    # Group key is `(q, d, target_id_or_anon)`: shutters at the same column
    # with the same target id form one slitlet. BUT a slitlet is a
    # CONTIGUOUS column of shutters, so each bucket is then split into
    # maximal runs of consecutive rows. Opens loaded from a mask CSV carry
    # no target id, so two physically separate slitlets in the same column
    # (e.g. Q4 d=303 with s=22-24 AND s=108-110) would otherwise merge into
    # one bogus group spanning s=22-110 — whose (min,max) row span then
    # (a) masks the whole empty gap between the clusters and (b) never
    # shrinks when you close a single shutter, so the Mask Conflict never
    # clears. Splitting on row gaps fixes both. A 1-row gap (a single
    # failed-closed shutter inside a real slitlet) keeps the run together.
    _raw_groups: dict[tuple, list[tuple[int, int, int]]] = {}
    for (q, s, d), sh in state["open_shutters"].items():
        tid = getattr(sh, "target_id", None)
        key = (q, d, tid if tid is not None else "_anon_")
        _raw_groups.setdefault(key, []).append((q, s, d))
    user_groups: dict[tuple, list[tuple[int, int, int]]] = {}
    for key, shts in _raw_groups.items():
        shts.sort(key=lambda t: t[1])            # by MSA row s
        run = [shts[0]]
        run_i = 0
        for prev, cur in zip(shts, shts[1:]):
            if cur[1] - prev[1] <= 1:            # adjacent row → same slitlet
                run.append(cur)
            else:                                # gap → separate slitlet
                user_groups[(*key, run_i)] = run
                run_i += 1
                run = [cur]
        user_groups[(*key, run_i)] = run

    stuck_flat = np.where(_FLAT_REASON == 2)[0]
    stuck_keys = [
        (int(f // (171 * 365)) + 1,
         int((f % (171 * 365)) // 365) + 1,
         int(f % 365) + 1)
        for f in stuck_flat
    ]
    # Stuck shutters are individual physical defects, not slitlets —
    # each one is its own dispersion source.
    stuck_groups: list[list[tuple[int, int, int]]] = [[k] for k in stuck_keys]

    # Per-affected-shutter conflict counts, split FOUR ways:
    #   user_direct, user_buffer, stuck_direct, stuck_buffer.
    #
    # DIRECT means the candidate row lies in the slitlet's actual
    # row range (i.e. the spectrum's true cross-dispersion path
    # lands on this row). BUFFER means the candidate is at the ±1
    # tolerance edge — a conservative safety margin, NOT an actual
    # contamination. Only direct hits count as "real" overlap; the
    # buffer rows are warnings for diagnostic display.
    #
    # The total counts (user_counts, stuck_counts) sum direct + buffer
    # and are used by the silver-edge filter and the stats bar.
    user_direct: dict[int, int] = {}
    user_buffer: dict[int, int] = {}
    stuck_direct: dict[int, int] = {}
    stuck_buffer: dict[int, int] = {}
    stuck_counts: dict[int, int] = {}
    user_counts: dict[int, int] = {}
    # Picks of slitlets in a cross-quadrant conflict — forced purple even when
    # not doubly-hit (built on a cache miss; empty otherwise).
    xq_conflict_pick_flats: set[int] = set()

    # Base alpha per category — slider-controlled in Settings →
    # Overlay appearance, persisted via the prefs system, stashed
    # on `state`. Per-polygon alpha = min(1, base × n_conflicts)
    # where n_conflicts is the total dispersing-source count
    # overlapping that shutter. Read here (before the accumulation)
    # because the spec-overlap cache signature includes them.
    base_alpha_stuck = float(state.get("overlap_base_alpha_stuck", 0.20))
    base_alpha_user = float(state.get("overlap_base_alpha_user", 0.20))
    base_alpha_both = float(state.get("overlap_base_alpha_both", 0.20))

    # ── Spec-overlap memoization ──────────────────────────────────────
    # The conflict identification + colour assignment below is purely a
    # function of the open-shutter mask, the disperser/filter, and the
    # alpha sliders — it does NOT depend on the pointing, the V3 PA, or
    # the zoom (those only enter at *projection* time, in the `_iv` /
    # `_partition_and_project` calls further down). But refresh_overlays
    # re-runs on every pan/zoom (RangesUpdate), and for dense masks the
    # global accumulation + colouring is hundreds of ms. So we memoise
    # the four view-independent outputs (the three alpha dicts + the
    # affected-index set) keyed by that signature, and on a pure pan/zoom
    # reuse them — only the cheap render-time cull (`_iv`) re-runs.
    _open_sig = frozenset(
        (q, s, d, getattr(sh, "target_id", None))
        for (q, s, d), sh in state["open_shutters"].items()
    )
    _overlap_sig = (
        _open_sig, state["disperser"], state["filter"],
        base_alpha_stuck, base_alpha_user, base_alpha_both,
    )
    _overlap_cache = state.get("_spec_overlap_cache")
    _overlap_cached = (
        _overlap_cache is not None and _overlap_cache[0] == _overlap_sig
    )

    if (not _overlap_cached) and (user_groups or stuck_keys):
        v2_overlap = float(v2_overlap_distance(state["disperser"], state["filter"]))
        s_arr = (np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) % (171 * 365)) // 365
        q_arr = np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) // (171 * 365) + 1

        # Per-shutter on-detector x/y range (loaded lazily by combo).
        # Two complementary uses:
        #
        #   1. SUBTRACTIVE x-range filter — drops candidates whose
        #      spectrum doesn't share any detector x with the open
        #      shutter (e.g. G140M/F070LP at ΔV2 ≈ 84″: spectra are
        #      only ~60″ wide on detector, so the V2-distance check's
        #      "potentially overlapping" verdict is wrong).
        #
        #   2. ADDITIVE cross-quadrant detector-y check — catches
        #      cross-quadrant pairs whose spectrum y-stripes coincide
        #      at the x-overlap region even though MSA s differs
        #      (e.g. G140M/F100LP Q4 s=34 ↔ Q2 s=33: MSA-row + tilt
        #      predicts a row offset that excludes the candidate, but
        #      detector y at the x-overlap is essentially identical).
        #      Same-quadrant pairs keep their existing MSA-row + tilt
        #      behaviour (so the tilt-jump-at-open-shutter visual is
        #      preserved).
        #
        # If the precomputed xy grids aren't shipped for this combo
        # the filter passes everything through unchanged (legacy V2
        # distance + MSA-row behaviour).
        from .wavelengths import shutter_xy_grids as _xy_grids
        from .wavelengths import tilt_slope_grid as _slope_grid
        _xy = _xy_grids(state["disperser"], state["filter"])
        if _xy is not None:
            _xlo_n1 = _xy["x_lo_nrs1"].reshape(-1)
            _xhi_n1 = _xy["x_hi_nrs1"].reshape(-1)
            _y_n1   = _xy["y_nrs1"  ].reshape(-1)
            _xlo_n2 = _xy["x_lo_nrs2"].reshape(-1)
            _xhi_n2 = _xy["x_hi_nrs2"].reshape(-1)
            _y_n2   = _xy["y_nrs2"  ].reshape(-1)
            _on_n1 = np.isfinite(_xlo_n1) & np.isfinite(_xhi_n1) & np.isfinite(_y_n1)
            _on_n2 = np.isfinite(_xlo_n2) & np.isfinite(_xhi_n2) & np.isfinite(_y_n2)
            _sg = _slope_grid(state["disperser"], state["filter"])
            # Convert slope (rows/arcsec) → detector px / detector px
            # using the NIRSpec geometry constants 5 px/row (cross-
            # dispersion) and ≈ 13 px/arcsec (dispersion).
            _slope_det_flat = (
                _sg.reshape(-1).astype(float) * (5.0 / 13.0)
                if _sg is not None else np.zeros_like(_xlo_n1)
            )
        else:
            _xlo_n1 = _xhi_n1 = _y_n1 = None
            _xlo_n2 = _xhi_n2 = _y_n2 = None
            _on_n1 = _on_n2 = None
            _slope_det_flat = None

        # Tilt-aware row check (MPT-faithful). NIRSpec spectral traces
        # drift slightly in the cross-dispersion direction as you move
        # along the dispersion direction — but in APT MPT the drift
        # snaps in discrete row steps, with a clear breakpoint at the
        # dispersing shutter: the spectrum stays on s_o until the
        # cumulative drift |slope·Δv2| exceeds 0.5, then jumps to
        # s_o ± 1, and so on.
        #
        # We implement that by rounding the predicted row offset to
        # the nearest integer BEFORE the row check. For PRISM the
        # rounded offset is 0 everywhere (sub-row tilt across the full
        # V2 extent), so PRISM matches the flat-row v1.0–v1.3.0
        # behaviour exactly. For the M and H gratings the discrete
        # 1-row steps become visible at the spectrum's V2 edges.
        #
        # The tilt slope (rows per arcsec V2) is bilinearly
        # interpolated from a 10×10-per-quadrant grid precomputed via
        # jwst.assign_wcs (see scripts/precompute_trace_tilt.py).
        # Missing combo → slope = 0 → identical to flat-row.
        from .wavelengths import tilt_slope_for_shutter

        v2_all = _V2_OFFSETS_ALL + MSA_V2_REF
        # 0-based column index per flat element (length 4*171*365).
        d_arr = np.arange(_V2_OFFSETS_ALL.size, dtype=np.int64) % 365

        # `hit_sources[i]` = the set of (source_type, source_idx)
        # tuples that contribute a direct or buffer hit to candidate
        # `i`. Used after both _accumulate passes to identify which
        # source slitlets are in active collision (their own shutters
        # have been hit by ANOTHER source) — those slitlets' entire
        # contamination bands then propagate purple.
        hit_sources: dict[int, set[tuple[str, int]]] = {}

        # Broadened candidate mask: include stuck-open shutters too so
        # we can detect user-pick ↔ stuck-open touching collisions. The
        # silver-edge layer computes its own operable-AND-in-view mask
        # further down; this `in_view_candidates` is just for the
        # contamination accumulator.
        #
        # IMPORTANT — the contamination computation is **view-INDEPENDENT**.
        # A shutter's overlap colour (purple / orange / pink) is a property
        # of the open mask + disperser + pointing, NOT of where the user is
        # currently looking. Culling the candidate set to the visible view
        # made the conflict identification — and therefore the colours —
        # change with zoom/pan (a conflicting partner scrolling off-screen
        # would silently demote its whole band purple→orange→pink). So we
        # accumulate over ALL operable + stuck shutters across the MSA here;
        # the result is culled to the view only at render time (the
        # `_partition_and_project` calls below), keeping the polygon count
        # bounded without making the colours depend on the viewport.
        in_view_candidates = _FLAT_REASON != 1

        def _accumulate(
            source_type: str,
            groups,
            dst_direct: dict[int, int],
            dst_buffer: dict[int, int],
        ) -> None:
            """Each `group` is a list of (q, s, d) shutters in one
            dispersing slitlet (same column, consecutive rows). Records
            per-candidate hits in TWO bins:

            * ``dst_direct`` — candidate row lies in the slitlet's
              actual row range (with row-offset tilt), i.e. the
              spectrum's true cross-dispersion path lands on this row.
              These are "real" contamination hits.

            * ``dst_buffer`` — candidate row is at the ±1 tolerance
              edge around the direct range (the safety buffer for the
              ±1-row collision conservatism). The spectrum's TRUE row
              doesn't actually hit this candidate; the buffer is for
              edge cases (tilt fluctuation, mis-calibration).

            The MPT-faithful purple ("Mask Conflict") classification
            requires DIRECT overlap — buffer-only overlaps from two
            sources at adjacent rows (e.g. stuck-opens at s=155 and
            s=157 both ±1 to row 156) must not be treated as a real
            collision.
            """
            for src_idx, group in enumerate(groups):
                if not group:
                    continue
                qs = [g[0] for g in group]
                ss = [g[1] for g in group]
                ds = [g[2] for g in group]
                q_o = qs[0]
                d_o = ds[0]
                s_center = int(np.round(np.mean(ss)))
                slope_k = tilt_slope_for_shutter(
                    state["disperser"], state["filter"],
                    q_o, s_center, d_o,
                )
                # Same detector half = both quadrants that fold onto this
                # detector (Q1/Q3 → NRS1, Q2/Q4 → NRS2). The spectrograph
                # optics map the partner quadrants onto the SAME detector
                # rows: Q1 s=k and Q3 s=k land within ~3 px in detector-y
                # (less than one 5-px shutter) — verified against the
                # `y_nrs*` grids — so matching the partner quadrant by `s`
                # here is correct (a real, same-row collision), NOT over-
                # firing. Do NOT restrict this to `q_arr == q_o`: the large
                # V2 separation between the quadrants is irrelevant because
                # `s` already encodes the (shared) detector row, and the
                # restriction would silently drop the real Q1↔Q3 / Q2↔Q4
                # collisions.
                partners = NRS1_QUADS if q_o in NRS1_QUADS else NRS2_QUADS
                same_det = np.isin(q_arr, list(partners))
                anchor_flat = (
                    (q_o - 1) * 171 * 365 + (s_center - 1) * 365 + (d_o - 1)
                )
                if not (0 <= anchor_flat < v2_all.size):
                    continue
                v2_o = float(v2_all[anchor_flat])
                v2_open_row = V2_MSA[q_arr - 1, s_center - 1, d_arr]
                dv2 = v2_open_row - v2_o
                # Tilt-aware row prediction was tried (round-half-away-
                # from-zero of slope×Δv2) but the resulting band shifts
                # by ±1 row when |drift| crosses 0.5 — at far-d edges
                # of v2_overlap the band's union becomes N+3 or N+4 rows
                # instead of the expected N+2. Users found the diagonal
                # jump visually distracting (see G140M Q1 s=115 d=37
                # screenshot — band sticks out at d ≈ 153 transition).
                # We clamp row_offset = 0 so the band stays exactly at
                # s_o ± (half_extent + SHVAL_S_TOLERANCE) across the
                # full v2_overlap range. Trade-off: at far-d candidates
                # where the actual spectrum drift exceeds 0.5 rows, the
                # rendered band doesn't follow the spectrum — but
                # within typical M/H grating slopes and v2_overlap
                # ranges, the dropped contamination at the wing edges
                # is small.
                row_offset = np.zeros_like(q_arr, dtype=np.int64)
                # `drift` retained so any downstream debug tools can
                # still read it; not used in the row check.
                _drift_for_debug = slope_k * dv2  # noqa: F841
                different_col = d_arr != (d_o - 1)
                near_v2 = different_col & (np.abs(dv2) < v2_overlap)
                # Cross-dispersion footprint = the slitlet's ACTUAL open rows
                # [s_lo, s_hi] ± the ±1 buffer — NOT s_center ± half_extent.
                # The centre-based form is symmetric around a *rounded* centre,
                # which skews even-shutter slitlets: a 2-shutter slitlet at
                # s=[112,113] rounds its centre to 112 and so reaches one row
                # LOWER than the 3-shutter slitlet [112,113,114] that contains
                # it — letting a subset flag a collision its own superset does
                # not. Anchoring on [s_lo, s_hi] makes the band monotonic in
                # the open set (superset band ⊇ subset band). For odd slitlets
                # it is identical to the old centre form. row_offset (tilt) = 0.
                s_lo0 = (min(ss) - 1) + row_offset   # lowest open row (0-based)
                s_hi0 = (max(ss) - 1) + row_offset   # highest open row (0-based)
                direct_row = (s_arr >= s_lo0) & (s_arr <= s_hi0)
                buffer_row = (
                    (s_arr >= s_lo0 - SHVAL_S_TOLERANCE)
                    & (s_arr <= s_hi0 + SHVAL_S_TOLERANCE)
                    & (~direct_row)
                )

                # === SUBTRACTIVE: x-range filter ===
                # Drops candidates whose spectrum doesn't physically
                # share any detector x with the open (e.g. F070LP
                # spectra at ΔV2≈84″).
                if _xlo_n1 is not None:
                    o_on_n1f = bool(_on_n1[anchor_flat])
                    o_on_n2f = bool(_on_n2[anchor_flat])
                    x_share_n1 = np.zeros_like(q_arr, dtype=bool)
                    x_share_n2 = np.zeros_like(q_arr, dtype=bool)
                    if o_on_n1f:
                        o_xlo, o_xhi = float(_xlo_n1[anchor_flat]), float(_xhi_n1[anchor_flat])
                        x_share_n1 = _on_n1 & (_xlo_n1 <= o_xhi) & (_xhi_n1 >= o_xlo)
                    if o_on_n2f:
                        o_xlo, o_xhi = float(_xlo_n2[anchor_flat]), float(_xhi_n2[anchor_flat])
                        x_share_n2 = _on_n2 & (_xlo_n2 <= o_xhi) & (_xhi_n2 >= o_xlo)
                    x_share = x_share_n1 | x_share_n2
                    if not (o_on_n1f or o_on_n2f):
                        x_share[:] = False
                else:
                    x_share = np.ones_like(q_arr, dtype=bool)

                # === ADDITIVE: cross-quadrant detector-y check ===
                # MSA s is local to each quadrant — candidates in
                # other quadrants can have very different s yet land
                # on essentially the same detector y as the open's
                # spectrum at the x-overlap region (tilt brings the
                # spectra together). The MSA-row + tilt check above
                # would miss these. Compute local y at the x-overlap
                # midpoint for each shutter using its own median y +
                # slope, and OR in any candidate whose stripe meets
                # the open's stripe in y at that x.
                cross_q_direct = np.zeros_like(q_arr, dtype=bool)
                cross_q_buffer = np.zeros_like(q_arr, dtype=bool)
                if _xlo_n1 is not None and _slope_det_flat is not None:
                    different_quad = q_arr != q_o
                    SLIT_HALF_PX = 2.5   # → direct  (stripes co-aligned)
                    SLIT_FULL_PX = 5.0   # → buffer  (stripes touching)
                    # Minimum on-detector x-overlap width before a
                    # cross-quadrant hit is reported. Filters out
                    # cases where the spectra just clip each other's
                    # edge by a few pixels (e.g. G235M Q4 d=335 s=34
                    # ↔ Q2 d=280 s=33: only ~3 px x-overlap and
                    # only ~2 nm of shared spectral content —
                    # truly negligible).
                    MIN_X_OVERLAP_PX = 10.0
                    o_sd = float(_slope_det_flat[anchor_flat])
                    for det_share, xlo_a, xhi_a, y_a, on_a in (
                        (x_share_n1, _xlo_n1, _xhi_n1, _y_n1, _on_n1),
                        (x_share_n2, _xlo_n2, _xhi_n2, _y_n2, _on_n2),
                    ):
                        if xlo_a is None:
                            continue
                        o_on = bool(on_a[anchor_flat])
                        if not o_on:
                            continue
                        o_xlo_f = float(xlo_a[anchor_flat])
                        o_xhi_f = float(xhi_a[anchor_flat])
                        o_y_f   = float(y_a[anchor_flat])
                        o_xc    = 0.5 * (o_xlo_f + o_xhi_f)
                        x_int_lo = np.maximum(xlo_a, o_xlo_f)
                        x_int_hi = np.minimum(xhi_a, o_xhi_f)
                        x_mid    = 0.5 * (x_int_lo + x_int_hi)
                        x_overlap_w = np.maximum(x_int_hi - x_int_lo, 0.0)
                        y_o_local = o_y_f + o_sd * (x_mid - o_xc)
                        c_xc      = 0.5 * (xlo_a + xhi_a)
                        y_c_local = y_a + _slope_det_flat * (x_mid - c_xc)
                        dy = np.abs(y_o_local - y_c_local)
                        valid = (
                            different_quad
                            & on_a & det_share
                            & (x_overlap_w >= MIN_X_OVERLAP_PX)
                        )
                        cross_q_direct |= valid & (dy <= SLIT_HALF_PX)
                        cross_q_buffer |= valid & (dy > SLIT_HALF_PX) & (dy <= SLIT_FULL_PX)
                    cross_q_buffer &= ~cross_q_direct
                cross_q_direct &= different_col
                cross_q_buffer &= different_col

                # Combine MSA-row + tilt with cross-quadrant additions
                final_direct = (same_det & direct_row & near_v2 & x_share) | cross_q_direct
                final_buffer = (
                    ((same_det & buffer_row & near_v2 & x_share) | cross_q_buffer)
                    & ~final_direct
                )

                # Candidate mask = `in_view_candidates` (operable AND
                # stuck-open) so we can also detect when a source's
                # spectrum touches a stuck-open shutter.
                idx_direct = np.where(in_view_candidates & final_direct)[0]
                idx_buffer = np.where(in_view_candidates & final_buffer)[0]
                slitlet_flat = {
                    (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
                    for (q, s, d) in group
                }
                source_tag = (source_type, int(src_idx))
                for i in idx_direct.tolist():
                    if i in slitlet_flat:
                        continue
                    dst_direct[i] = dst_direct.get(i, 0) + 1
                    hit_sources.setdefault(i, set()).add(source_tag)
                for i in idx_buffer.tolist():
                    if i in slitlet_flat:
                        continue
                    dst_buffer[i] = dst_buffer.get(i, 0) + 1
                    hit_sources.setdefault(i, set()).add(source_tag)

        user_groups_list = list(user_groups.values())
        _accumulate("user", user_groups_list, user_direct, user_buffer)
        _accumulate("stuck", stuck_groups, stuck_direct, stuck_buffer)

        # Identify conflicting source PAIRS: any open whose own shutters
        # are hit by a DIFFERENT source. We record the actual PAIRS (not
        # just a flat "conflicted" set) so each conflict's purple band can
        # be bounded to the rows where the two slitlets actually crowd
        # each other, rather than promoting a conflicted slitlet's WHOLE
        # contamination band to purple.
        conflicted_user: set[int] = set()
        conflicted_stuck: set[int] = set()
        conflict_pairs: set = set()

        def _scan_conflicts(groups, my_type):
            for src_idx, group in enumerate(groups):
                me = (my_type, src_idx)
                others_all: set = set()
                for (q, s, d) in group:
                    flat = (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
                    others_all |= hit_sources.get(flat, set()) - {me}
                for o in others_all:
                    (conflicted_user if my_type == "user"
                     else conflicted_stuck).add(src_idx)
                    (conflicted_user if o[0] == "user"
                     else conflicted_stuck).add(o[1])
                    conflict_pairs.add(frozenset((me, o)))

        _scan_conflicts(user_groups_list, "user")
        _scan_conflicts(stuck_groups, "stuck")

        # Per-conflicted-slitlet PURPLE WINDOW. A conflict promotes shutters
        # to purple (Mask Conflict) ONLY within ±2 rows of where the two
        # slitlets crowd each other — [upper-slitlet bottom − 2 …
        # lower-slitlet top + 2] in MSA rows — NOT along the whole band.
        # (Confirmed against APT MPT: two adjacent N=3 slitlets give a
        # 2-orange / 4-purple / 2-orange stack.) Beyond the window the band
        # stays orange/pink. Because we only ever re-bucket shutters that
        # are ALREADY contaminated (below), purple ⊆ orange∪pink — a clean
        # silver shutter is never promoted, and a lone stuck-open's ±1 pink
        # band can't be inflated past where contamination actually exists.
        #
        # The window lives in MSA rows, so it only applies to a SAME-quadrant
        # pair (shared row frame). A cross-quadrant conflict (spectra fold
        # together from different quadrants — no "buffer row" concept between
        # them) keeps the full-band promotion (window = None = unbounded).
        CONFLICT_PURPLE_BUFFER = 2
        # Grating diagonal-step relaxation: under a grating, demote a
        # no-buffer ADJACENCY conflict (rows exactly 1 apart, no real overlap)
        # to orange once the two slitlets are ≥ this many columns apart. That
        # threshold is 1 (any nonzero column offset → a deliberate diagonal
        # step → orange; only exact same-column stacking Δd=0 stays purple) —
        # see wavelengths.grating_adjacency_min_colsep for why magnitude is
        # irrelevant. None for PRISM / non-gratings → adjacency always purple.
        _adj_min_colsep = grating_adjacency_min_colsep(state.get("disperser"))

        def _slit_rows(tag):
            t, si = tag
            grp = user_groups_list[si] if t == "user" else stuck_groups[si]
            ss0 = [s for (_q, s, _d) in grp]
            return grp[0][0], min(ss0), max(ss0), grp[0][2]

        slitlet_window: dict[tuple, set] = {}    # same-quad ±2 row windows
        slitlet_quad: dict[tuple, int] = {}
        # Cross-quadrant conflicting pairs share no MSA-row frame (the row
        # index isn't comparable across quadrants), so a row window is
        # meaningless. Instead remember each tag's cross-quad partners; a
        # candidate is then cross-quad PURPLE only where it is contaminated
        # by BOTH a tag and one of its partners — the genuine doubly-masked
        # overlap of the two spectra — never the tag's whole band. (The old
        # behaviour promoted each cross-quad-conflicted slitlet's ENTIRE
        # contamination band, in BOTH quadrants, to purple — tens of
        # thousands of shutters for a full mask, and impossible to clear by
        # closing same-quadrant shutters. This keeps purple ⊆ orange/pink.)
        cross_partners: dict[tuple, set] = {}
        for pair in conflict_pairs:
            tag1, tag2 = tuple(pair)
            q1, lo1, hi1, d1 = _slit_rows(tag1)
            q2, lo2, hi2, d2 = _slit_rows(tag2)
            slitlet_quad[tag1] = q1
            slitlet_quad[tag2] = q2
            if q1 != q2:
                cross_partners.setdefault(tag1, set()).add(tag2)
                cross_partners.setdefault(tag2, set()).add(tag1)
                continue
            # GRATING diagonal-step relaxation (same-quadrant only): a pure
            # ADJACENCY conflict (rows disjoint, 0 buffer rows) under a grating
            # with a large enough COLUMN offset isn't a real collision — skip
            # its purple window so those shutters stay orange (Masked). Real
            # row overlap and PRISM are unaffected.
            if _grating_adjacency_relaxed(
                    lo1, hi1, d1, lo2, hi2, d2, _adj_min_colsep):
                continue
            inner_bot = max(lo1, lo2)            # bottom of the upper slitlet
            inner_top = min(hi1, hi2)            # top of the lower slitlet
            w_lo = inner_bot - CONFLICT_PURPLE_BUFFER
            w_hi = inner_top + CONFLICT_PURPLE_BUFFER
            rows = set(range(w_lo, w_hi + 1)) if w_hi >= w_lo else set()
            for tg in (tag1, tag2):
                slitlet_window.setdefault(tg, set())
                slitlet_window[tg] |= rows
        # A user PICK that belongs to a slitlet in a CROSS-quadrant conflict is
        # itself a Mask Conflict: its own spectrum is one of the two colliding
        # spectra, even though the partner's light lands on neighbouring
        # shutters rather than this exact one (so it isn't "doubly-hit"). Mark
        # every open shutter of such a slitlet purple so the pick reads as a
        # conflict, not just Masked/clean. (Same-quadrant picks are already
        # coloured via the ±2 window, which intentionally leaves the slitlet's
        # far edges orange — so this is cross-quadrant only.)
        xq_conflict_pick_flats: set[int] = set()
        for _si, _g in enumerate(user_groups_list):
            if ("user", _si) in cross_partners:
                for (_q, _s, _d) in _g:
                    xq_conflict_pick_flats.add(
                        (_q - 1) * 171 * 365 + (_s - 1) * 365 + (_d - 1))
        # Total counts (direct + buffer) drive the silver-edge filter
        # and the stats bar's "Conflict shutters" cell.
        for i in set(user_direct) | set(user_buffer):
            user_counts[i] = user_direct.get(i, 0) + user_buffer.get(i, 0)
        for i in set(stuck_direct) | set(stuck_buffer):
            stuck_counts[i] = stuck_direct.get(i, 0) + stuck_buffer.get(i, 0)

    def _partition_and_project(
        idx_to_alpha: dict[int, float],
    ) -> dict:
        """Project the chosen indices to polygons + emit a fill_alpha
        column parallel to xs / ys / q / s / d."""
        if not idx_to_alpha:
            return dict(xs=[], ys=[], q=[], s=[], d=[], fill_alpha=[])
        idx_arr = np.fromiter(idx_to_alpha.keys(), dtype=np.int64)
        cds = _project_indices_to_cds(idx_arr, pa_v3, fid_pix, jinv)
        cds["fill_alpha"] = [
            float(idx_to_alpha[i]) for i in idx_arr.tolist()
        ]
        return cds

    # Colour rule (alpha always stacks by the number of contaminating
    # sources, so a shutter hit by N sources is N× as opaque):
    #
    # PURPLE (Mask Conflict) — a shutter that is contaminated AND sits in
    # a conflict's bounded ±2-row window (`_purple_here`, below). This
    # covers both the user's own picks and operable shutters caught
    # between two crowding slitlets, but ONLY within that window — beyond
    # it the same band reverts to the warning colours. Purple is therefore
    # always a re-bucketing of an already-contaminated shutter, never an
    # escalation of a clean one.
    #
    # ORANGE (Masked) — contaminated by ≥1 user-open spectrum, outside any
    # conflict window.
    # PINK (Mask Stuck) — contaminated only by stuck-open spectra, outside
    # any conflict window.
    # silver — no contamination (handled by the operable layer, not here).
    user_open_flat = {
        (q - 1) * 171 * 365 + (s - 1) * 365 + (d - 1)
        for (q, s, d) in state["open_shutters"].keys()
    }

    def _purple_here(i: int) -> bool:
        """True iff candidate ``i`` is a genuine Mask Conflict — by either of
        two routes:

        * SAME-quadrant: ``i`` is contaminated by a conflicted slitlet and
          its row falls in that conflict's bounded ±2-row window.
        * CROSS-quadrant: ``i`` is contaminated by BOTH a slitlet and one of
          that slitlet's cross-quadrant conflict partners — two spectra from
          different quadrants genuinely land on it. Bounded to that real
          overlap, never either slitlet's whole band.

        ``hit_sources``, ``slitlet_window``, ``slitlet_quad`` and
        ``cross_partners`` only exist on a cache MISS, which is also the only
        time ``all_idx`` is non-empty — so this is never reached on a hit."""
        if i in xq_conflict_pick_flats:
            return True                           # own pick in a cross-quad conflict
        q_i = i // (171 * 365) + 1
        s_i = (i // 365) % 171 + 1
        srcs = hit_sources.get(i, ())            # sources contaminating i
        for tag in srcs:
            win = slitlet_window.get(tag)
            if win and slitlet_quad.get(tag) == q_i and s_i in win:
                return True                       # same-quad ±2 window
            partners = cross_partners.get(tag)
            if partners and not partners.isdisjoint(srcs):
                return True                       # cross-quad genuine overlap
        return False

    stuck_only_alpha: dict[int, float] = {}
    user_only_alpha: dict[int, float] = {}
    both_alpha: dict[int, float] = {}
    all_idx = (
        set(user_direct) | set(user_buffer)
        | set(stuck_direct) | set(stuck_buffer)
        | xq_conflict_pick_flats
    )
    for i in all_idx:
        n_total_u = user_direct.get(i, 0) + user_buffer.get(i, 0)
        n_total_s = stuck_direct.get(i, 0) + stuck_buffer.get(i, 0)
        n_total = n_total_u + n_total_s
        if _purple_here(i):
            # max(n_total, 1): a cross-quad-conflicted pick may not itself be
            # hit (n_total = 0) but must still render at the base purple alpha.
            both_alpha[i] = min(1.0, base_alpha_both * max(n_total, 1))
        elif i in user_open_flat:
            # A picked shutter that's masked but not in a conflict window
            # reads as orange (Masked), like any user-contaminated shutter.
            if n_total >= 1:
                user_only_alpha[i] = min(1.0, base_alpha_user * n_total)
        elif n_total_u >= 1:
            user_only_alpha[i] = min(1.0, base_alpha_user * n_total)
        elif n_total_s >= 1:
            stuck_only_alpha[i] = min(1.0, base_alpha_stuck * n_total_s)

    # On a cache hit the accumulation + colouring above were skipped (the
    # `if (not _overlap_cached) …` gate and the empty `all_idx` loop), so
    # restore the four view-independent outputs from the cache. On a miss
    # we just computed them — stash them for the next pan/zoom.
    if _overlap_cached:
        (stuck_only_alpha, user_only_alpha,
         both_alpha, all_idx) = _overlap_cache[1]
    else:
        state["_spec_overlap_cache"] = (
            _overlap_sig,
            (stuck_only_alpha, user_only_alpha, both_alpha, all_idx),
        )

    # The colours above are computed globally (view-independent); cull to
    # the visible view ONLY here, at render time, so the polygon count
    # stays bounded while a shutter's colour never depends on the zoom.
    def _iv(alpha):
        return {i: a for i, a in alpha.items() if in_view[i]}
    src_spec_overlap_stuck.data = _partition_and_project(_iv(stuck_only_alpha))
    src_spec_overlap_user.data = _partition_and_project(_iv(user_only_alpha))
    src_spec_overlap_both.data = _partition_and_project(_iv(both_alpha))

    # Aggregate index set for the "operable silver-edge" filter
    # below — every shutter that's affected by ANY overlap (any of
    # the three colours) gets excluded from the silver layer.
    overlap_idx = (
        np.fromiter(all_idx, dtype=np.int64) if all_idx
        else np.empty(0, dtype=np.int64)
    )

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
        # Cull to the visible view (NOT the image cutout) so sources that
        # project beyond the image edge — at x<0 / y<0 or past W/H — still
        # render once that area is on screen; `on_ranges_update` re-runs
        # this on pan/zoom. The margin keeps a marker whose centre sits
        # just outside the view from popping in/out at the edge.
        # Non-finite projections (degenerate WCS far off-field) are dropped.
        _vx0, _vx1, _vy0, _vy1 = view_bbox
        _cm = 50.0
        mask = (
            np.isfinite(x) & np.isfinite(y)
            & (x >= _vx0 - _cm) & (x <= _vx1 + _cm)
            & (y >= _vy0 - _cm) & (y <= _vy1 + _cm)
        )
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
        # shutter's target_id: green ring, thicker line. Unmatched
        # entries take their per-catalog colour (cycling palette so
        # users can tell which catalog a marker belongs to when
        # several are loaded at once).
        matched_target_ids = {
            str(sh.target_id) for sh in state["open_shutters"].values()
            if sh.target_id is not None
        }
        per_source_colors = state.get("catalog_colors")
        if per_source_colors is not None and len(per_source_colors) == len(cat.ra_deg):
            unmatched_colors = np.asarray(per_source_colors)[mask].tolist()
        else:
            unmatched_colors = ["#ffd200"] * len(ids)
        per_source_alphas = state.get("catalog_alphas")
        if per_source_alphas is not None and len(per_source_alphas) == len(cat.ra_deg):
            base_alphas = np.asarray(per_source_alphas)[mask].tolist()
        else:
            base_alphas = [1.0] * len(ids)
        line_colors = [
            "#2e9b3f" if tid in matched_target_ids else fallback
            for tid, fallback in zip(ids, unmatched_colors)
        ]
        line_widths = [
            2.5 if tid in matched_target_ids else 1.5
            for tid in ids
        ]
        # Matched sources stay fully opaque so a "picked" marker is
        # never visually demoted by its catalog's z-depth.
        line_alphas = [
            1.0 if tid in matched_target_ids else a
            for tid, a in zip(ids, base_alphas)
        ]
        # Per-target hover fields beyond ID / RA / Dec / Pr. NaN /
        # missing → empty string in the tooltip rather than "nan".
        def _fmt_num(v):
            try:
                f = float(v)
                if not np.isfinite(f):
                    return ""
                # Trim trailing zeros — keeps `Mag 23` rather than
                # `Mag 23.000000`.
                return f"{f:g}"
            except (TypeError, ValueError):
                return ""
        wt_arr = (np.asarray(cat.weight, dtype=float)
                  if cat.weight is not None
                  and len(cat.weight) == len(cat.ra_deg)
                  else np.full(len(cat.ra_deg), np.nan))
        mag_arr = (np.asarray(cat.mag, dtype=float)
                   if cat.mag is not None else
                   np.full(len(cat.ra_deg), np.nan))
        z_arr = (np.asarray(cat.z, dtype=float)
                 if cat.z is not None
                 else np.full(len(cat.ra_deg), np.nan))
        label_arr = (np.asarray(cat.label, dtype=object)
                     if cat.label is not None
                     else np.array([""] * len(cat.ra_deg), dtype=object))
        # Compact "constraint" indicator: a comma-list of single-
        # letter flags (G=no_gap, B=extend_blue, R=extend_red,
        # P=protect) plus "λ:N" when N required-λ ranges are set.
        # Blank when no constraint applies to this row.
        req_lam = getattr(cat, "required_lam", None)
        no_gap = np.asarray(getattr(cat, "no_gap", []), dtype=bool)
        ext_b = np.asarray(getattr(cat, "extend_blue", []), dtype=bool)
        ext_r = np.asarray(getattr(cat, "extend_red", []), dtype=bool)
        prot = np.asarray(getattr(cat, "protect", []), dtype=bool)
        cent = np.asarray(
            getattr(cat, "centration", []), dtype=object,
        )

        # Single-letter shorthand for the per-target centration label
        # used in the hover-summary glyph. Anything that doesn't match
        # falls back to "?" — defensive; the loader normalises so this
        # branch shouldn't fire in practice.
        _CENTR_SHORT = {
            "UNCONSTRAINED": "U", "ENTIRE_OPEN": "E",
            "MIDPOINT": "M", "CONSTRAINED": "C",
            "TIGHTLY_CONSTRAINED": "T",
        }

        def _constr_glyph(i: int) -> str:
            tags: list[str] = []
            n_req = (len(req_lam[i]) if req_lam is not None
                     and i < len(req_lam)
                     and req_lam[i] is not None
                     else 0)
            if n_req > 0:
                tags.append(f"λ:{n_req}")
            if no_gap.size > i and bool(no_gap[i]):
                tags.append("G")
            if ext_b.size > i and bool(ext_b[i]):
                tags.append("B")
            if ext_r.size > i and bool(ext_r[i]):
                tags.append("R")
            if prot.size > i and bool(prot[i]):
                tags.append("🛡")
            if cent.size > i:
                c = str(cent[i]).strip()
                if c:
                    tags.append("C:" + _CENTR_SHORT.get(c.upper(), "?"))
            return "·".join(tags)

        rows_mask = np.where(mask)[0]
        src_targets.data = dict(
            x=x[mask].tolist(),
            y=y[mask].tolist(),
            id=ids,
            ra=cat.ra_deg[mask].tolist(),
            dec=cat.dec_deg[mask].tolist(),
            pr=cat.priority[mask].tolist(),
            wt=[_fmt_num(wt_arr[i]) for i in rows_mask],
            mag=[_fmt_num(mag_arr[i]) for i in rows_mask],
            z=[_fmt_num(z_arr[i]) for i in rows_mask],
            label=[str(label_arr[i]) for i in rows_mask],
            constr=[_constr_glyph(i) for i in rows_mask],
            line_color=line_colors,
            line_width=line_widths,
            line_alpha=line_alphas,
        )
    else:
        src_targets.data = dict(
            x=[], y=[], id=[], ra=[], dec=[], pr=[],
            wt=[], mag=[], z=[], label=[], constr=[],
            line_color=[], line_width=[], line_alpha=[],
        )

    # Update stats panel and status with the current configuration.
    n_op = len(state["open_shutters"])
    n_tgt_open = len({sh.target_id for sh in state["open_shutters"].values() if sh.target_id})
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

    # Each named cell as an independent HTML span so the user can pick
    # which to show and in what order via Settings → Top stats bar.
    # Keys here must match those in STATS_BAR_CELLS (defined near the
    # Settings UI) so the order/visibility picker stays in sync.
    cells: dict[str, str] = {
        "image": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">Image</span>'
            f'  <span style="{val_style}">{img_label}</span>'
            f'</span>'
        ),
        "radec": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">RA · Dec</span>'
            f'  <span style="{val_style}">{fiducial.ra.deg:.5f} · {fiducial.dec.deg:.5f}</span>'
            f'  <span style="color:#7c8aa0; font-size:11px; margin-left:6px;">'
            f'    ({ra_hms} · {dec_dms})'
            f'  </span>'
            f'</span>'
        ),
        "pa": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">V3 PA</span>'
            f'  <span style="{val_style}">{pa_v3:.2f}°</span>'
            f'  <span style="{label_style}; margin-left:10px;">APA</span>'
            f'  <span style="{val_style}">{apa:.2f}°</span>'
            f'</span>'
        ),
        "disperser": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">Disperser</span>'
            f'  <span style="{val_style}">{state["disperser"]} / {state["filter"]}</span>'
            f'</span>'
        ),
        "open": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">Open</span>'
            f'  <span style="{val_style}">{n_op}</span>'
            f'  <span style="color:#7c8aa0; font-size:11px; margin-left:4px;">'
            f'    across {n_tgt_open} target{"s" if n_tgt_open != 1 else ""}'
            f'  </span>'
            f'</span>'
        ),
        "conflicts": (
            f'<span style="{cell}">'
            f'  <span style="{label_style}">Conflicts</span>'
            f'  <span style="color:{oc}; font-weight:700; font-size:14px; margin-left:4px;">{n_overlap}</span>'
            f'  <span style="color:#7c8aa0; font-size:11px; margin-left:4px;">shutters</span>'
            f'</span>'
        ),
    }
    # Resolved cell order = whatever the Settings → Top stats bar
    # picker says (user-customisable, v1.3.0+). Defaults to the
    # canonical order if state hasn't been touched.
    order = state.get("stats_bar_order") or list(STATS_BAR_DEFAULT_ORDER)
    ordered_html = [cells[k] for k in order if k in cells]
    # The right-most visible cell drops its right border so the bar
    # ends cleanly.
    if ordered_html:
        last = ordered_html[-1]
        ordered_html[-1] = last.replace(
            f'<span style="{cell}">',
            f'<span style="{cell} border-right:none;">',
            1,
        )
    # Multi-config: a bright "CONFIG k/N" chip pinned to the left of the
    # always-visible top bar so the user never loses track of which config
    # a manual open lands in (accent matches the Pointing-tab banner).
    # The chip is also a CONTROL: clicking it cycles the active config
    # (1→2→…→1) so you can switch without leaving the canvas. The onclick
    # calls the global installed by `_config_chip_install_js`.
    n_cfg = int(state.get("n_configs", 1))
    config_chip = ""
    if n_cfg > 1:
        act = int(state.get("active_config", 0)) + 1
        accent = _config_color(act - 1)
        config_chip = (
            f'<span onclick="window.__vmpt_cycle_config && '
            f'window.__vmpt_cycle_config()" '
            f'title="Click to switch the active config '
            f'(now Config {act} of {n_cfg})" '
            f'style="display:inline-flex; align-items:center; cursor:pointer; '
            f'user-select:none; background:{accent}; color:white; '
            f'font-weight:700; padding:3px 10px; border-radius:4px; '
            f'margin:0 8px 0 2px; font-size:12px; white-space:nowrap;">'
            f'⇄ CONFIG {act} / {n_cfg}</span>'
        )
    stats_div.text = (
        '<div style="display:flex; flex-wrap:wrap; align-items:center; gap:0;">'
        + config_chip
        + "".join(ordered_html)
        + '</div>'
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
    # Lock the canvas frame pixel dimensions to the user-set values
    # so 1 data unit renders the same in x and y. Default 800x800;
    # the legacy "frame_max" key (pre-split) still drives both axes
    # if present, for session-reload backwards compatibility.
    if W > 0 and H > 0:
        # User-adjustable canvas dimensions (state-backed so values
        # survive image reloads). Default 800x600; legacy
        # "frame_max" key (pre-X/Y-split) still drives both axes if
        # present, for session-reload backwards compatibility.
        legacy = state.get("frame_max")
        frame_x = int(state.get("frame_x", legacy or 800))
        frame_y = int(state.get("frame_y", legacy or 600))
        fig.frame_width = max(100, frame_x)
        fig.frame_height = max(100, frame_y)
        # **Per-pixel square** invariant: every data-unit-in-x must
        # render at the same screen-pixel size as every
        # data-unit-in-y, regardless of the (frame_x, frame_y) the
        # user chose. `match_aspect=True` on the figure can only
        # enforce this when the data range aspect equals the frame
        # aspect — so we PRE-COMPUTE the data ranges to match the
        # frame aspect by inflating the short axis symmetrically
        # around the image centre.
        #
        # The image stays anchored at (x=0, y=0) and renders at its
        # native W:H aspect inside the inflated data area. The
        # inflated "empty" padding shows on whichever axis was
        # short for the canvas — user can pan into it. The MSA
        # outline + catalog overlay sit at their WCS-derived
        # coords (in the same pixel-coord frame as the image), so
        # they keep their correct geometric ratios at any (X, Y).
        #
        # NB: We need explicit Float values here — Bokeh 3.7's
        # DataRange1d.start/end rejects None even though they're
        # "Nullable". The previous fix using start=None blew up
        # at runtime; this version pins explicit floats that are
        # geometrically consistent with the chosen frame so
        # match_aspect's check passes trivially.
        canvas_aspect = float(frame_x) / float(frame_y)
        image_aspect = float(W) / float(H)
        if image_aspect < canvas_aspect:
            # Image narrower than canvas — pad x_range on both sides.
            target_x_span = canvas_aspect * H
            pad_x = (target_x_span - W) / 2.0
            x_start, x_end = -pad_x, W + pad_x
            y_start, y_end = 0.0, float(H)
        else:
            # Image wider than canvas — pad y_range on both sides.
            target_y_span = W / canvas_aspect
            pad_y = (target_y_span - H) / 2.0
            x_start, x_end = 0.0, float(W)
            y_start, y_end = -pad_y, H + pad_y
        fig.x_range.update(start=x_start, end=x_end)
        fig.y_range.update(start=y_start, end=y_end)
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


def _load_fits_from_path(path: str, force_recenter: bool = False,
                         on_complete=None) -> None:
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
        if on_complete is not None:
            curdoc().add_next_tick_callback(on_complete)


def _load_jpg_pair_from_paths(
    jpg_path: str, sidecar_path: str, force_recenter: bool = False,
    on_complete=None,
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
        if on_complete is not None:
            curdoc().add_next_tick_callback(on_complete)


def _catalog_alpha_for_depth(depth: int) -> float:
    """Per-catalog marker alpha based on z-order depth (0 = top).

    Earlier-loaded catalogs sit at the top with full opacity; each
    additional layer below fades slightly so a dense reference
    catalog underneath a curated target list doesn't drown it out.
    The floor of 0.35 keeps even deeply-buried catalogs visible
    enough to be useful."""
    return max(0.35, 1.0 - 0.20 * depth)


def _rebuild_merged_catalog() -> None:
    """Combine all ENABLED entries in state['catalogs'] into the single
    `state['catalog']` cache that the rest of the app (overlay, source
    matching, export) reads. Disabled / removed catalogs drop out
    automatically.

    Z-order: Bokeh's scatter renders later indices ON TOP, so we put
    the earliest-loaded catalog LAST in the merged arrays. The order
    of `state['catalogs']` is the visual stack (index 0 = top, index
    n−1 = bottom). Reordering the list via the ▲/▼ buttons in
    `_render_catalog_list` changes the on-screen z-order.

    Also builds `state['catalog_colors']` (per-source marker colour)
    and `state['catalog_alphas']` (per-source line_alpha, decayed by
    z-depth). Both are parallel arrays the overlay uses without
    needing to walk `state['catalogs']` for every row."""
    enabled = [e for e in state["catalogs"] if e["enabled"]]
    if not enabled:
        state["catalog"] = None
        state["catalog_colors"] = None
        state["catalog_alphas"] = None
        try:
            _refresh_opt_status_div()
        except NameError:
            pass
        return

    # Depth is the catalog's index in state['catalogs'] (NOT the
    # enabled-only list) — so the alpha stays stable when the user
    # toggles a layer off and on without reordering.
    depth_of_entry = {
        id(entry): idx for idx, entry in enumerate(state["catalogs"])
    }

    # Single-catalog fast path: don't re-allocate the merged arrays,
    # just hand the original Catalog through (preserves the int dtype
    # of `ids`, which downstream callers depend on for lookups).
    if len(enabled) == 1:
        only = enabled[0]
        cat = only["catalog"]
        a = _catalog_alpha_for_depth(depth_of_entry[id(only)])
        state["catalog"] = cat
        state["catalog_colors"] = np.full(len(cat.ra_deg), only["color"], dtype=object)
        state["catalog_alphas"] = np.full(len(cat.ra_deg), a, dtype=float)
        return

    # Walk catalogs in REVERSE state-list order so the merged data
    # arrays have last-loaded sources first (drawn at the back) and
    # first-loaded sources last (drawn on top).
    ordered_for_draw = list(reversed(enabled))

    cats = [e["catalog"] for e in ordered_for_draw]
    ids = np.concatenate([np.asarray(c.ids, dtype=object) for c in cats])
    ra = np.concatenate([c.ra_deg for c in cats])
    dec = np.concatenate([c.dec_deg for c in cats])
    pri = np.concatenate([c.priority for c in cats])
    # Weight: backfill NaN for catalogs without a weight column so the
    # merged array is always length-aligned to ra. Needed by Meritocracy
    # mode and by the collision-protection "By weight ≥" rule.
    weight = np.concatenate([
        (np.asarray(c.weight, dtype=float) if getattr(c, "weight", None)
         is not None and len(c.weight) == len(c.ra_deg)
         else np.full(len(c.ra_deg), np.nan, dtype=float))
        for c in cats
    ])
    mag = np.concatenate([c.mag for c in cats])
    z = np.concatenate([c.z for c in cats])
    label = np.concatenate([
        np.asarray(c.label, dtype=object) if c.label is not None
        else np.array([""] * len(c.ra_deg), dtype=object)
        for c in cats
    ])
    # Per-target spectral constraints (v1.3.0+). Backfill defaults
    # when a catalog was loaded before these fields existed (e.g.
    # older session round-trips) so the merged arrays always match
    # the row count.
    def _pad_bool_per_cat(name: str) -> np.ndarray:
        chunks = []
        for c in cats:
            arr = getattr(c, name, None)
            arr = np.asarray(arr if arr is not None else [], dtype=bool)
            if arr.size != len(c.ra_deg):
                arr = np.zeros(len(c.ra_deg), dtype=bool)
            chunks.append(arr)
        return np.concatenate(chunks) if chunks else np.array([], dtype=bool)
    def _pad_lam_req() -> np.ndarray:
        chunks = []
        for c in cats:
            arr = getattr(c, "required_lam", None)
            if arr is None or len(arr) != len(c.ra_deg):
                # 1D object array of [] — `np.array([[]…], dtype=O)`
                # would build a 2D shape=(n, 0) array instead, which
                # then trips bool() ambiguity downstream.
                tmp = np.empty(len(c.ra_deg), dtype=object)
                for i in range(len(c.ra_deg)):
                    tmp[i] = []
                arr = tmp
            else:
                arr = np.asarray(arr, dtype=object)
            chunks.append(arr)
        return np.concatenate(chunks) if chunks else np.array([], dtype=object)
    required_lam_merged = _pad_lam_req()
    no_gap_merged = _pad_bool_per_cat("no_gap")
    extend_blue_merged = _pad_bool_per_cat("extend_blue")
    extend_red_merged = _pad_bool_per_cat("extend_red")
    protect_merged = _pad_bool_per_cat("protect")

    # Per-target centration override (v1.3.1+). Strings (not bools),
    # so it gets its own merge helper — pad missing/short entries
    # with "" (= use the global optimizer setting).
    def _pad_centration_per_cat() -> np.ndarray:
        chunks = []
        for c in cats:
            arr = getattr(c, "centration", None)
            arr = (np.asarray(arr, dtype=object) if arr is not None
                   else np.array([], dtype=object))
            if arr.size != len(c.ra_deg):
                arr = np.array([""] * len(c.ra_deg), dtype=object)
            chunks.append(arr)
        return (np.concatenate(chunks) if chunks
                else np.array([], dtype=object))
    centration_merged = _pad_centration_per_cat()

    # Per-target multi-config cap (v1.4.0+). Float with NaN = unset; pad
    # missing/short catalogs with NaN so they inherit the global default.
    def _pad_max_configs_per_cat() -> np.ndarray:
        chunks = []
        for c in cats:
            arr = getattr(c, "max_configs", None)
            arr = (np.asarray(arr, dtype=float) if arr is not None
                   else np.array([], dtype=float))
            if arr.size != len(c.ra_deg):
                arr = np.full(len(c.ra_deg), np.nan, dtype=float)
            chunks.append(arr)
        return (np.concatenate(chunks) if chunks
                else np.array([], dtype=float))
    max_configs_merged = _pad_max_configs_per_cat()
    colors = np.concatenate([
        np.full(len(e["catalog"].ra_deg), e["color"], dtype=object)
        for e in ordered_for_draw
    ])
    alphas = np.concatenate([
        np.full(
            len(e["catalog"].ra_deg),
            _catalog_alpha_for_depth(depth_of_entry[id(e)]),
            dtype=float,
        )
        for e in ordered_for_draw
    ])
    state["catalog"] = Catalog(
        ids=ids, ra_deg=ra, dec_deg=dec, priority=pri, weight=weight,
        mag=mag, z=z, label=label,
        required_lam=required_lam_merged,
        no_gap=no_gap_merged,
        extend_blue=extend_blue_merged,
        extend_red=extend_red_merged,
        protect=protect_merged,
        centration=centration_merged,
        max_configs=max_configs_merged,
        source_path=" + ".join(c.source_path for c in cats),
    )
    state["catalog_colors"] = colors
    state["catalog_alphas"] = alphas
    # Keep the optimizer's status line in sync with what the catalog
    # actually has — was previously a static "Load a catalog with
    # priorities…" message that confused users with catalogs already
    # loaded. Defined later in the file; gated on existence so the
    # initial autoload doesn't hit a forward-reference error.
    try:
        _refresh_opt_status_div()
    except NameError:
        pass
    # Refresh the per-target centration-override hint in the optimizer
    # modal. Gated for the same forward-reference reason — defined
    # later in the file.
    try:
        _refresh_centration_override_hint()
    except NameError:
        pass


def _assign_catalog_color(index: int) -> str:
    """Pick a palette entry for the catalog at the given list index.
    Cycles when more catalogs are loaded than palette entries."""
    return CATALOG_COLOR_PALETTE[index % len(CATALOG_COLOR_PALETTE)]


def _render_catalog_list() -> None:
    """Rebuild the catalog_list_column from state['catalogs'].

    Each row carries: a colour chip (marker colour key), the
    catalog's name + source-count, an enable/disable checkbox, ▲/▼
    buttons for reorder, and × to delete.

    The list order IS the visual z-order — index 0 sits on top of the
    canvas, the last index is on the bottom (and slightly fainter).
    Click ▲ to move a catalog up the stack (closer to the top) or ▼
    to push it down. The first row's ▲ and the last row's ▼ are
    disabled."""
    rows = []
    n_entries = len(state["catalogs"])
    for idx, entry in enumerate(state["catalogs"]):
        name = entry["name"]
        n = len(entry["catalog"].ra_deg)
        color = entry.get("color", "#ffd200")
        # Coloured chip — small Div with a solid background that matches
        # the marker colour drawn on the canvas. Gives the user a
        # visual key linking each list row to its on-image markers.
        chip = Div(
            text=(
                f'<div style="width:14px; height:14px; border-radius:50%; '
                f'background:{color}; border:1px solid #333; '
                f'display:inline-block; vertical-align:middle;"></div>'
            ),
            width=22, height=24,
        )
        cb = CheckboxGroup(
            labels=[f"{name} ({n})"],
            active=[0] if entry["enabled"] else [],
            width=SIDEBAR_W - 170,
        )
        up_btn = Button(
            label="▲", button_type="default",
            width=28, height=28, disabled=(idx == 0),
        )
        down_btn = Button(
            label="▼", button_type="default",
            width=28, height=28, disabled=(idx == n_entries - 1),
        )
        del_btn = Button(
            label="×", button_type="warning",
            width=30, height=28,
        )

        def _make_toggle(i):
            def _cb(attr, old, new):
                if not (0 <= i < len(state["catalogs"])):
                    return
                state["catalogs"][i]["enabled"] = (0 in new)

                def _apply():
                    try:
                        _rebuild_merged_catalog()
                        _rebuild_shutter_catalog_index()
                        refresh_overlays()
                    finally:
                        _hide_loading()
                # Re-rendering thousands of catalog circles takes a beat;
                # show the spinner first (deferred build) for big catalogs.
                n = len(state["catalogs"][i]["catalog"].ra_deg)
                if n > _LARGE_CATALOG_N:
                    _show_loading(f"Rendering {n} sources…")
                    _deferred(_apply)
                else:
                    _apply()
            return _cb

        def _make_delete(i, label_name):
            def _cb():
                if 0 <= i < len(state["catalogs"]):
                    del state["catalogs"][i]
                    _rebuild_merged_catalog()
                    _rebuild_shutter_catalog_index()
                    _render_catalog_list()
                    refresh_overlays()
                    _set_status(f"Removed catalog: {label_name}", "ok")
            return _cb

        def _make_move(i, direction):
            """direction: -1 = up (toward index 0, top of stack);
            +1 = down (toward end, bottom of stack)."""
            def _cb():
                j = i + direction
                cats = state["catalogs"]
                if 0 <= i < len(cats) and 0 <= j < len(cats):
                    cats[i], cats[j] = cats[j], cats[i]
                    _rebuild_merged_catalog()
                    _rebuild_shutter_catalog_index()
                    _render_catalog_list()
                    refresh_overlays()
            return _cb

        cb.on_change("active", _make_toggle(idx))
        up_btn.on_click(_make_move(idx, -1))
        down_btn.on_click(_make_move(idx, +1))
        del_btn.on_click(_make_delete(idx, name))
        rows.append(row(chip, cb, up_btn, down_btn, del_btn))
    catalog_list_column.children = rows


def _load_catalog_from_path(path: str, on_complete=None) -> None:
    # Guard against a non-file path (empty, or a directory like '.') before
    # handing it to the parser — otherwise astropy raises IsADirectoryError
    # and dumps a scary traceback for what is really just "nothing to load".
    # Belt-and-braces: the UI handlers also pre-check, but a stray value
    # (e.g. a path input that ended up '.') could still reach here.
    if not path or not Path(path).is_file():
        if path:
            _set_status(f"Not a catalog file — skipped: {path}", "warn")
        _hide_loading()
        if on_complete is not None:
            curdoc().add_next_tick_callback(on_complete)
        return
    try:
        cat = load_catalog(path)
        name = Path(path).name
        # Refuse duplicates: if the same path is already loaded, just
        # re-enable it rather than appending a redundant copy.
        for entry in state["catalogs"]:
            if entry["catalog"].source_path == cat.source_path:
                entry["enabled"] = True
                _rebuild_merged_catalog()
                _rebuild_shutter_catalog_index()
                _render_catalog_list()
                refresh_overlays()
                _set_status(
                    f"Catalog already loaded — re-enabled {name} "
                    f"({len(cat.ra_deg)} targets).", "ok",
                )
                return
        state["catalogs"].append({
            "name": name,
            "catalog": cat,
            "enabled": True,
            "color": _assign_catalog_color(len(state["catalogs"])),
        })
        _rebuild_merged_catalog()
        _rebuild_shutter_catalog_index()
        _render_catalog_list()
        refresh_overlays()
        n_total = sum(len(e["catalog"].ra_deg) for e in state["catalogs"]
                      if e["enabled"])
        _set_status(
            f"Catalog added: {name} ({len(cat.ra_deg)} targets). "
            f"{len(state['catalogs'])} catalog(s) loaded, {n_total} "
            f"active targets total.", "ok",
        )
    except Exception as e:  # noqa: BLE001
        _set_status(f"Catalog load failed: {e}", "err")
        traceback.print_exc()
    finally:
        _hide_loading()
        if on_complete is not None:
            curdoc().add_next_tick_callback(on_complete)


# Path-based callbacks (primary input for a local tool).
# Slow loads are deferred to the next tick so the loading banner renders first.
def on_fits_path(attr, old, new):
    if state.get("_autoload_active"):
        return
    p = fits_path_input.value.strip()
    if not p:
        return
    if not Path(p).exists():
        _set_status(f"FITS path not found: {p}", "err")
        return
    _show_loading(f"Loading FITS: {Path(p).name}…")
    _close_load_image_modal()
    _deferred(_load_fits_from_path, p)


def on_sidecar_path(attr, old, new):
    if state.get("_autoload_active"):
        return
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
        _close_load_image_modal()
        _deferred(_load_jpg_pair_from_paths, jpg_p, p)
    else:
        _set_status("Sidecar set. Now enter a JPG path.", "info")


def on_jpg_path(attr, old, new):
    if state.get("_autoload_active"):
        return
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
    _close_load_image_modal()
    _deferred(_load_jpg_pair_from_paths, jpg_p, side_p)


def on_catalog_path(attr, old, new):
    """Path change → load + append the catalog. Convenient for the
    pastes-a-path flow. The explicit "Add" button is a parallel
    trigger so the user can re-load a catalog that was previously
    removed without first clearing the path."""
    if state.get("_autoload_active"):
        return
    p = catalog_path_input.value.strip()
    if not p:
        return
    if not Path(p).is_file():
        _set_status(f"Catalog path is not a file: {p}", "err")
        return
    _show_loading(f"Loading catalog: {Path(p).name}…")
    _close_load_catalog_modal()
    _deferred(_load_catalog_from_path, p)


def on_catalog_add():
    """Explicit Add button — re-uses whatever's in catalog_path_input."""
    p = catalog_path_input.value.strip()
    if not p:
        _set_status("Set a catalog path first.", "warn")
        return
    if not Path(p).is_file():
        _set_status(f"Catalog path is not a file: {p}", "err")
        return
    _show_loading(f"Loading catalog: {Path(p).name}…")
    _close_load_catalog_modal()
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


def _rebuild_shutter_catalog_index_core() -> None:
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


def _rebuild_shutter_catalog_index() -> None:
    """Rebuild the source↔shutter footprint index, then re-derive the
    catalog tags of any *raw-mask* open shutters (loaded from a shutter
    CSV, role 'manual') so the MPT catalog reflects what those shutters
    observe at the current pointing — no optimizer run required, and kept
    in sync as the pointing / PA / catalog changes. Deliberate picks
    (manual clicks, optimizer results) are left untouched."""
    _rebuild_shutter_catalog_index_core()
    _retag_manual_opens_from_catalog()


def _shutter_source_id(q: int, s: int, d: int) -> str | None:
    """Return the catalog source id (as a string) that falls inside this
    shutter, or None if none does. If multiple sources land in the same
    shutter we return the first one (catalog order)."""
    bucket = state.get("shutter_to_catids") or {}
    ids = bucket.get((int(q), int(s), int(d)))
    if not ids:
        return None
    return str(ids[0])


def _shutter_source_ids(q: int, s: int, d: int) -> list[str]:
    """Every catalog source id (as strings, catalog order) whose footprint
    falls inside this shutter at the current pointing — usually 0 or 1, but
    two very close sources can share one shutter (v1.4.0 records them all)."""
    bucket = state.get("shutter_to_catids") or {}
    return [str(t) for t in (bucket.get((int(q), int(s), int(d))) or [])]


def _open_shutter_ids(sh) -> list[str]:
    """All source ids attributed to an OpenShutter. Prefers the snapshotted
    ``target_ids`` list; falls back to the scalar ``target_id`` so shutters
    loaded from pre-1.4.0 sessions still report their source."""
    ids = [str(t) for t in (getattr(sh, "target_ids", None) or []) if str(t) != ""]
    if not ids:
        tid = getattr(sh, "target_id", None)
        if tid is not None and str(tid) != "":
            ids = [str(tid)]
    return ids


def _retag_manual_opens_from_catalog() -> int:
    """Re-derive the catalog source ids of the active config's *raw-mask*
    open shutters from the current footprint index.

    A shutter mask loaded from CSV carries no source ids (every shutter is
    anonymous, ``role='manual'``), so the MPT catalog would be empty until
    the optimizer ran. This attributes each such shutter the catalog
    source(s) that fall inside it at the current pointing — exactly the
    match a hand-pick or the optimizer would make — so the MPT catalog
    populates straight from a loaded mask and stays in sync when the
    pointing / PA / catalog changes.

    Only ``role='manual'`` shutters are touched; deliberate picks (manual
    clicks and optimizer results, role 'target'/'sky', which carry their
    own ids) are never clobbered. Returns the number of raw-mask shutters
    that now carry at least one catalog id.
    """
    opens = state.get("open_shutters") or {}
    if not opens:
        return 0
    matched = 0
    for key, sh in list(opens.items()):
        if getattr(sh, "role", None) != "manual":
            continue
        q, s, d = key
        ids = _shutter_source_ids(q, s, d)
        primary = ids[0] if ids else None
        # Replace (don't mutate in place) so undo snapshots holding the old
        # object are unaffected; skip if nothing actually changed.
        if _open_shutter_ids(sh) != ids or getattr(sh, "target_id", None) != primary:
            opens[key] = OpenShutter(q=q, s=s, d=d, target_id=primary,
                                     role="manual", target_ids=list(ids))
        if ids:
            matched += 1
    # Keep an open MPT-catalog viewer live as the pointing/catalog changes.
    try:
        if mpt_view_modal_card.visible:
            _mpt_view_refresh()
    except NameError:
        pass
    return matched


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
        key = (q, s, d)
        # Every catalog source physically inside THIS shutter (snapshot).
        shutter_ids = _shutter_source_ids(q, s, d)
        # Make sure the anchor's chosen primary is recorded on the shutter it
        # sits in, even if a pointing-rounding mismatch kept it out of the
        # footprint bucket.
        if role == "target" and target_id is not None and str(target_id) not in shutter_ids:
            shutter_ids = [str(target_id)] + shutter_ids
        existing = state["open_shutters"].get(key)
        if existing is not None:
            # Shutter already open (e.g. an adjacent slitlet, or two optimizer
            # targets that mapped to the same shutter): MERGE sources rather
            # than clobber. Replace (don't mutate) so undo snapshots that hold
            # the old object are unaffected.
            merged = list(dict.fromkeys(_open_shutter_ids(existing) + shutter_ids))
            state["open_shutters"][key] = OpenShutter(
                q=q, s=s, d=d,
                target_id=(existing.target_id
                           if existing.target_id is not None else target_id),
                role=existing.role,
                target_ids=merged,
            )
        else:
            state["open_shutters"][key] = OpenShutter(
                q=q, s=s, d=d, target_id=target_id, role=role,
                target_ids=shutter_ids,
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


def _toggle_single_shutter_at(q: int, s: int, d: int) -> None:
    """Toggle exactly ONE shutter open↔closed, independent of the
    N-shutter slitlet setting. Used by the hover + spacebar shortcut.

    - If the shutter is open, close just that shutter (its slitlet
      siblings, if any, stay open — this targets the single hovered cell).
    - If it's closed and operable, open just that one shutter (N=1), with
      the usual auto-match to a catalog source inside it.
    - A closed, non-operable shutter (stuck-open / failed-closed) is left
      alone with a warning.
    """
    key = (q, s, d)
    is_open = key in state["open_shutters"]
    if not is_open and not OPERABLE[q - 1, s - 1, d - 1]:
        _set_status(f"Shutter ({q},{s},{d}) isn't operable — can't open it.",
                    "warn")
        return
    _push_history()
    if is_open:
        state["open_shutters"].pop(key)
        _set_status(f"Closed shutter ({q},{s},{d}).", "ok")
    else:
        # Force a single-shutter open regardless of the slitlet-height
        # setting by borrowing _add_slitlet with N temporarily pinned to 1
        # (keeps its catalog auto-match + merge behaviour).
        saved_h = state["slitlet_height"]
        state["slitlet_height"] = 1
        try:
            _add_slitlet(q, s, d, target_id=None)
        finally:
            state["slitlet_height"] = saved_h
        opened = state["open_shutters"].get(key)
        auto_id = opened.target_id if opened else None
        if auto_id:
            _set_status(f"Opened single shutter ({q},{s},{d}); "
                        f"matched catalog source {auto_id}.", "ok")
        else:
            _set_status(f"Opened single shutter ({q},{s},{d}).", "ok")
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
    _set_open_shutters(state["history"].pop())
    refresh_overlays()
    _set_status("Undid last action.", "ok")


def on_clear():
    if not state["open_shutters"]:
        return
    _push_history()
    _set_open_shutters({})
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
    used_int_ids: set[int] = set()
    # Preserve the original (possibly string) catalog token in the Label
    # column. Output No_cat is always an int derived via _to_int_id.
    # Match unsigned digit runs only — hyphens in catalog IDs like
    # "RJ0600-10274-P0" are separators, NOT minus signs (taking them
    # as signs would turn the source number negative).
    _digit_run_re = re.compile(r"\d+")

    def _to_int_id(raw, taken: set[int]) -> int:
        """Coerce an arbitrary source ID into a unique positive integer.

        Strategy:
          1. If the raw value parses directly as an int, use it.
          2. Otherwise pick the **largest** digit run in the string —
             in JWST catalog IDs like "RJ0600-10274-P0" the unique
             source number (10274) is always the biggest integer in
             the token, dwarfing the field prefix (0600 = 600) and
             priority-class suffix (P0 = 0). "Largest" beats "longest"
             when two runs share a length (0600 vs 8846 → 8846 wins).
          3. On collision with a previously-used ID, walk forward to
             the next free integer.
        """
        s = str(raw).strip()
        candidates: list[int] = []
        try:
            candidates.append(int(s))
        except (ValueError, TypeError):
            pass
        for r in _digit_run_re.findall(s):
            try:
                n = int(r)
                if n > 0:
                    candidates.append(n)
            except ValueError:
                pass
        # Largest first — that's the source ID in conventional catalogs.
        candidates.sort(reverse=True)
        for c in candidates:
            if c not in taken:
                return c
        # Last resort: next free positive integer.
        n = max((c for c in taken if c > 0), default=0) + 1
        while n in taken:
            n += 1
        return n

    def _cat_num(arr, k) -> float | None:
        """A finite float from catalog column `arr` at row `k`, else None."""
        try:
            v = float(arr[k])
        except (TypeError, ValueError, IndexError):
            return None
        return v if np.isfinite(v) else None

    def _cat_lookup(tid: str) -> dict | None:
        """Look up a catalog row by id. Returns a dict with ra/dec, priority
        (`pr`), `weight`, `mag`, `z`, and `label`. `label` is the catalog's
        `label`/`name` value when available, else the raw ID (preserving the
        trace from the MPT integer ID back to the user's source name)."""
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
        label_val = label_val.strip() or tid
        return {
            "ra": float(cat.ra_deg[k]), "dec": float(cat.dec_deg[k]),
            "pr": pr,
            "weight": _cat_num(getattr(cat, "weight", []), k),
            "mag": _cat_num(getattr(cat, "mag", []), k),
            "z": _cat_num(getattr(cat, "z", []), k),
            "label": label_val,
        }

    def _push_target_row(
        raw_tid: str | None,
        ra_d: float,
        dec_d: float,
        *,
        pr: int,
        primary: int,
        weight: float | None = None,
        mag: float | None = None,
        z: float | None = None,
        label: str,
    ) -> int:
        """Append a target row. Always returns the assigned integer ID.

        If `raw_tid` is None, a fresh sequential integer is generated.
        Otherwise the integer is derived from `raw_tid` via _to_int_id
        (e.g. "RJ0600-10274-P0" → 10274).

        `Pr` (priority) feeds the eMPT observed_targets.cat; the APT `.cat`
        uses `Weight` (the source's weight, falling back to priority) plus
        the `Primary` flag (1 = catalog source, 0 = vMPT-synthesised) and
        optional `Magnitude`/`Redshift` carried from the input catalog.
        """
        if raw_tid is None:
            target_no = max(used_int_ids, default=0) + 1
            while target_no in used_int_ids:
                target_no += 1
        else:
            target_no = _to_int_id(raw_tid, used_int_ids)
        used_int_ids.add(target_no)
        targets_rows.append({
            "No_cat": target_no,
            "Pr": int(pr),
            "Weight": (weight if weight is not None else pr),
            "Primary": int(primary),
            "Magnitude": mag,
            "Redshift": z,
            "ra_deg": ra_d,
            "dec_deg": dec_d,
            "label": label,
        })
        return target_no

    # Step 1: open shutters with a known catalog source. One output row
    # per source — when two very close sources share one shutter (v1.4.0)
    # BOTH are recorded, not just the slitlet's primary. The output `.cat`
    # row's integer ID is derived from the original catalog token (digit-
    # run extraction); the original token survives in the Label column.
    for (q, s, d), sh in state["open_shutters"].items():
        tids = _open_shutter_ids(sh)
        if not tids:
            fallback = _shutter_source_id(q, s, d)
            tids = [fallback] if fallback else []
        for tid in tids:
            if not tid or str(tid) in real_ids_seen:
                continue
            info = _cat_lookup(str(tid))
            if info is None:
                continue  # synthesise later from geometry
            real_ids_seen.add(str(tid))
            _push_target_row(
                str(tid), info["ra"], info["dec"],
                pr=info["pr"], primary=1, weight=info["weight"],
                mag=info["mag"], z=info["z"], label=info["label"],
            )

    # Step 1b: append ALL OTHER sources from the loaded input catalog —
    # whether or not they're inside any open shutter. Output is a
    # strict superset of the input list so collaborators see the full
    # context (rejected / unobserved targets included).
    if cat is not None:
        for i in range(len(cat.ra_deg)):
            tid = str(cat.ids[i])
            if tid in real_ids_seen:
                continue
            real_ids_seen.add(tid)
            pr = int(cat.priority[i]) if np.isfinite(cat.priority[i]) else 1
            label_val = ""
            try:
                if cat.label is not None and i < len(cat.label):
                    label_val = str(cat.label[i]).strip()
            except (AttributeError, IndexError, TypeError):
                label_val = ""
            if not label_val:
                label_val = tid  # preserve the original ID string as Label
            _push_target_row(
                tid,
                float(cat.ra_deg[i]),
                float(cat.dec_deg[i]),
                pr=pr, primary=1,
                weight=_cat_num(getattr(cat, "weight", []), i),
                mag=_cat_num(getattr(cat, "mag", []), i),
                z=_cat_num(getattr(cat, "z", []), i),
                label=label_val,
            )

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
            fake_id = str(_push_target_row(
                None, ra_d, dec_d, pr=5, primary=0, label="vMPT_synth",
            ))
            # Tag every shutter in the run with this fake id so later
            # exporters (MPT plan primaryIds) see consistent target IDs.
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
    # APT-facing filenames: target-prefixed + role-suffixed so it's obvious
    # which file APT/MPT loads. The .cat basename equals the plan JSON's
    # catalog.name (apt_catalog_basename), so APT's filename-derived default
    # lines up with the plan automatically.
    _cat_src = getattr(cat, "source_path", None) if cat is not None else None
    mpt_catalog_name = apt_catalog_basename(_cat_src) + ".cat"
    mpt_plan_name = apt_plan_basename(_cat_src) + ".json"

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
        # 2b) Per-config eMPT artifacts (v1.4.0). The top-level files
        # above describe the ACTIVE config; when a plan has >1 config,
        # each also gets a config_N/ subdir with its own one-row-per-
        # source observed list + shutter mask (eMPT outputs are
        # per-pointing). Single-config plans skip this entirely so their
        # bundle is unchanged from ≤1.3.x.
        n_cfg_out = int(state.get("n_configs", 1))
        n_sub = 0
        if n_cfg_out > 1:
            _save_active_config_pointing()
            for ci in range(min(n_cfg_out, len(state["configs"]))):
                cobj = state["configs"][ci]
                copen = list(cobj["open_shutters"].values())
                sub = out_dir / f"config_{ci + 1}"
                sub.mkdir(parents=True, exist_ok=True)
                crows: list = []
                cused: set = set()
                cseen: set = set()
                for sh in copen:
                    for tid in _open_shutter_ids(sh):
                        if tid in cseen:
                            continue
                        cseen.add(tid)
                        info = _cat_lookup(str(tid))
                        if info is None:
                            continue
                        no = _to_int_id(str(tid), cused)
                        cused.add(no)
                        crows.append({"No_cat": no, "Pr": info["pr"],
                                      "ra_deg": info["ra"],
                                      "dec_deg": info["dec"],
                                      "label": info["label"]})
                write_observed_targets_cat(
                    str(sub / EMPT_OBSERVED_FILENAME), crows)
                write_shutter_mask_csv(
                    str(sub / EMPT_SHUTTER_MASK_FILENAME),
                    copen, OPERABLE, REASON)
                cra, cdec, cpa = (cobj.get("ra_deg"), cobj.get("dec_deg"),
                                  cobj.get("pa_v3"))
                if cra is not None and cdec is not None and cpa is not None:
                    cpa_ap = (float(cpa) + V3_IDL_Y_ANGLE) % 360.0
                    write_pointing_summary_txt(
                        str(sub / EMPT_POINTING_FILENAME),
                        Pointing(ra_deg=float(cra), dec_deg=float(cdec),
                                 apa_v3_deg=float(cpa), pa_ap_deg=cpa_ap),
                        state["disperser"], state["filter"],
                        n_targets_total=(len(cat.ra_deg)
                                         if cat is not None else 0),
                        n_targets_accepted=len(crows))
                n_sub += 1
        # 3) MPT plan + vMPT workspace sidecar so the bundle round-trips
        # cleanly: Session → Load on either file restores everything.
        export_session_json(_build_current_session(), str(out_dir / mpt_plan_name))
        # 4) README telling the user exactly how to load this into APT/MPT.
        try:
            (out_dir / "README.md").write_text(bundle_readme_text(
                catalog_filename=mpt_catalog_name,
                catalog_name=apt_catalog_basename(_cat_src),
                plan_filename=mpt_plan_name,
                n_configs=n_cfg_out,
            ))
        except OSError:
            pass  # README is a convenience; never fail the export over it
        sub_note = (f" + {n_sub} config_N/ subdir(s)" if n_sub else "")
        _set_status(
            f"Wrote bundle to {out_dir} — {mpt_plan_name} + "
            f"{mpt_catalog_name} (APT import) + README.md + "
            f"{WORKSPACE_FILENAME} (vMPT state) + 3 eMPT_* files{sub_note}. "
            f"{len(targets_rows)} targets, {len(open_list)} open shutters.",
            "ok", clear_after=18,
        )
        # Pre-fill the Session-load input so the user can re-load with one click.
        session_load_path_input.value = str(out_dir / mpt_plan_name)
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
    # `catalog_paths` is the multi-catalog source of truth; the
    # legacy single-`catalog_path` field is kept in lockstep (set to
    # the first enabled entry) so vMPT 1.0 bundles still round-trip.
    catalog_entries = [
        {
            "path": e["catalog"].source_path,
            "enabled": bool(e.get("enabled", True)),
        }
        for e in state["catalogs"]
        if getattr(e.get("catalog"), "source_path", None)
    ]
    first_enabled = next(
        (e["path"] for e in catalog_entries if e["enabled"]),
        catalog_entries[0]["path"] if catalog_entries else None,
    )
    # v1.4.0: snapshot the live pointing into the active config, then
    # serialize every live config so the bundle round-trips the whole
    # multi-config plan. Single-config bundles leave `configs` empty.
    _save_active_config_pointing()
    n_cfg = int(state.get("n_configs", 1))
    configs_payload = []
    for ci in range(min(n_cfg, len(state["configs"]))):
        c = state["configs"][ci]
        configs_payload.append({
            "name": c.get("name") or f"Config {ci + 1}",
            "ra_deg": c.get("ra_deg"),
            "dec_deg": c.get("dec_deg"),
            "pa_v3": c.get("pa_v3"),
            "open_shutters": list(c["open_shutters"].values()),
            "highlighted": list(c.get("highlighted") or []),
        })
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
        catalog_path=first_enabled,
        catalog_paths=catalog_entries,
        configs=configs_payload,
        active_config=int(state.get("active_config", 0)),
    )


def on_session_save():
    """Save the current session JSON to the user-supplied path.

    Wrapped in `_confirm_overwrite_if_exists` so the user is asked
    before clobbering an existing session file. v1.3.3+.
    """
    path = session_save_path_input.value.strip()
    if not path:
        _set_status("Set a session save path first.", "warn")
        return
    p = Path(path).expanduser()

    def _do_save():
        p.parent.mkdir(parents=True, exist_ok=True)
        try:
            export_session_json(_build_current_session(), str(p))
            _set_status(f"Session saved → {p}", "ok", clear_after=10)
        except Exception as e:  # noqa: BLE001
            _set_status(f"Session save failed: {e}", "err")
            traceback.print_exc()

    _confirm_overwrite_if_exists(p, _do_save, what="session JSON")


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

    # Loading a session can take several seconds (image decode + catalog
    # loads + overlay rebuild), so always show the spinner immediately,
    # then do the heavy work on the next tick once the overlay has painted.
    _show_loading("Loading session…")
    _close_import_modal()  # session load lives in the Import dialog
    _deferred(_apply_loaded_session, sess)


def _apply_loaded_session(sess) -> None:
    # If we trigger an image load below (by changing a path input), that
    # loader's own `finally` hides the spinner when its slow decode ends —
    # so we leave it up in that case; otherwise we hide it ourselves at
    # the end. (`_show_loading`'s 60 s safety timeout is the backstop.)
    triggered_image = False

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

    # ── Rebuild the config list (v1.4.0) ───────────────────────────────
    # A multi-config bundle carries `sess.configs`; a single-config /
    # legacy bundle describes one config via the top-level pointing +
    # open_shutters. Either way we reset state["configs"] cleanly so a
    # single-config load after a multi-config session clears stale picks.
    sess_configs = getattr(sess, "configs", None)
    if sess_configs and len(sess_configs) > 1:
        new_configs = []
        for ci, entry in enumerate(sess_configs):
            new_configs.append({
                "name": entry.get("name") or f"Config {ci + 1}",
                "open_shutters": {(sh.q, sh.s, sh.d): sh
                                  for sh in entry.get("open_shutters", [])},
                "highlighted": set(entry.get("highlighted") or []),
                "history": [],
                "ra_deg": entry.get("ra_deg"),
                "dec_deg": entry.get("dec_deg"),
                "pa_v3": entry.get("pa_v3"),
            })
        state["configs"] = new_configs
        state["n_configs"] = len(new_configs)
        state["active_config"] = max(
            0, min(int(getattr(sess, "active_config", 0) or 0),
                   len(new_configs) - 1))
    else:
        c0 = _new_config("Config 1")
        c0["open_shutters"] = {(sh.q, sh.s, sh.d): sh
                               for sh in sess.open_shutters}
        c0["highlighted"] = set(sess.highlighted)
        c0["ra_deg"] = float(sess.pointing_ra_deg)
        c0["dec_deg"] = float(sess.pointing_dec_deg)
        c0["pa_v3"] = float(sess.pa_v3_deg)
        state["configs"] = [c0]
        state["n_configs"] = 1
        state["active_config"] = 0
    # Re-point the legacy aliases at the active config.
    _ac = state["configs"][state["active_config"]]
    state["open_shutters"] = _ac["open_shutters"]
    state["highlighted"] = _ac["highlighted"]
    state["history"] = _ac["history"]
    _refresh_config_select_options()
    _prefs_save_suppress["flag"] = True
    try:
        mpt_num_configs_spinner.value = int(state["n_configs"])
    finally:
        _prefs_save_suppress["flag"] = False
    # Load the active config's saved pointing into the widgets (overrides
    # the top-level pointing applied above for a multi-config bundle).
    if _ac["ra_deg"] is not None:
        state["ra_deg"] = float(_ac["ra_deg"])
        ra_input.value = f"{float(_ac['ra_deg']):.6f}"
    if _ac["dec_deg"] is not None:
        state["dec_deg"] = float(_ac["dec_deg"])
        dec_input.value = f"{float(_ac['dec_deg']):.6f}"
    if _ac["pa_v3"] is not None:
        _sync_pa_widgets(float(_ac["pa_v3"]))

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
                _new_jpg = str(img_path)
                if jpg_path_input.value != _new_jpg:
                    jpg_path_input.value = _new_jpg  # triggers JPG load + spinner
                    triggered_image = True
        elif ext in (".fits", ".fit", ".fts"):
            _new_fits = str(img_path)
            if fits_path_input.value != _new_fits:
                fits_path_input.value = _new_fits  # triggers FITS load + spinner
                triggered_image = True
        else:
            image_note = f" Unknown image extension {ext!r}; load manually."
    elif sess.image_path:
        image_note = f" Image not found at {sess.image_path}; load manually."

    # Restore catalogs. New bundles (vMPT 1.1+) record an ordered list
    # in `catalog_paths`; older bundles store a single `catalog_path`,
    # which session_io.py normalises into a one-entry list. We replace
    # the in-memory catalog list with the session's entries so a load
    # is a clean reset, not an additive merge.
    state["catalogs"] = []
    catalog_notes: list[str] = []
    for entry in (sess.catalog_paths or []):
        path = entry.get("path")
        if not path:
            continue
        if not Path(path).exists():
            catalog_notes.append(f"missing {Path(path).name}")
            continue
        try:
            cat = load_catalog(path)
        except Exception as e:  # noqa: BLE001
            catalog_notes.append(f"{Path(path).name}: {e}")
            continue
        state["catalogs"].append({
            "name": Path(path).name,
            "catalog": cat,
            "enabled": bool(entry.get("enabled", True)),
            "color": _assign_catalog_color(len(state["catalogs"])),
        })
    _rebuild_merged_catalog()
    _render_catalog_list()
    if state["catalogs"]:
        # Surface the first loaded catalog's path in the input so the
        # user has something to edit / re-add. Pick the first enabled
        # entry, falling back to the first if all are disabled.
        first = next(
            (e for e in state["catalogs"] if e["enabled"]),
            state["catalogs"][0],
        )
        catalog_path_input.value = first["catalog"].source_path

    refresh_overlays()
    if state["image"] is None and not image_note:
        image_note = " Load an image to see the overlay."
    catalog_note = (
        " Catalog issues: " + "; ".join(catalog_notes) if catalog_notes else ""
    )
    _set_status(
        f"Session loaded: {len(state['open_shutters'])} open shutters."
        f"{image_note}{catalog_note}",
        "warn" if (image_note or catalog_note) else "ok", clear_after=14,
    )
    # An in-flight image load hides the spinner itself; otherwise hide now.
    if not triggered_image:
        _hide_loading()


session_save_btn.on_click(on_session_save)
session_load_btn.on_click(on_session_load)


# ---------------------------------------------------------------------------
# Example data quick-load
# ---------------------------------------------------------------------------


# Locations searched (in priority order) for the example_a370/ and
# example_r0600/ folders. Built lazily on every click so the user can
# `vmpt examples download` AFTER opening the app and the buttons pick
# them up without restarting.
_EX_DEV_DIR = Path(__file__).resolve().parent.parent  # source checkout layout
_EX_USER_DIR = Path.home() / ".vmpt" / "examples"     # standard cache location
                                                       # (matches vmpt cli default)


def _find_example_root() -> Path | None:
    """Return the first directory under which ``example_a370/`` or
    ``example_r0600/`` is found, or ``None`` if neither location has
    them. Order: ``~/.vmpt/examples`` → source-checkout parent → CWD.
    """
    candidates = [
        _EX_USER_DIR,
        _EX_DEV_DIR,
        Path.cwd(),
    ]
    for d in candidates:
        if (d / "example_a370").exists() or (d / "example_r0600").exists():
            return d
    return None


def _example_missing_msg(which: str) -> str:
    """Format a helpful status line when an example folder isn't found,
    pointing the user at the CLI command they need to run."""
    return (
        f"Example {which} not found. Run `vmpt examples download` in "
        f"your terminal to fetch it into {_EX_USER_DIR}, then click "
        "the button again (no restart needed)."
    )


def _reset_pointing_inputs() -> None:
    """Clear RA/Dec inputs so the next image-load auto-recenters."""
    ra_input.value = ""
    dec_input.value = ""


def on_example_a370():
    root = _find_example_root()
    p = (root / "example_a370" / "a370_f182m_f200w_f210m.fits") if root else None
    if p is None or not p.exists():
        _set_status(_example_missing_msg("Abell 370"), "err")
        return
    # Example buttons are a hard reset: clear pointing so the new image
    # auto-recenters (otherwise loading R0600 after A370 leaves pointing
    # at A370 and the MSA disappears off-screen).
    _reset_pointing_inputs()
    fits_path_input.value = str(p)


def on_example_r0600():
    root = _find_example_root()
    jpg = (root / "example_r0600" / "JWST_F090W_F200W_F444W.jpg") if root else None
    wcs = (root / "example_r0600" / "wcs.fits") if root else None
    if jpg is None or not (jpg.exists() and wcs.exists()):
        _set_status(_example_missing_msg("RXCJ0600"), "err")
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
    _close_import_modal()  # plan loads (JSON / .aptx) come from the Import dialog
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
    _set_open_shutters({
        (sh.q, sh.s, sh.d): sh for sh in plan.to_open_shutters()
    })
    n_open = len(state["open_shutters"])
    img: LoadedImage | None = state["image"]
    if img is None:
        _set_status(
            f"Loaded plan '{plan.name}': {n_open} open shutters at "
            f"APA={plan.aperture_pa_deg:.2f}°, V3 PA={plan.v3_pa_deg:.2f}°. "
            f"Load an image (Input tab) to see and edit the overlay.",
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
    _close_import_modal()
    _push_history()
    _set_open_shutters({(sh.q, sh.s, sh.d): sh for sh in opens})
    # Auto-match the loaded mask against the source catalog so the MPT
    # catalog populates without a pre-run of the optimizer. (Rebuild the
    # footprint index first in case the catalog changed since the last
    # pointing update; the rebuild also performs the re-tag, but we call it
    # again to capture the match count for the status line.)
    _rebuild_shutter_catalog_index()
    matched = _retag_manual_opens_from_catalog()
    refresh_overlays()
    _mpt_view_refresh()
    cat = state.get("catalog")
    if cat is not None and len(getattr(cat, "ra_deg", [])):
        extra = (f" — {matched} matched to catalog sources "
                 f"(MPT catalog populated)")
    else:
        extra = " — load a catalog to populate the MPT catalog"
    _set_status(
        f"Loaded shutter CSV: {len(opens)} open shutters{extra}.",
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
_bind_browse(cat_edit_csv_browse_btn, cat_edit_csv_path_input,
             "Save edited catalog as…", _CSV_TYPES, mode="save")
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
# Spec-overlap base alpha per category. Each per-polygon alpha is
# `min(1, base × n_conflicts)`, where n_conflicts is the total
# count of dispersing sources overlapping that shutter. The base
# values live on `state` so the Settings → Overlay appearance
# slider can rewrite them and trigger a re-render via
# `refresh_overlays()`.
state["overlap_base_alpha_stuck"] = 0.20
state["overlap_base_alpha_user"] = 0.20
state["overlap_base_alpha_both"] = 0.20


def _set_overlap_base_alpha(category: str, value: float) -> None:
    """Update the base alpha for one spec-overlap colour and refresh
    the canvas so the per-polygon alphas pick up the new base."""
    state[f"overlap_base_alpha_{category}"] = float(value)
    try:
        refresh_overlays()
    except Exception:  # noqa: BLE001
        pass


_OVERLAY_LAYER_CONFIG = {
    "Operable shutters":     {
        "glyph": bg_shutters_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
        "default_alpha": 0.20, "default_stroke": 1.0,
    },
    # APT MPT-style overlap layers (v1.3.1+). Each carries a
    # field-referenced fill_alpha (`fill_alpha="fill_alpha"` on the
    # glyph), so the per-layer alpha slider can't write directly to
    # the glyph attribute — it goes through `_set_overlap_base_alpha`
    # which stashes the base on state and triggers refresh_overlays
    # to recompute the per-polygon alphas.
    "Mask Stuck (pink)":     {
        "glyph": spec_overlap_stuck_glyph,
        "alpha_attr": None,
        "alpha_state_key": "overlap_base_alpha_stuck",
        "alpha_setter": lambda v: _set_overlap_base_alpha("stuck", v),
        "default_alpha": 0.20, "default_stroke": 0.5,
        "stroke_attr": "line_width",
        "alpha_range": (0.0, 0.8, 0.02),
        "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
        "stroke_extra": lambda v: setattr(
            spec_overlap_stuck_glyph.glyph, "line_alpha",
            0.6 if v > 0 else 0.0,
        ),
    },
    "Masked (overlapping warning)": {
        "glyph": spec_overlap_user_glyph,
        "alpha_attr": None,
        "alpha_state_key": "overlap_base_alpha_user",
        "alpha_setter": lambda v: _set_overlap_base_alpha("user", v),
        "default_alpha": 0.20, "default_stroke": 0.5,
        "stroke_attr": "line_width",
        "alpha_range": (0.0, 0.8, 0.02),
        "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
        "stroke_extra": lambda v: setattr(
            spec_overlap_user_glyph.glyph, "line_alpha",
            0.6 if v > 0 else 0.0,
        ),
    },
    "Mask Conflict (purple)":{
        "glyph": spec_overlap_both_glyph,
        "alpha_attr": None,
        "alpha_state_key": "overlap_base_alpha_both",
        "alpha_setter": lambda v: _set_overlap_base_alpha("both", v),
        "default_alpha": 0.20, "default_stroke": 0.5,
        "stroke_attr": "line_width",
        "alpha_range": (0.0, 0.8, 0.02),
        "stroke_range": (0.0, 3.0, 0.05),
        "stroke_label": "Stroke (px)",
        "stroke_extra": lambda v: setattr(
            spec_overlap_both_glyph.glyph, "line_alpha",
            0.6 if v > 0 else 0.0,
        ),
    },
    "Picked shutters":       {
        "glyph": open_shutters_glyph,
        "alpha_attr": "fill_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 4.0, 0.1),
        "stroke_label": "Stroke (px)",
        "default_alpha": 0.35, "default_stroke": 1.5,
    },
    "Stuck open":            {
        "glyph": stuck_open_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "line_width",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (0.0, 5.0, 0.1),
        "stroke_label": "Stroke (px)",
        "default_alpha": 1.0, "default_stroke": 2.5,
    },
    "Catalog sources":       {
        "glyph": target_glyph,
        "alpha_attr": "line_alpha", "stroke_attr": "size",
        "alpha_range": (0.0, 1.0, 0.05), "stroke_range": (4.0, 30.0, 1.0),
        "stroke_label": "Marker size (px)",
        # line_alpha is field-driven (per-catalog z-depth decay); reset
        # restores the field reference rather than a scalar. size is a
        # plain scalar (the marker diameter in px).
        "default_alpha": "line_alpha", "default_stroke": 10,
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
        # Read the current alpha value either from the glyph (scalar
        # property) or from `state` (when the layer uses a field-
        # referenced fill_alpha — see the MPT-style overlap layers).
        if cfg.get("alpha_attr") is not None:
            cur_alpha = getattr(cfg["glyph"].glyph, cfg["alpha_attr"])
        else:
            cur_alpha = state.get(cfg["alpha_state_key"], a_lo)
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
    # Layers with `alpha_attr=None` use a field-referenced fill_alpha;
    # the slider's value goes to `alpha_setter` (which updates state
    # and re-renders) instead of writing the glyph attribute directly.
    if cfg.get("alpha_attr") is None:
        setter = cfg.get("alpha_setter")
        if setter is not None:
            setter(float(new))
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


# How long to keep the "Resizing canvas…" overlay up after the
# Python-side resize completes. The Python work (mutating
# fig.frame_width / x_range / y_range / tick formatters) takes ~ms,
# but the BROWSER then needs to actually re-render the image — for
# an 8000x12000 JPG this can take a noticeable fraction of a second.
# Without a buffer the spinner blinks away before the user sees the
# resize take effect; with too much buffer it lingers awkwardly.
# 1.2 s covers a typical R2211-sized redraw with margin to spare.
_CANVAS_RESIZE_OVERLAY_MS = 1200


def _canvas_resize_apply(axis: str, value: int) -> None:
    """Apply a canvas-axis resize and schedule the loading overlay
    to fade out after a short buffer.

    Split out so :func:`_on_canvas_x` / :func:`_on_canvas_y` can
    show the loading overlay BEFORE deferring the actual work to a
    next-tick callback — Bokeh batches every Python callback's
    model changes into one browser update, so without the deferral
    the show→resize→hide sequence collapses into a single frame
    and the spinner never appears.
    """
    if axis == "x":
        state["frame_x"] = value
    elif axis == "y":
        state["frame_y"] = value
    if state.get("image") is not None:
        refresh_image_glyph()
    # Hide on a TIMED callback so the spinner stays up while the
    # browser actually paints the resized image + catalog overlay
    # (the heavy lifting happens client-side after the model
    # changes are flushed). add_timeout_callback fires after the
    # given delay in ms.
    try:
        curdoc().add_timeout_callback(
            _hide_loading, _CANVAS_RESIZE_OVERLAY_MS,
        )
    except Exception:  # noqa: BLE001
        _hide_loading()


def _on_canvas_x(attr, old, new):
    """Resize the figure frame's WIDTH in response to the Settings
    slider. Fires only on slider RELEASE (value_throttled) so a
    long drag doesn't trigger a redraw on every mouse tick.

    Shows the same full-page loading overlay (gold spinner on a
    blurred backdrop) that file loads use, so the visual is
    consistent. The overlay stays up for
    `_CANVAS_RESIZE_OVERLAY_MS` after the Python-side resize so the
    browser has time to actually paint the new image at the new
    dimensions before the spinner fades away.
    """
    try:
        v = int(new)
    except (TypeError, ValueError):
        return
    v = max(400, min(1600, v))
    _show_loading("Resizing canvas…")
    curdoc().add_next_tick_callback(
        lambda v=v: _canvas_resize_apply("x", v)
    )


def _on_canvas_y(attr, old, new):
    """Same as :func:`_on_canvas_x` for the Y (height) axis."""
    try:
        v = int(new)
    except (TypeError, ValueError):
        return
    v = max(400, min(1600, v))
    _show_loading("Resizing canvas…")
    curdoc().add_next_tick_callback(
        lambda v=v: _canvas_resize_apply("y", v)
    )


# Spinner commits its `value` only on blur / Enter / arrow click —
# that's naturally throttled, so no `value_throttled` distinction
# is needed (Spinner doesn't expose one).
canvas_x_spinner.on_change("value", _on_canvas_x)
canvas_y_spinner.on_change("value", _on_canvas_y)


# ---------------------------------------------------------------------------
# Drag-pointing handle
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
# Pointing-optimizer callbacks
# ---------------------------------------------------------------------------
# Wrapped imports — heavy (scipy.interpolate Delaunay), so we keep them
# local to delay the first-build cost until the user actually clicks
# Run optimization.
def _open_advanced_modal():
    opt_advanced_modal_backdrop.visible = True
    opt_advanced_modal_card.visible = True


def _close_advanced_modal():
    opt_advanced_modal_backdrop.visible = False
    opt_advanced_modal_card.visible = False


opt_advanced_btn.on_click(_open_advanced_modal)
opt_advanced_modal_close_btn.on_click(_close_advanced_modal)
opt_advanced_modal_top_close_btn.on_click(_close_advanced_modal)


def _refresh_opt_status_div() -> None:
    """Update the status line under the Run button to match the
    current catalog + method state.

    Replaces the static "Load a catalog with priorities, then click
    Run." message that used to show even when a catalog *with*
    priorities was loaded. Triggers: catalog load / remove, method
    change, weight/priority computation, protect-mode toggle.
    """
    cat = state.get("catalog") if state else None
    if cat is None or len(cat.ra_deg) == 0:
        opt_status_div.text = (
            "<small style='color:#5a6b85'>"
            "Load a catalog (Input tab) before running.</small>"
        )
        return
    n = len(cat.ra_deg)
    method = opt_method_select.value or "Democracy"
    pri = np.asarray(getattr(cat, "priority", []), dtype=float)
    wgt = np.asarray(getattr(cat, "weight", []), dtype=float)
    has_pri = pri.size == n and np.isfinite(pri).any()
    has_wgt = wgt.size == n and np.isfinite(wgt).any()
    if method == "Meritocracy" and not has_wgt:
        opt_status_div.text = (
            "<small style='color:#a05a30'>"
            "⚠ Meritocracy needs a <code>weight</code> column. "
            "Use the catalog editor to add weights or "
            "<b>Compute w from p</b>.</small>"
        )
        return
    if method == "Hierarchy" and not has_pri:
        opt_status_div.text = (
            "<small style='color:#a05a30'>"
            "⚠ Hierarchy needs a <code>priority</code> column. "
            "Use the catalog editor to add priorities.</small>"
        )
        return
    opt_status_div.text = (
        f"<small style='color:#1a3b66'>"
        f"Ready · <b>{n:,}</b> sources · {method}.</small>"
    )


# ── Optimizer config modal handlers ─────────────────────────────────
def _open_opt_config_modal():
    # Refresh the status text every time the modal opens so it
    # reflects the catalog + method state the user is about to see.
    _refresh_opt_status_div()
    # Same idea for the per-target centration-override hint — the
    # catalog may have changed since the modal was last opened (editor
    # Apply, session reload, multi-catalog Add/Remove).
    try:
        _refresh_centration_override_hint()
    except NameError:
        pass
    opt_config_modal_backdrop.visible = True
    opt_config_modal_card.visible = True


def _close_opt_config_modal():
    opt_config_modal_backdrop.visible = False
    opt_config_modal_card.visible = False


opt_open_btn.on_click(_open_opt_config_modal)
opt_config_close_btn.on_click(_close_opt_config_modal)
opt_config_top_close_btn.on_click(_close_opt_config_modal)


# ── MPT configurations — switching, count, and the read-only viewer ──────
_mpt_cfg_suppress = {"flag": False}


def _config_label(i: int) -> str:
    return f"Config {i + 1}"


def _save_active_config_pointing() -> None:
    """Snapshot the live widget pointing into the active config dict so it
    is restored when the user switches back to this config."""
    cfg = _active_config()
    cfg["ra_deg"] = state.get("ra_deg")
    cfg["dec_deg"] = state.get("dec_deg")
    cfg["pa_v3"] = state.get("pa_v3")


def _refresh_active_config_banner() -> None:
    """Show/hide the 'Editing Config k of N' banner (multi-config only)."""
    n = int(state.get("n_configs", 1))
    if n <= 1:
        mpt_active_config_div.text = ""
        mpt_active_config_div.visible = False
        return
    act = int(state.get("active_config", 0)) + 1
    # Distinct accent per config so the banner colour matches the active grid.
    accent = _config_color(act - 1)
    mpt_active_config_div.text = (
        f"<div style='background:#eef4ff; border-left:5px solid {accent}; "
        f"border-radius:4px; padding:5px 9px; font-size:12px; color:#23324d'>"
        f"✏️ Editing <b>Config {act}</b> of {n} — manual shutter opens land "
        f"here. Other configs show as faint dashed outlines.</div>"
    )
    mpt_active_config_div.visible = True


def _refresh_config_select_options() -> None:
    """Rebuild the 'Working on' dropdown to list exactly the live configs."""
    n = max(1, int(state.get("n_configs", 1)))
    opts = [_config_label(i) for i in range(n)]
    cur = _config_label(state["active_config"])
    if cur not in opts:
        cur = opts[0]
    _mpt_cfg_suppress["flag"] = True
    try:
        mpt_config_select.options = opts
        mpt_config_select.value = cur
    finally:
        _mpt_cfg_suppress["flag"] = False
    _refresh_active_config_banner()


def _switch_active_config(idx: int) -> None:
    """Make config `idx` the active one: rebind the open_shutters /
    highlighted / history aliases and swap the pointing widgets."""
    idx = max(0, min(int(idx), len(state["configs"]) - 1,
                     int(state["n_configs"]) - 1))
    if idx == state["active_config"]:
        return
    _save_active_config_pointing()
    state["active_config"] = idx
    cfg = _active_config()
    state["open_shutters"] = cfg["open_shutters"]
    state["highlighted"] = cfg["highlighted"]
    state["history"] = cfg["history"]
    # Load this config's saved pointing (None = keep current widgets).
    if cfg["ra_deg"] is not None and cfg["dec_deg"] is not None:
        try:
            ra_input.value = f'{float(cfg["ra_deg"]):.6f}'
            dec_input.value = f'{float(cfg["dec_deg"]):.6f}'
        except (TypeError, ValueError):
            pass
        if cfg["pa_v3"] is not None:
            try:
                _sync_pa_widgets(float(cfg["pa_v3"]))
            except (TypeError, ValueError):
                pass
    _refresh_config_select_options()
    _rebuild_shutter_catalog_index()
    refresh_overlays()
    _set_status(
        f"Now editing {cfg['name']} — "
        f"{len(cfg['open_shutters'])} open shutter(s).",
        "ok", clear_after=8,
    )


def _ensure_n_configs(n: int) -> None:
    """Grow/cap the live config list to `n` (1.._MAX_CONFIGS).

    New configs are born with pointing = None ("inherit the live pointing
    when first activated"), NOT a creation-time snapshot — otherwise a
    config created before an image loads (e.g. via a persisted
    default_num_configs pref) would freeze at (0, 0) and switching to it
    would jump the MSA off the image. `_switch_active_config` keeps the
    current widgets when the target config's pointing is None, so a
    never-positioned config simply shares wherever you are now."""
    n = max(1, min(_MAX_CONFIGS, int(n)))
    while len(state["configs"]) < n:
        state["configs"].append(
            _new_config(_config_label(len(state["configs"])))
        )
    state["n_configs"] = n
    if state["active_config"] >= n:
        _switch_active_config(0)
    _refresh_config_select_options()


def _on_mpt_config_select(attr, old, new):
    if _mpt_cfg_suppress["flag"]:
        return
    try:
        idx = int(str(new).split()[-1]) - 1
    except (ValueError, IndexError):
        idx = 0
    _switch_active_config(idx)


def _on_mpt_num_configs(attr, old, new):
    try:
        n = int(new)
    except (TypeError, ValueError):
        return
    _ensure_n_configs(n)
    # Repaint so the top-bar CONFIG chip + idle outlines reflect the count.
    try:
        refresh_overlays()
    except Exception:  # noqa: BLE001
        pass
    if int(state["n_configs"]) > 1:
        _set_status(
            f"{int(state['n_configs'])} MPT configs active — pick which to "
            f"work on with the 'Working on' dropdown.",
            "ok", clear_after=10,
        )


def _fmt_gap(glo, ghi, has_spectrum: bool) -> str:
    """Label for the viewer's Gap column: ``"lo–hi"`` (μm) when the
    NRS1/NRS2 detector gap falls inside the spectrum, ``"none"`` when the
    spectrum is gap-free, ``""`` when there's no spectrum at all."""
    try:
        if glo is not None and ghi is not None:
            return f"{float(glo):.3f}–{float(ghi):.3f}"
    except (TypeError, ValueError):
        return ""
    return "none" if has_spectrum else ""


def _shutter_lambda_cov(q: int, s: int, d: int) -> tuple[float, float, str]:
    """(λ_blue, λ_red) in μm (NaN where no spectrum) and a Gap-range label
    for shutter (q, s, d) under the current Disperser / Filter, via the
    accurate per-shutter ``wavelengths.cutoffs`` table path."""
    try:
        v2 = float(V2_MSA[q - 1, s - 1, d - 1])
        v3 = float(V3_MSA[q - 1, s - 1, d - 1])
        cut = cutoffs(v2, v3, state["disperser"], state["filter"],
                      q=q, s=s, d=d)
    except Exception:  # noqa: BLE001
        return float("nan"), float("nan"), ""
    b, r = cut.get("lam_blue"), cut.get("lam_red")
    glo, ghi = cut.get("lam_gap_lo"), cut.get("lam_gap_hi")
    has_spec = (b is not None) or (r is not None)
    return (b if b is not None else float("nan"),
            r if r is not None else float("nan"),
            _fmt_gap(glo, ghi, has_spec))


def _mpt_view_collect_rows() -> list:
    """One dict per (config, selected source). A source counted once per
    config even if it spans several shutters (reported at the shutter it
    sits in). Catalog fields filled from the merged active catalog.
    Each row also carries the shutter's wavelength coverage (λ_blue,
    λ_red) + detector-gap range under the current Disperser / Filter."""
    cat = state["catalog"]
    id_to_row: dict = {}
    if cat is not None:
        for i in range(len(cat.ra_deg)):
            id_to_row[str(cat.ids[i])] = i
    rows: list = []
    n = min(int(state.get("n_configs", 1)), len(state["configs"]))
    for ci in range(n):
        cfg = state["configs"][ci]
        seen: dict = {}
        for key in sorted(cfg["open_shutters"].keys()):
            sh = cfg["open_shutters"][key]
            for sid in _open_shutter_ids(sh):
                # Report each source at the shutter where it is the slitlet's
                # PRIMARY (centred) target if such a shutter is open; only
                # fall back to a flanking/neighbour shutter otherwise. This
                # keeps the reported (Q, s, d) the source's real target
                # shutter rather than an arbitrary first-sorted one.
                is_primary = str(getattr(sh, "target_id", "")) == str(sid)
                if sid not in seen or (is_primary and not seen[sid][3]):
                    seen[sid] = (key[0], key[1], key[2], is_primary)
        for sid, (q, s, d, _is_primary) in seen.items():
            r = id_to_row.get(sid)
            row = {"config": ci + 1, "id": sid, "q": q, "s": s, "d": d}
            if r is not None:
                lbl = ""
                try:
                    lbl = str(cat.label[r]).strip()
                except (AttributeError, IndexError, TypeError):
                    lbl = ""
                row.update({
                    "ra": float(cat.ra_deg[r]),
                    "dec": float(cat.dec_deg[r]),
                    "priority": (cat.priority[r]
                                 if r < len(cat.priority) else float("nan")),
                    "weight": (cat.weight[r]
                               if r < len(cat.weight) else float("nan")),
                    "mag": (cat.mag[r]
                            if r < len(cat.mag) else float("nan")),
                    "z": cat.z[r] if r < len(cat.z) else float("nan"),
                    "label": lbl or sid,
                })
            else:
                row.update({"ra": float("nan"), "dec": float("nan"),
                            "priority": float("nan"), "weight": float("nan"),
                            "mag": float("nan"),
                            "z": float("nan"), "label": "(not in catalog)"})
            lam_b, lam_r, gap_str = _shutter_lambda_cov(q, s, d)
            row["lam_blue"], row["lam_red"], row["gap"] = lam_b, lam_r, gap_str
            rows.append(row)
    return rows


def _mpt_view_refresh() -> None:
    rows = _mpt_view_collect_rows()

    # Numeric columns are stored as floats (NaN = missing) so the
    # DataTable sorts them numerically; the column formatters render
    # them (int / fixed decimals) and blank NaN. config / id / label
    # stay strings.
    def _f(v) -> float:
        try:
            f = float(v)
            return f if np.isfinite(f) else float("nan")
        except (TypeError, ValueError):
            return float("nan")

    data = dict(config=[], id=[], ra=[], dec=[], priority=[], weight=[],
                mag=[], z=[], lam_blue=[], lam_red=[], gap=[],
                q=[], s=[], d=[], label=[])
    for r in rows:
        data["config"].append(str(r["config"]))
        data["id"].append(str(r["id"]))
        data["ra"].append(_f(r.get("ra")))
        data["dec"].append(_f(r.get("dec")))
        data["priority"].append(_f(r.get("priority")))
        data["weight"].append(_f(r.get("weight")))
        data["mag"].append(_f(r.get("mag")))
        data["z"].append(_f(r.get("z")))
        data["lam_blue"].append(_f(r.get("lam_blue")))
        data["lam_red"].append(_f(r.get("lam_red")))
        data["gap"].append(str(r.get("gap", "")))
        data["q"].append(_f(r.get("q")))
        data["s"].append(_f(r.get("s")))
        data["d"].append(_f(r.get("d")))
        data["label"].append(str(r.get("label", "")))
    _mpt_view_source.data = data

    if not rows:
        mpt_view_summary_div.text = (
            "<span style='color:#5a6b85'>No sources selected yet. Load a "
            "shutter mask, open shutters by hand, or run the optimizer — "
            "any of these populates this list (with a catalog loaded).</span>"
        )
        return
    per_cfg: dict = {}
    for r in rows:
        per_cfg[r["config"]] = per_cfg.get(r["config"], 0) + 1
    parts = ", ".join(f"Config&nbsp;{c}:&nbsp;{per_cfg[c]}"
                      for c in sorted(per_cfg))
    uniq = len({r["id"] for r in rows})
    mpt_view_summary_div.text = (
        f"<b>{len(rows)}</b> selected row(s) across "
        f"<b>{min(int(state['n_configs']), len(state['configs']))}</b> "
        f"config(s) — {parts}. <b>{uniq}</b> unique source(s)."
    )


def _open_mpt_view():
    # Computing per-shutter wavelength coverage for every selected source
    # can take a moment; show the spinner first, then build on the next
    # tick so the overlay paints before the (blocking) compute + render.
    _show_loading("Computing wavelength coverage…")

    def _build():
        try:
            _mpt_view_rebuild_columns()
            _mpt_view_refresh()
            mpt_view_modal_backdrop.visible = True
            mpt_view_modal_card.visible = True
        finally:
            _hide_loading()
    _deferred(_build)


def _mpt_view_close():
    mpt_view_modal_backdrop.visible = False
    mpt_view_modal_card.visible = False


def _mpt_view_save_csv():
    path = mpt_view_csv_path_input.value.strip()
    if not path:
        _set_status("Set a CSV path for the MPT catalog first.", "warn")
        return
    rows = _mpt_view_collect_rows()
    if not rows:
        _set_status("Nothing selected to save.", "warn")
        return
    p = Path(path).expanduser()

    def _do_save():
        import csv

        def g(r, k, nd=6, as_int=False) -> str:
            v = r.get(k)
            try:
                f = float(v)
                if not np.isfinite(f):
                    return ""
                if as_int:
                    return str(int(round(f)))
                return (f"{f:.{nd}f}".rstrip("0").rstrip(".")) or "0"
            except (TypeError, ValueError):
                return "" if v is None else str(v)

        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["config", "ID", "RA", "DEC", "priority", "weight",
                        "z", "lam_blue", "lam_red", "gap",
                        "q", "s", "d", "label"])
            for r in rows:
                w.writerow([
                    r["config"], r["id"], g(r, "ra"), g(r, "dec"),
                    g(r, "priority", as_int=True), g(r, "weight", as_int=True),
                    g(r, "z", 4), g(r, "lam_blue", 4), g(r, "lam_red", 4),
                    r.get("gap", ""),
                    r["q"], r["s"], r["d"],
                    r.get("label", ""),
                ])
        _set_status(f"Saved {len(rows)} MPT-catalog rows → {p}", "ok",
                    clear_after=12)

    _confirm_overwrite_if_exists(p, _do_save, what="MPT catalog CSV")


mpt_num_configs_spinner.on_change("value", _on_mpt_num_configs)
mpt_config_select.on_change("value", _on_mpt_config_select)
mpt_view_btn.on_click(_open_mpt_view)
mpt_view_close_btn.on_click(_mpt_view_close)
mpt_view_top_close_btn.on_click(_mpt_view_close)
mpt_view_csv_save_btn.on_click(_mpt_view_save_csv)

# Keep the status line live as the user changes the method or
# touches the catalog. The catalog-change hook is set up later, in
# the catalog-load handlers, where state["catalog"] is mutated.
opt_method_select.on_change(
    "value", lambda attr, old, new: _refresh_opt_status_div(),
)


def _refresh_centration_override_hint() -> None:
    """Update the per-target centration-override hint Div under the
    optimizer modal's Source-centering Select.

    Counts how many catalog rows have a non-empty ``centration`` value;
    if ≥1, shows ``N sources have per-target centering overrides —
    those rows ignore this setting.`` Empty otherwise. Same triggers
    as ``_refresh_opt_status_div`` (catalog load / remove, editor
    Apply, session reload).
    """
    cat = state.get("catalog") if state else None
    if cat is None:
        opt_centration_override_hint.text = ""
        return
    cent = np.asarray(
        getattr(cat, "centration", []), dtype=object,
    )
    n_overrides = sum(1 for v in cent if str(v).strip())
    if n_overrides == 0:
        opt_centration_override_hint.text = ""
        return
    word_source = "source" if n_overrides == 1 else "sources"
    word_verb = "has" if n_overrides == 1 else "have"
    opt_centration_override_hint.text = (
        "<small style='color:#5a6b85; font-style:italic'>"
        f"{n_overrides} {word_source} {word_verb} a per-target "
        "centering override — those rows ignore this setting."
        "</small>"
    )


# ── Catalog editor handlers ──────────────────────────────────────────────

# Undo / redo stacks hold prior-state snapshots (copies of the
# ColumnDataSource data dict). `_cat_edit_suppress` muzzles the
# on_change("data") listener while we're applying snapshots — otherwise
# undo / redo / populate operations would recursively push to the stack.
_cat_edit_undo_stack: list[dict] = []
_cat_edit_redo_stack: list[dict] = []
_cat_edit_suppress: dict = {"flag": False}
_CAT_EDIT_MAX_HISTORY = 100


def _cat_edit_snapshot() -> dict:
    """Take a deep copy of the current working-copy source.data."""
    return {k: list(v) for k, v in _cat_edit_source.data.items()}


def _cat_edit_set_data_silently(data: dict) -> None:
    """Replace source.data without triggering an undo push."""
    _cat_edit_suppress["flag"] = True
    try:
        _cat_edit_source.data = {k: list(v) for k, v in data.items()}
    finally:
        _cat_edit_suppress["flag"] = False


def _cat_edit_render_history():
    n_u = len(_cat_edit_undo_stack)
    n_r = len(_cat_edit_redo_stack)
    cat_edit_undo_btn.disabled = (n_u == 0)
    cat_edit_redo_btn.disabled = (n_r == 0)
    cat_edit_history_div.text = (
        f"<small style='color:#5a6b85'>"
        f"{n_u} undo · {n_r} redo</small>"
    )


def _cat_edit_fmt_optional(v) -> str:
    """Render a possibly-NaN numeric as a clean string ("" for NaN)."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(f):
        return ""
    return _fmt_num(f)


# Standard columns the loader always populates. Extras live under
# `Catalog.extras` and are surfaced through the same picker.
_CAT_STD_COLS = ("id", "ra", "dec", "priority", "weight", "mag", "z", "label")
_CAT_STD_TITLES = {
    "id": "ID", "ra": "RA (deg)", "dec": "Dec (deg)",
    "priority": "Priority", "weight": "Weight",
    "mag": "Mag", "z": "z", "label": "Label",
}
_CAT_STD_WIDTHS = {
    "id": 110, "ra": 110, "dec": 110,
    "priority": 80, "weight": 80,
    "mag": 80, "z": 80, "label": 150,
}


def _cat_edit_populate_table(idx: int) -> None:
    """Copy the indexed catalog's rows into the editor's source.

    All numeric columns are PRE-FORMATTED to strings (NaN → "") so the
    DataTable's StringEditor never sees a NaN it can't render. Extras
    columns from the original CSV/FITS file are stuffed in alongside
    the standard ones — the column picker (`cat_edit_columns_choice`)
    controls which are visible in the table."""
    if not (0 <= idx < len(state["catalogs"])):
        _cat_edit_set_data_silently(dict(
            id=[], ra=[], dec=[], priority=[], weight=[],
            mag=[], z=[], label=[], _idx=[],
            lam_req=[], no_gap=[], extend_blue=[], extend_red=[],
            protect=[], centration=[], max_configs=[], _has_constraint=[],
        ))
        cat_edit_columns_choice.options = []
        cat_edit_columns_choice.value = []
    else:
        cat = state["catalogs"][idx]["catalog"]
        n = len(cat.ra_deg)
        ra_arr = np.asarray(cat.ra_deg, dtype=float)
        dec_arr = np.asarray(cat.dec_deg, dtype=float)
        # `weight` may be missing (older sessions / catalogs loaded
        # before the field existed). Pad to the row count with NaN.
        weight_arr = np.asarray(getattr(cat, "weight", []), dtype=float)
        if weight_arr.size != n:
            weight_arr = np.full(n, np.nan, dtype=float)
        # Per-target constraints (v1.3.0+). All optional — pad with the
        # zero / empty default when the Catalog dataclass field is
        # missing or has a mismatched length (older sessions etc.).
        from vmpt.catalog import _format_lam_req as _fmt_lam_req
        def _pad_bool(arr_name: str) -> list[int]:
            arr = np.asarray(getattr(cat, arr_name, []), dtype=bool)
            if arr.size != n:
                arr = np.zeros(n, dtype=bool)
            return [int(bool(v)) for v in arr]
        no_gap_list = _pad_bool("no_gap")
        extend_blue_list = _pad_bool("extend_blue")
        extend_red_list = _pad_bool("extend_red")
        protect_list = _pad_bool("protect")
        req_lam_arr = getattr(cat, "required_lam", None)
        lam_req_list: list[str] = []
        if req_lam_arr is not None and len(req_lam_arr) == n:
            for entry in req_lam_arr:
                lam_req_list.append(_fmt_lam_req(entry))
        else:
            lam_req_list = ["" for _ in range(n)]
        # Per-target centration override (v1.3.1+) — empty string means
        # "use the optimizer's global Source-centering setting". Pad to
        # length n with "" when the Catalog dataclass field is missing
        # or has the wrong length.
        centration_arr = np.asarray(
            getattr(cat, "centration", []), dtype=object,
        )
        if centration_arr.size != n:
            centration_list = ["" for _ in range(n)]
        else:
            centration_list = [str(v).strip() for v in centration_arr]
        # Per-target multi-config cap (v1.4.0+). Stored in the CDS as a
        # string ("" = unset, "1"/"2") to round-trip with the popover
        # Select; the Catalog field is float (NaN = unset).
        mc_arr = np.asarray(getattr(cat, "max_configs", []), dtype=float)
        if mc_arr.size != n:
            max_configs_list = ["" for _ in range(n)]
        else:
            max_configs_list = [
                ("" if not np.isfinite(v) else str(int(round(v))))
                for v in mc_arr
            ]
        has_constraint_list = [
            int(bool(lam_req_list[i]) or no_gap_list[i] or
                extend_blue_list[i] or extend_red_list[i]
                or protect_list[i] or bool(centration_list[i])
                or bool(max_configs_list[i]))
            for i in range(n)
        ]
        # Priority and weight are stored as FLOATS (NaN for missing)
        # — that makes SlickGrid sort them numerically instead of
        # lexicographically ("10" < "2" with string sort). The
        # HTMLTemplateFormatter in the column definition renders the
        # float as an integer and blanks NaN. StringEditor still
        # works for editing — `_on_cat_edit_data_change` coerces user
        # input back to float.
        data: dict = dict(
            id=[str(v) for v in np.asarray(cat.ids).tolist()],
            ra=[f"{v:.7f}" for v in ra_arr.tolist()],
            dec=[f"{v:.7f}" for v in dec_arr.tolist()],
            priority=[float(v) for v in cat.priority],
            weight=[float(v) for v in weight_arr],
            mag=[_cat_edit_fmt_optional(v) for v in cat.mag],
            # `z` stored as FLOAT (NaN for missing) so SlickGrid
            # sorts the column numerically — same pattern priority +
            # weight follow. The HTMLTemplateFormatter at the table
            # column level renders 3 decimals + blanks NaN; the
            # data-change handler coerces user-typed strings back to
            # float on edit commit. v1.3.2+.
            z=[float(v) if (isinstance(v, (int, float))
                            and np.isfinite(float(v)))
               else float("nan")
               for v in cat.z],
            label=[str(v) for v in (cat.label if cat.label is not None
                                    else [""] * n)],
            _idx=list(range(n)),
            # --- per-target constraints ---
            lam_req=lam_req_list,
            no_gap=no_gap_list,
            extend_blue=extend_blue_list,
            extend_red=extend_red_list,
            protect=protect_list,
            centration=centration_list,
            max_configs=max_configs_list,
            _has_constraint=has_constraint_list,
        )
        # Add every extras column verbatim (already object arrays of
        # str values, courtesy of `load_catalog`).
        extras = getattr(cat, "extras", {}) or {}
        for ex_name, ex_vals in extras.items():
            # Avoid clobbering a standard column if some catalog uses
            # the same name (it shouldn't, since `claimed` in the
            # loader prevents that, but defend in depth).
            if ex_name in data:
                continue
            data[ex_name] = list(ex_vals)
        _cat_edit_set_data_silently(data)

        # Refresh the column picker. Default selection: the 7 standard
        # columns. Extras are listed but off by default so the table
        # opens to its familiar layout.
        std_titled = [(c, _CAT_STD_TITLES[c]) for c in _CAT_STD_COLS]
        ex_titled = [(k, k) for k in extras.keys()]
        all_opts = std_titled + ex_titled
        cat_edit_columns_choice.options = [t for _, t in all_opts]
        cat_edit_columns_choice.value = [t for c, t in std_titled]
    _cat_edit_undo_stack.clear()
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()
    _cat_edit_rebuild_columns()


# Integer-display formatter for priority + weight. Source data is
# stored as FLOAT (NaN for missing) so SlickGrid sorts the column
# numerically. The template renders ints (rounds floats) and blanks
# NaN/null/empty; if a value is still a string (transient state
# right after the user commits an edit, before `_on_cat_edit_data_change`
# coerces it back to float), it's shown verbatim so the user sees
# their typed text.
_INT_OR_BLANK_TEMPLATE = (
    "<%= (value === null || value === undefined || value === '' || "
    "     (typeof value === 'number' && isNaN(value))) "
    "    ? '' "
    "    : (typeof value === 'number' ? Math.round(value) : value) %>"
)
# Like `_INT_OR_BLANK_TEMPLATE` but renders 3 decimal places instead
# of rounding to int — used for the redshift `z` column (v1.3.2+),
# which is stored as float so SlickGrid sorts numerically. A
# transient string-typed value (right after a user commits an edit
# but before the data-change handler coerces it back to float) is
# shown verbatim so the user sees their typed text.
_FLOAT_OR_BLANK_TEMPLATE = (
    "<%= (value === null || value === undefined || value === '' || "
    "     (typeof value === 'number' && isNaN(value))) "
    "    ? '' "
    "    : (typeof value === 'number' ? value.toFixed(3) : value) %>"
)


def _cat_edit_rebuild_columns() -> None:
    """Rebuild `cat_edit_table.columns` according to the picker. The
    Constraints + trash columns always stay at the end. Priority +
    Weight use a NaN-safe integer formatter so the sort is numeric,
    not lexicographic."""
    visible_titles = set(cat_edit_columns_choice.value or [])
    title_to_field = {**{_CAT_STD_TITLES[c]: c for c in _CAT_STD_COLS}}
    data = _cat_edit_source.data
    # Source-data keys that the editor manages internally — never
    # exposed in the column picker, never rendered as a normal column.
    _INTERNAL = (
        "_idx", "_has_constraint",
        # Constraint values surface via the dedicated Constraints…
        # popover, not as inline-editable columns.
        "lam_req", "no_gap", "extend_blue", "extend_red", "protect",
        "centration", "max_configs",
    )
    # Extras: any source-data key not in the standard set + not
    # internal.
    extras = [k for k in data.keys()
              if k not in _CAT_STD_COLS and k not in _INTERNAL]
    for k in extras:
        title_to_field[k] = k

    cols: list = []
    int_fmt = HTMLTemplateFormatter(template=_INT_OR_BLANK_TEMPLATE)
    # 3-decimal float formatter for the redshift column.
    float_fmt = HTMLTemplateFormatter(template=_FLOAT_OR_BLANK_TEMPLATE)
    # Standard columns first, in their fixed order.
    for c in _CAT_STD_COLS:
        if _CAT_STD_TITLES[c] in visible_titles:
            col_kwargs = dict(
                field=c, title=_CAT_STD_TITLES[c],
                editor=StringEditor(), width=_CAT_STD_WIDTHS[c],
            )
            if c in ("priority", "weight"):
                col_kwargs["formatter"] = int_fmt
            elif c == "z":
                col_kwargs["formatter"] = float_fmt
            cols.append(TableColumn(**col_kwargs))
    # Extras in CSV order.
    for k in extras:
        if k in visible_titles:
            cols.append(TableColumn(
                field=k, title=k, editor=StringEditor(), width=120,
            ))
    # Constraints column: always visible, sits just before the trash.
    cols.append(TableColumn(
        field="_idx", title="Constraints", width=82,
        formatter=HTMLTemplateFormatter(template=_CONSTRAINT_TEMPLATE),
        sortable=False,
    ))
    # Trash column always last.
    cols.append(TableColumn(
        field="_idx", title="🗑", width=34,
        formatter=HTMLTemplateFormatter(template=_TRASH_TEMPLATE),
        sortable=False,
    ))
    cat_edit_table.columns = cols


def _on_cat_edit_columns_change(attr, old, new):
    _cat_edit_rebuild_columns()


def _on_cat_edit_add_column():
    """Append a user-named string column to the working copy.

    The new column starts empty and is auto-ticked in the picker so
    it appears in the table immediately. Edits flow into it via
    Bokeh's normal cell-edit path; Apply pushes it back into
    `Catalog.extras` so it round-trips through Save-as-CSV and
    survives session reload."""
    name = (cat_edit_new_col_input.value or "").strip()
    if not name:
        _set_status("Type a column name first.", "warn")
        return
    # Disallow names that would shadow internal or standard fields.
    if name == "_idx":
        _set_status("`_idx` is reserved.", "err")
        return
    data = dict(_cat_edit_source.data)
    n = len(data.get("ra", []))
    if name in data:
        # Already exists — just make sure it's ticked + visible.
        title = _CAT_STD_TITLES.get(name, name)
        if title not in (cat_edit_columns_choice.value or []):
            cat_edit_columns_choice.value = [
                *(cat_edit_columns_choice.value or []), title,
            ]
        _set_status(
            f"Column {name!r} already exists; made it visible.",
            "ok", clear_after=8,
        )
        cat_edit_new_col_input.value = ""
        return
    # Add as an empty-string column.
    data[name] = ["" for _ in range(n)]
    _cat_edit_set_data_silently(data)
    # Push to the picker as a new option and tick it.
    cat_edit_columns_choice.options = [
        *(cat_edit_columns_choice.options or []), name,
    ]
    cat_edit_columns_choice.value = [
        *(cat_edit_columns_choice.value or []), name,
    ]
    _cat_edit_rebuild_columns()
    _set_status(
        f"Added column {name!r}. Fill values then click Apply.",
        "ok", clear_after=10,
    )
    cat_edit_new_col_input.value = ""


from vmpt.catalog_ops import (
    compute_priorities_from_weights as _compute_priorities_from_weights,
    compute_weights_from_priorities as _compute_weights_from_priorities,
)


def _str_or_blank_to_float(v) -> float:
    """Helper: empty/NaN-ish strings → NaN; numeric strings → float."""
    if isinstance(v, (int, float)):
        return float(v)
    s = ("" if v is None else str(v)).strip()
    if not s or s.lower() in ("nan", "none", "null", "--"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _on_cat_edit_compute_w_from_p():
    data = dict(_cat_edit_source.data)
    pris = list(data.get("priority", []))
    new_w_str = _compute_weights_from_priorities(pris)
    if new_w_str is None:
        _set_status(
            "Compute w from p: no finite priorities found. "
            "Fill the Priority column first.", "warn",
        )
        return
    # Helper returns strings (per its signature). Convert to floats so
    # the column stays numerically sortable.
    data["weight"] = [_str_or_blank_to_float(s) for s in new_w_str]
    _cat_edit_set_data_silently(data)
    # Make sure Weight is visible in the picker.
    titles = cat_edit_columns_choice.value or []
    if _CAT_STD_TITLES["weight"] not in titles:
        cat_edit_columns_choice.value = [*titles, _CAT_STD_TITLES["weight"]]
    # Record on undo stack so this isn't a one-way action.
    snap = {k: list(v) for k, v in _cat_edit_source.data.items()}
    _cat_edit_undo_stack.append(snap)
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()
    _set_status("Computed weights from priorities. Click Apply to commit.",
                "ok", clear_after=10)


def _on_cat_edit_compute_p_from_w():
    data = dict(_cat_edit_source.data)
    weights = list(data.get("weight", []))
    new_p_str = _compute_priorities_from_weights(weights)
    if new_p_str is None:
        _set_status(
            "Compute p from w: no finite weights found. "
            "Fill the Weight column first.", "warn",
        )
        return
    data["priority"] = [_str_or_blank_to_float(s) for s in new_p_str]
    _cat_edit_set_data_silently(data)
    titles = cat_edit_columns_choice.value or []
    if _CAT_STD_TITLES["priority"] not in titles:
        cat_edit_columns_choice.value = [*titles, _CAT_STD_TITLES["priority"]]
    snap = {k: list(v) for k, v in _cat_edit_source.data.items()}
    _cat_edit_undo_stack.append(snap)
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()
    _set_status("Computed priorities from weights. Click Apply to commit.",
                "ok", clear_after=10)


def _on_cat_rule_apply():
    """Bulk-set ``max_configs`` for every row matching a boolean column
    expression (the 'condition + value' rule). The expression is validated
    by :func:`evaluate_catalog_condition` FIRST — on any error nothing is
    changed and the message is shown inline. Matching rows get the chosen
    value; non-matching rows keep theirs. Recorded on the undo stack."""
    data = dict(_cat_edit_source.data)
    n = len(data.get("id", []))
    if n == 0:
        cat_rule_status_div.text = (
            "<span style='color:#b00020'>Load a catalog first.</span>")
        return
    expr = (cat_rule_condition_input.value or "").strip()
    # Columns the rule may reference: user-facing data only — skip the
    # editor's internal (_idx, _has_constraint) and constraint-bookkeeping
    # columns so they can't shadow a real column name.
    skip = {"_idx", "_has_constraint", "protect", "centration",
            "max_configs", "lam_req", "no_gap", "extend_blue", "extend_red"}
    cols = {k: list(v) for k, v in data.items()
            if k not in skip and not k.startswith("_")}
    try:
        mask = evaluate_catalog_condition(expr, cols)
    except ValueError as exc:
        msg = (str(exc).replace("&", "&amp;")
               .replace("<", "&lt;").replace(">", "&gt;"))
        cat_rule_status_div.text = f"<span style='color:#b00020'>⚠ {msg}</span>"
        return
    sel = (cat_rule_value_select.value or "").strip()
    set_val = "" if sel == "(use global)" else sel
    mc = list(data.get("max_configs", [""] * n))
    if len(mc) != n:
        mc = [""] * n
    count = 0
    for i in range(n):
        if i < len(mask) and bool(mask[i]):
            mc[i] = set_val
            count += 1
    data["max_configs"] = mc
    # Recompute the per-row "has any constraint" flag (drives the
    # Constraints-column highlight) from the live constraint columns.
    def _col(key):
        return list(data.get(key, [""] * n))
    lam, ng, eb = _col("lam_req"), _col("no_gap"), _col("extend_blue")
    er, pr, ce = _col("extend_red"), _col("protect"), _col("centration")
    data["_has_constraint"] = [
        int(bool(lam[i]) or bool(ng[i]) or bool(eb[i]) or bool(er[i])
            or bool(pr[i]) or bool(ce[i]) or bool(mc[i]))
        for i in range(n)
    ]
    _cat_edit_set_data_silently(data)
    snap = {k: list(v) for k, v in _cat_edit_source.data.items()}
    _cat_edit_undo_stack.append(snap)
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()
    label = "(use global)" if set_val == "" else set_val
    if count == 0:
        cat_rule_status_div.text = (
            "<span style='color:#9a6a00'>No sources matched — nothing "
            "changed.</span>")
    else:
        cat_rule_status_div.text = (
            f"<span style='color:#1a7f37'>✓ Set max configs = {label} "
            f"for {count} source(s). Click <b>Apply changes</b> to commit "
            "to the live catalog.</span>")


def _on_cat_edit_data_change(attr, old, new):
    """Listener on source.data — fires after any cell edit. Pushes the
    PRIOR state onto the undo stack so the user can revert.

    Also coerces priority + weight back to float — they're stored as
    floats so SlickGrid sorts numerically, but StringEditor writes
    a string when the user commits. Without this coercion the column
    becomes mixed-type and sorting flips back to lexicographic."""
    if _cat_edit_suppress["flag"]:
        return
    snap = {k: list(v) for k, v in old.items()}
    _cat_edit_undo_stack.append(snap)
    if len(_cat_edit_undo_stack) > _CAT_EDIT_MAX_HISTORY:
        _cat_edit_undo_stack.pop(0)
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()

    # Coerce any string entries in the numeric columns back to float.
    # priority + weight + z are all float-stored so SlickGrid sorts
    # them numerically; StringEditor commits a string on every cell
    # edit and we coerce here so the column doesn't drift to mixed
    # types (which silently re-enables lexicographic sort).
    #
    # Apply the coercion as a targeted `source.patch()` of just the
    # string-typed cells (the one the user edited, plus any stragglers)
    # rather than reassigning the whole `source.data` — otherwise an
    # ~1800-row table gets re-serialised on every numeric edit. The patch
    # fires a `patching` event that `_cat_edit_scroll_restore_js` (below)
    # listens for to pin the viewport back, so editing a cell doesn't
    # fling the table to the top when a column sort is active.
    patches: dict = {}
    for col in ("priority", "weight", "z"):
        if col not in new:
            continue
        col_patches = []
        for i, v in enumerate(new[col]):
            if not isinstance(v, str):
                continue
            s = v.strip()
            if not s or s.lower() in ("nan", "none", "null", "--"):
                fv = float("nan")
            else:
                try:
                    fv = float(s)
                except ValueError:
                    fv = float("nan")
            col_patches.append((i, fv))
        if col_patches:
            patches[col] = col_patches
    if patches:
        # Suppress recursion: patch() re-triggers this listener.
        _cat_edit_suppress["flag"] = True
        try:
            _cat_edit_source.patch(patches)
        finally:
            _cat_edit_suppress["flag"] = False


def _on_cat_edit_undo():
    if not _cat_edit_undo_stack:
        return
    _cat_edit_redo_stack.append(_cat_edit_snapshot())
    snap = _cat_edit_undo_stack.pop()
    _cat_edit_set_data_silently(snap)
    _cat_edit_render_history()


def _on_cat_edit_redo():
    if not _cat_edit_redo_stack:
        return
    _cat_edit_undo_stack.append(_cat_edit_snapshot())
    snap = _cat_edit_redo_stack.pop()
    _cat_edit_set_data_silently(snap)
    _cat_edit_render_history()


def _on_cat_edit_delete_signal(attr, old, new):
    """JS-side `window.__vmpt_delete_row(idx)` writes the row index here.
    We delete that row from the working copy AND push the prior state
    onto the undo stack."""
    try:
        target = int(new.get("idx", [-1])[0])
    except (TypeError, ValueError, IndexError):
        return
    if target < 0:
        return
    data = _cat_edit_source.data
    idxs = list(data.get("_idx", []))
    if target not in idxs:
        return
    # Capture the prior state for undo.
    snap = {k: list(v) for k, v in data.items()}
    _cat_edit_undo_stack.append(snap)
    if len(_cat_edit_undo_stack) > _CAT_EDIT_MAX_HISTORY:
        _cat_edit_undo_stack.pop(0)
    _cat_edit_redo_stack.clear()

    # Drop the row whose `_idx` == target. (We use _idx, not the
    # positional index in the list, because the user may have sorted.)
    keep = [i for i, v in enumerate(idxs) if v != target]
    new_data = {k: [v[i] for i in keep] for k, v in data.items()}
    # Re-stamp _idx so future rows have unique non-colliding indices.
    new_data["_idx"] = list(range(len(keep)))
    _cat_edit_set_data_silently(new_data)
    _cat_edit_render_history()


def _cat_edit_open():
    if not state["catalogs"]:
        _set_status("Load a catalog first.", "warn")
        return
    options = [
        f"#{i + 1}: {e['name']} ({len(e['catalog'].ra_deg)} rows)"
        for i, e in enumerate(state["catalogs"])
    ]
    cat_edit_select.options = options
    # Default: edit the first enabled catalog (or first overall).
    default_idx = next(
        (i for i, e in enumerate(state["catalogs"]) if e["enabled"]),
        0,
    )
    cat_edit_select.value = options[default_idx]
    _cat_edit_populate_table(default_idx)
    # Pre-fill the save-as path with "<source>_edited.csv" — only if
    # the input box is currently empty.
    if not cat_edit_csv_path_input.value.strip():
        src = state["catalogs"][default_idx]["catalog"].source_path or ""
        if src:
            stem, _, ext = src.rpartition(".")
            suggested = (stem or src) + "_edited.csv"
            cat_edit_csv_path_input.value = suggested
    cat_edit_modal_backdrop.visible = True
    cat_edit_modal_card.visible = True


def _cat_edit_close():
    cat_edit_modal_backdrop.visible = False
    cat_edit_modal_card.visible = False


def _cat_edit_selected_idx() -> int:
    """Return the index of the catalog currently picked in the dropdown."""
    try:
        return int(cat_edit_select.value.split(":", 1)[0].lstrip("#")) - 1
    except (ValueError, AttributeError, IndexError):
        return -1


def _on_cat_edit_select(attr, old, new):
    _cat_edit_populate_table(_cat_edit_selected_idx())


def _str_to_float_or_nan(v) -> float:
    """Tolerantly parse a cell value back to float. Empty / unparsable
    → NaN. Used when committing string-edited columns to numeric
    Catalog arrays."""
    s = str(v).strip()
    if not s or s.lower() in ("nan", "none", "null", "--"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _on_cat_edit_apply():
    """Write the working copy back to state['catalogs'][idx] and refresh."""
    idx = _cat_edit_selected_idx()
    if not (0 <= idx < len(state["catalogs"])):
        _set_status("No catalog selected.", "warn")
        return
    data = _cat_edit_source.data
    n = len(data["ra"])
    if n == 0:
        _set_status(
            "Catalog is empty after edits. Delete the catalog instead "
            "if you don't need it.", "err",
        )
        return
    # Coerce string columns back to the Catalog dataclass schema.
    ids_in = list(data["id"])
    ids_out: list = []
    for v in ids_in:
        s = str(v).strip()
        try:
            ids_out.append(int(float(s)))
        except (ValueError, TypeError):
            ids_out.append(s)
    cat = state["catalogs"][idx]["catalog"]
    # Round-trip extras + any user-added columns. Treat the working
    # copy as the source of truth: every column in source.data that
    # isn't a standard field or the internal _idx becomes / updates
    # an extras entry. Columns in the original `cat.extras` that the
    # user removed from the picker still survive (we never strip
    # source.data on column-hide), but if a user explicitly added a
    # new column it lands here for the first time.
    new_extras: dict = {}
    src_keys = set(_cat_edit_source.data.keys())
    # Editor-internal columns that are NOT extras: standard fields,
    # the row-index helper, and the per-target constraint fields
    # (which are first-class Catalog attributes since v1.3.0).
    _NON_EXTRAS = set(_CAT_STD_COLS) | {
        "_idx", "_has_constraint",
        "lam_req", "no_gap", "extend_blue", "extend_red", "protect",
        "centration",
    }
    for k in src_keys:
        if k in _NON_EXTRAS:
            continue
        new_extras[k] = np.asarray(list(data[k]), dtype=object)
    # Carry through any extras that for some reason weren't in source
    # (defensive: should never happen with the current flow).
    for k, v in (getattr(cat, "extras", {}) or {}).items():
        new_extras.setdefault(k, np.asarray(list(v), dtype=object))
    weight_in = data.get("weight", ["" for _ in range(n)])

    # Per-target constraint columns. The editor stores them as:
    #   lam_req      : list[str]   ("1.0-1.3; 1.5-1.8" or "")
    #   no_gap/etc.  : list[int]   (0 or 1)
    # We parse the lam_req strings back to list[(lo, hi)] and cast
    # the bool fields. Missing-from-source = default (False / []).
    from vmpt.catalog import _parse_lam_req_str
    lam_strs = list(data.get("lam_req", [""] * n))
    while len(lam_strs) < n:
        lam_strs.append("")
    new_required_lam = np.empty(n, dtype=object)
    for i, s in enumerate(lam_strs[:n]):
        new_required_lam[i] = _parse_lam_req_str(s)

    def _bool_col(name: str) -> np.ndarray:
        col = data.get(name, [0] * n)
        return np.asarray([bool(int(v)) if str(v).strip() not in ("", "nan",
                                                                  "none",
                                                                  "null")
                           else False
                           for v in col[:n]], dtype=bool)

    def _centration_col_from_editor(
        editor_data: dict, n_rows: int,
    ) -> np.ndarray:
        """Build a length-N object array of centration labels from the
        editor's working source. Empty / unrecognised cells normalise
        to ``""``. Used by `_cat_edit_apply_to_catalog`."""
        from vmpt.catalog import _normalise_centration
        col = editor_data.get("centration", [""] * n_rows)
        out = np.empty(n_rows, dtype=object)
        for i in range(n_rows):
            v = col[i] if i < len(col) else ""
            out[i] = _normalise_centration(v)
        return out

    new_cat = Catalog(
        ids=np.asarray(ids_out, dtype=object),
        ra_deg=np.asarray([_str_to_float_or_nan(v) for v in data["ra"]],
                          dtype=float),
        dec_deg=np.asarray([_str_to_float_or_nan(v) for v in data["dec"]],
                           dtype=float),
        priority=np.asarray([_str_to_float_or_nan(v) for v in data["priority"]],
                            dtype=float),
        weight=np.asarray([_str_to_float_or_nan(v) for v in weight_in],
                          dtype=float),
        mag=np.asarray([_str_to_float_or_nan(v) for v in data["mag"]],
                       dtype=float),
        z=np.asarray([_str_to_float_or_nan(v) for v in data["z"]],
                     dtype=float),
        label=np.asarray([str(v) for v in data["label"]], dtype=object),
        required_lam=new_required_lam,
        no_gap=_bool_col("no_gap"),
        extend_blue=_bool_col("extend_blue"),
        extend_red=_bool_col("extend_red"),
        protect=_bool_col("protect"),
        # Centration override — coerce to a length-N object array of
        # strings. Anything not in the canonical five-value set turns
        # into "" via :func:`vmpt.catalog._normalise_centration`.
        centration=_centration_col_from_editor(data, n),
        max_configs=np.asarray(
            [_str_to_float_or_nan(v)
             for v in (list(data.get("max_configs", [""] * n)) + [""] * n)[:n]],
            dtype=float,
        ),
        source_path=cat.source_path,
        extras=new_extras,
    )
    state["catalogs"][idx]["catalog"] = new_cat
    _rebuild_merged_catalog()
    _rebuild_shutter_catalog_index()
    _render_catalog_list()
    refresh_overlays()
    _set_status(
        f"Applied edits → catalog #{idx + 1} now has {n} rows.",
        "ok", clear_after=10,
    )
    _cat_edit_close()


def _fmt_int_or_blank(v) -> str:
    """Render a float/int/string cell value as an integer string, or "".
    NaN / None / empty → "". Used by CSV save for priority + weight,
    which are stored as floats internally."""
    try:
        if isinstance(v, str):
            s = v.strip()
            if not s or s.lower() in ("nan", "none", "null", "--"):
                return ""
            f = float(s)
        else:
            f = float(v)
        if not np.isfinite(f):
            return ""
        return str(int(round(f)))
    except (ValueError, TypeError):
        return "" if v is None else str(v).strip()


def _fmt_float_or_blank(v) -> str:
    """Render a float/int/string cell value as a clean float string.
    NaN / None / empty → "". Trailing zeros / trailing decimal point
    are stripped so the CSV stays tidy for hand-editing. Used by CSV
    save for `z` (v1.3.2+) which is float-stored in the editor."""
    try:
        if isinstance(v, str):
            s = v.strip()
            if not s or s.lower() in ("nan", "none", "null", "--"):
                return ""
            f = float(s)
        else:
            f = float(v)
        if not np.isfinite(f):
            return ""
        s = f"{f:.6f}".rstrip("0").rstrip(".")
        return s or "0"
    except (ValueError, TypeError):
        return "" if v is None else str(v).strip()


def _on_cat_edit_save_csv():
    """Write the working copy to CSV at the user-supplied path.

    Wrapped in `_confirm_overwrite_if_exists` so the user gets a
    confirmation dialog when the target path already exists. v1.3.3+.
    """
    path = (cat_edit_csv_path_input.value or "").strip()
    if not path:
        _set_status("Enter a save path first.", "warn")
        return
    data = _cat_edit_source.data
    n = len(data["ra"])
    if n == 0:
        _set_status("Nothing to save — table is empty.", "warn")
        return
    p = Path(path).expanduser()

    def _do_save():
        _cat_edit_save_csv_to(p, data, n)

    _confirm_overwrite_if_exists(p, _do_save, what="CSV")


def _cat_edit_save_csv_to(p: Path, data: dict, n: int) -> None:
    """Inner CSV writer extracted from `_on_cat_edit_save_csv` so the
    overwrite-confirmation flow can call it either directly (when the
    path is new) or after the user confirms (when it exists)."""
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        import csv as _csv
        # Column order: the eight standard columns first, then the
        # per-target constraint columns (only when at least one row
        # has them set — otherwise the CSV stays unchanged from v1.2.x
        # for users who never touch the Constraints… popover), then
        # any user-added extras. Skip the internal _idx /
        # _has_constraint columns. Priority + weight are stored as
        # floats in source.data so we coerce them to int-or-blank.
        _INTERNAL = {"_idx", "_has_constraint", "lam_req",
                     "no_gap", "extend_blue", "extend_red", "protect",
                     "centration", "max_configs"}
        extras_cols = [k for k in data.keys()
                       if k not in _CAT_STD_COLS and k not in _INTERNAL]
        # Are any per-target constraints set anywhere in the catalog?
        # If not, omit the constraint columns so v1.2.x users get the
        # same CSV format they had before.
        has_constraints = any(
            int(v) for v in (data.get("_has_constraint") or [])
        )
        constraint_cols = (["lam_req", "no_gap", "extend_blue",
                            "extend_red", "protect", "centration",
                            "max_configs"]
                           if has_constraints else [])
        header = ["ID", "RA", "DEC", "priority", "weight", "mag", "z",
                  "label", *constraint_cols, *extras_cols]
        weights = data.get("weight", [""] * n)
        with open(p, "w", newline="") as f:
            w = _csv.writer(f)
            w.writerow(header)
            for i in range(n):
                row_vals = [
                    data["id"][i], data["ra"][i], data["dec"][i],
                    _fmt_int_or_blank(data["priority"][i]),
                    _fmt_int_or_blank(weights[i]),
                    # mag is still string-stored in the editor; z
                    # became float in v1.3.2 so we format it like a
                    # float (3+ decimals, NaN → blank).
                    data["mag"][i], _fmt_float_or_blank(data["z"][i]),
                    data["label"][i],
                ]
                for k in constraint_cols:
                    val = (data.get(k) or [""])[i]
                    if k == "lam_req":
                        row_vals.append("" if not val else str(val))
                    elif k in ("centration", "max_configs"):
                        # String label / count or empty — write verbatim
                        # (already normalised by the popover/loader).
                        row_vals.append(str(val) if val else "")
                    else:
                        row_vals.append("1" if int(val or 0) else "")
                for k in extras_cols:
                    row_vals.append(data[k][i])
                w.writerow(row_vals)
        _set_status(f"Saved {n} rows → {p}", "ok", clear_after=10)
    except Exception as e:  # noqa: BLE001
        _set_status(f"CSV save failed: {e}", "err")
        traceback.print_exc()


def _fmt_num(v: float) -> str:
    """Strip trailing zeros from a float so the CSV is tidy."""
    s = f"{float(v):.6f}".rstrip("0").rstrip(".")
    return s if s else "0"


catalog_edit_btn.on_click(_cat_edit_open)
cat_edit_close_btn.on_click(_cat_edit_close)
cat_edit_top_close_btn.on_click(_cat_edit_close)
cat_edit_select.on_change("value", _on_cat_edit_select)
cat_edit_columns_choice.on_change("value", _on_cat_edit_columns_change)
cat_edit_new_col_btn.on_click(_on_cat_edit_add_column)
cat_edit_compute_w_btn.on_click(_on_cat_edit_compute_w_from_p)
cat_edit_compute_p_btn.on_click(_on_cat_edit_compute_p_from_w)
cat_rule_apply_btn.on_click(_on_cat_rule_apply)
cat_edit_undo_btn.on_click(_on_cat_edit_undo)
cat_edit_redo_btn.on_click(_on_cat_edit_redo)
cat_edit_apply_btn.on_click(_on_cat_edit_apply)
cat_edit_csv_save_btn.on_click(_on_cat_edit_save_csv)

# Push every cell-edit's prior state onto the undo stack.
_cat_edit_source.on_change("data", _on_cat_edit_data_change)
# JS-side trash icons write to the signal source; Python deletes the row.
_cat_edit_delete_signal.on_change("data", _on_cat_edit_delete_signal)


def _cat_constraints_close(_event=None) -> None:
    cat_constraints_modal_backdrop.visible = False
    cat_constraints_modal_card.visible = False


def _cat_constraints_warn_for_lam(text: str) -> str:
    """Return an HTML warning string if any range in `text` falls
    outside the current (disperser, filter); empty string otherwise.

    Used as both the live-validation hint inside the popover and the
    save-time check. The warning is informational only — we still
    accept the save (user might be pre-staging for a future
    disperser).
    """
    from vmpt.catalog import _parse_lam_req_str
    from vmpt.wavelengths import disperser_range
    ranges = _parse_lam_req_str(text)
    if not ranges:
        return ""
    disp = (state.get("disperser") or "").upper()
    filt = (state.get("filter") or "").upper()
    rng = disperser_range(disp, filt) if disp and filt else None
    if rng is None:
        return ""
    lo_d, hi_d = rng
    out_of: list[str] = []
    for lo, hi in ranges:
        if lo < lo_d - 1e-3 or hi > hi_d + 1e-3:
            out_of.append(f"{lo:g}-{hi:g}")
    if not out_of:
        return ""
    return (
        f"<small style='color:#a05a30'>"
        f"⚠ {', '.join(out_of)} μm "
        f"is outside {disp} / {filt} ({lo_d:g}–{hi_d:g} μm). The "
        f"optimizer will never satisfy this constraint under the "
        f"current Disperser / Filter — change disperser or fix the "
        f"range. Save is still accepted in case you're pre-staging.</small>"
    )


def _on_cat_constraints_signal(attr, old, new) -> None:
    """JS-side click on Constraints… writes the row index here.

    Open the popover pre-filled with that row's current values.
    """
    global _cat_constraints_row_idx
    try:
        idx = int(new["idx"][0])
    except (KeyError, IndexError, TypeError, ValueError):
        return
    data = _cat_edit_source.data
    n = len(data.get("_idx", []))
    if not (0 <= idx < n):
        return
    _cat_constraints_row_idx = idx
    # Pre-fill the popover from source.data.
    id_str = str(data["id"][idx]) if "id" in data else f"row {idx}"
    cat_constraints_row_label.text = (
        f"<small style='color:#5a6b85'>Editing row "
        f"<b>#{idx + 1}</b> · ID <code>{id_str}</code></small>"
    )
    lam_val = str((data.get("lam_req") or [""] * n)[idx])
    cat_constraints_lam_input.value = lam_val
    cat_constraints_lam_warn.text = _cat_constraints_warn_for_lam(lam_val)
    active: list[int] = []
    if int((data.get("no_gap") or [0] * n)[idx]):
        active.append(0)
    if int((data.get("extend_blue") or [0] * n)[idx]):
        active.append(1)
    if int((data.get("extend_red") or [0] * n)[idx]):
        active.append(2)
    if int((data.get("protect") or [0] * n)[idx]):
        active.append(3)
    cat_constraints_checks.active = active
    # Centration override Select. Stored values are the canonical
    # five labels or ""; the Select's "(use global)" option maps to "".
    centration_stored = str(
        (data.get("centration") or [""] * n)[idx]
    ).strip()
    if centration_stored and centration_stored in (
        "UNCONSTRAINED", "ENTIRE_OPEN", "MIDPOINT",
        "CONSTRAINED", "TIGHTLY_CONSTRAINED",
    ):
        cat_constraints_centration_select.value = centration_stored
    else:
        cat_constraints_centration_select.value = "(use global)"
    # Max-configs cap Select. Stored "" / "1" / "2".
    mc_stored = str((data.get("max_configs") or [""] * n)[idx]).strip()
    cat_constraints_max_configs_select.value = (
        mc_stored if mc_stored in ("1", "2") else "(use global)"
    )
    cat_constraints_modal_backdrop.visible = True
    cat_constraints_modal_card.visible = True


def _on_cat_constraints_lam_change(attr, old, new) -> None:
    """Live-warn the user if any range falls outside the disperser."""
    cat_constraints_lam_warn.text = _cat_constraints_warn_for_lam(new or "")


def _on_cat_constraints_apply() -> None:
    """Write the popover values back into the editor's source for the
    current row, recompute its _has_constraint flag, and close."""
    global _cat_constraints_row_idx
    idx = _cat_constraints_row_idx
    if idx < 0:
        _cat_constraints_close()
        return
    data = dict(_cat_edit_source.data)
    n = len(data.get("_idx", []))
    if not (0 <= idx < n):
        _cat_constraints_close()
        return

    lam_val = (cat_constraints_lam_input.value or "").strip()
    active = set(cat_constraints_checks.active)
    new_no_gap = int(0 in active)
    new_blue = int(1 in active)
    new_red = int(2 in active)
    new_protect = int(3 in active)
    # Centration override — "(use global)" → empty string in storage.
    cent_val = (cat_constraints_centration_select.value or "").strip()
    if cent_val == "(use global)":
        cent_val = ""
    mc_val = (cat_constraints_max_configs_select.value or "").strip()
    if mc_val == "(use global)":
        mc_val = ""
    has_any = int(
        bool(lam_val) or new_no_gap or new_blue or new_red or new_protect
        or bool(cent_val) or bool(mc_val)
    )

    # Bokeh ColumnDataSource needs the WHOLE column replaced — mutating
    # an existing list in place doesn't fire the change notification.
    def _replace(field: str, value, default):
        col = list(data.get(field) or [default] * n)
        while len(col) < n:
            col.append(default)
        col[idx] = value
        data[field] = col

    _replace("lam_req", lam_val, "")
    _replace("no_gap", new_no_gap, 0)
    _replace("extend_blue", new_blue, 0)
    _replace("extend_red", new_red, 0)
    _replace("protect", new_protect, 0)
    _replace("centration", cent_val, "")
    _replace("max_configs", mc_val, "")
    _replace("_has_constraint", has_any, 0)
    # Use the silent setter so this edit doesn't push an extra undo
    # step on top of whatever the user was already doing.
    _cat_edit_set_data_silently(data)
    # But DO add it to the undo stack, manually — single push so undo
    # reverts the constraint change atomically.
    _cat_edit_undo_stack.append({k: list(v) for k, v in data.items()})
    if len(_cat_edit_undo_stack) > _CAT_EDIT_MAX_HISTORY:
        _cat_edit_undo_stack.pop(0)
    _cat_edit_redo_stack.clear()
    _cat_edit_render_history()
    _cat_constraints_row_idx = -1
    _cat_constraints_close()


cat_constraints_apply_btn.on_click(_on_cat_constraints_apply)
cat_constraints_cancel_btn.on_click(_cat_constraints_close)
cat_constraints_top_close_btn.on_click(_cat_constraints_close)
cat_constraints_lam_input.on_change("value", _on_cat_constraints_lam_change)
_cat_edit_constraint_signal.on_change("data", _on_cat_constraints_signal)


# Install the JS-side `window.__vmpt_delete_row(idx)` — called from
# the trash-icon column's `onclick`. The function closes over `sig`
# (our Python-side signal ColumnDataSource) so it can publish a
# change Python listens to. We bind to TWO triggers for reliability:
#
#   1. `DocumentReady` on the document — covers the case where the
#      table is populated by autoload (sys.argv → catalog → editor).
#   2. `source.js_on_change("data")` — fires every time the editor's
#      data is repopulated. The `if (!window.__vmpt_delete_row)` guard
#      makes re-installation a no-op.
_cat_edit_install_js = CustomJS(
    args=dict(sig=_cat_edit_delete_signal,
              csig=_cat_edit_constraint_signal),
    code="""
    if (!window.__vmpt_delete_row) {
        window.__vmpt_delete_row = function(idx) {
            // Coerce to int (template substitutes a numeric literal).
            const n = parseInt(idx, 10);
            // Stamp ensures every click is a fresh `data` value,
            // even if the same row is deleted twice in a session.
            sig.data = {idx: [n], stamp: [Date.now() + Math.random()]};
            sig.properties.data.change.emit();
        };
    }
    // Mirror installation for the per-row Constraints… button (the
    // template's onclick calls `window.__vmpt_open_constraints(idx)`).
    if (!window.__vmpt_open_constraints) {
        window.__vmpt_open_constraints = function(idx) {
            const n = parseInt(idx, 10);
            csig.data = {idx: [n], stamp: [Date.now() + Math.random()]};
            csig.properties.data.change.emit();
        };
    }
    """,
)
curdoc().js_on_event(DocumentReady, _cat_edit_install_js)
_cat_edit_source.js_on_change("data", _cat_edit_install_js)

# Hover tooltips on the main action buttons. A Bokeh Button renders its
# <button> inside the widget's shadow root, so the document-level <style>
# can't add a `title`; instead we walk the shadow DOM on DocumentReady
# (+ a few retries for async-rendered widgets) and set the native title
# attribute by matching the button label. (Section *titles* get their
# tooltips directly via the `title=` attr in `_section_header`.)
_btn_tips_install_js = CustomJS(code="""
const TIPS = [
  ["Open optimizer", "Search (RA, Dec, V3 PA) for the placement that captures the most targets in operable, well-centred shutters. Inside: count / weight / priority-tier ranking, the ΔRA/ΔDec/ΔPA box, collision protection, and per-config budgets."],
  ["View MPT catalog", "List the sources placed in open shutters across all configs, with each source's wavelength coverage and any detector gap."],
  ["Compute allowed V3 PA", "Query the JWST visibility tool for the V3 PA range allowed on the entered date."],
  ["Load image", "Open a background image: an example field, a FITS file, or a JPG/PNG with a WCS sidecar."],
  ["Load catalog", "Add a target catalog (CSV / ASCII / FITS with ID, RA, Dec). Add several to layer them."],
  ["Edit catalog", "Spreadsheet-style catalog editor: edit / sort / delete rows, set per-source constraints, save as CSV."],
  ["Import", "Bring in an APT/MPT plan (JSON), a shutter-mask CSV, an APT program (.aptx / ID), or a saved vMPT session."],
  ["Save session", "Write a vMPT session bundle you can reopen later or share with a collaborator."],
  ["Export to APT", "Write the eMPT bundle + the APT-importable plan + .cat for this design."],
  ["Reset display to defaults", "Restore layers, slitlet size, overlay alpha / stroke, and canvas size to their defaults."]
];
function setTitles(root){
  let bs; try { bs = root.querySelectorAll("button"); } catch (e) { return; }
  bs.forEach(function(b){
    const t = (b.textContent || "").replace(/\\s+/g, " ").trim();
    if (!t) return;
    for (const kv of TIPS){ if (t.indexOf(kv[0]) !== -1){ if (b.title !== kv[1]) b.title = kv[1]; break; } }
  });
  root.querySelectorAll("*").forEach(function(e){ if (e.shadowRoot) setTitles(e.shadowRoot); });
}
function applyTips(){ try { setTitles(document); } catch (e) {} }
applyTips();
[600, 1800, 4000].forEach(function(ms){ setTimeout(applyTips, ms); });
""")
curdoc().js_on_event(DocumentReady, _btn_tips_install_js)


# Install the global the top-bar CONFIG chip's onclick calls. Stamps the
# `_config_chip_signal` data so the Python `_on_config_chip_signal`
# listener fires and advances the active config. Bound to DocumentReady
# AND to every repaint of the stats bar (`stats_div.text`) so the global
# always exists by the time the chip is rendered; the guard makes
# re-installation a no-op.
_config_chip_install_js = CustomJS(
    args=dict(sig=_config_chip_signal),
    code="""
    if (!window.__vmpt_cycle_config) {
        window.__vmpt_cycle_config = function() {
            sig.data = {stamp: [Date.now() + Math.random()]};
            sig.properties.data.change.emit();
        };
    }
    """,
)
curdoc().js_on_event(DocumentReady, _config_chip_install_js)
stats_div.js_on_change("text", _config_chip_install_js)


def _on_config_chip_signal(attr, old, new):
    """Top-bar CONFIG chip was clicked → cycle to the next config
    (wraps 1→2→…→1). No-op in single-config mode. The dropdown in the
    Pointing tab and the chip's own colour stay in sync because
    `_switch_active_config` repaints both."""
    n = int(state.get("n_configs", 1))
    if n <= 1:
        return
    nxt = (int(state.get("active_config", 0)) + 1) % n
    _switch_active_config(nxt)


_config_chip_signal.on_change("data", _on_config_chip_signal)


# ── Draggable modal dialogs (v1.3.3+) ────────────────────────────────────
# Every pop-up card has TWO css classes that work together:
#   * .vmpt-modal-card    on the modal column itself
#   * .vmpt-modal-header  on a header row at the top of the card
#                         (contains the title + a clearly-visible ×
#                         close button)
#
# Only the **header** is a drag handle — `cursor: move` is on the
# header only, and the drag JS targets ONLY `.vmpt-modal-header`. The
# body of the modal is fully interactive: form controls, SlickGrid
# cells, sliders, buttons all receive their clicks normally with no
# skip-list logic to maintain.
#
# On the FIRST drag of a session the card's centred transform
# (`translate(-50%, -50%)`) is converted to absolute top/left so
# subsequent drags do straight arithmetic. The transform stays in
# place on modal open (so dialogs open centred on first display),
# but is cleared on the first header-mousedown — they stay where
# the user puts them after that.
_vmpt_drag_init_js = CustomJS(code=r"""
    if (window.__vmpt_drag_installed) return;
    window.__vmpt_drag_installed = true;

    let isDragging = false;
    let target = null;
    let startX = 0, startY = 0;
    let origLeft = 0, origTop = 0;

    // Walk the event's composedPath — that's the only way to "see"
    // through Bokeh's per-model shadow DOMs from a document-level
    // listener. Returns the first ancestor matching `selector` or
    // null when none of the path elements match.
    function pathFind(path, predicate) {
        for (const el of path) {
            if (el && el.classList && predicate(el)) return el;
        }
        return null;
    }

    function onMouseDown(e) {
        if (e.button !== 0) return;  // left-click only
        const path = e.composedPath ? e.composedPath() : [];
        // Don't drag when the user clicks on a Bokeh button (× close,
        // Done, Apply, etc.) or any focusable form control.
        const btn = pathFind(path, el => el.classList.contains("bk-btn"));
        if (btn) return;
        // The header is the drag handle. Find the header in the
        // event path; bail if the click started outside any header.
        const header = pathFind(path, el =>
            el.classList.contains("vmpt-modal-header"));
        if (!header) return;
        // The MODAL CARD is what we move. Find it in the same path —
        // it's an ancestor of the header in the shadow tree.
        const card = pathFind(path, el =>
            el.classList.contains("vmpt-modal-card"));
        if (!card) return;
        const rect = card.getBoundingClientRect();
        card.style.top = rect.top + "px";
        card.style.left = rect.left + "px";
        card.style.transform = "none";
        isDragging = true;
        target = card;
        startX = e.clientX;
        startY = e.clientY;
        origLeft = rect.left;
        origTop = rect.top;
        e.preventDefault();
    }

    function onMouseMove(e) {
        if (!isDragging || target === null) return;
        const newLeft = origLeft + (e.clientX - startX);
        const newTop = origTop + (e.clientY - startY);
        // Loose off-screen clamp — keep the top-left corner inside
        // the viewport so the user can always grab the dialog back.
        target.style.left = Math.max(
            -Math.max(0, target.offsetWidth - 80),
            Math.min(window.innerWidth - 80, newLeft)
        ) + "px";
        target.style.top = Math.max(
            0,
            Math.min(window.innerHeight - 32, newTop)
        ) + "px";
    }

    function onMouseUp() {
        isDragging = false;
        target = null;
    }

    // Document-level capture so we see mousedowns on elements inside
    // every Bokeh shadow root. Capture phase + composedPath is the
    // only way to "look into" shadow trees from outside.
    document.addEventListener("mousedown", onMouseDown, true);
    document.addEventListener("mousemove", onMouseMove);
    document.addEventListener("mouseup", onMouseUp);
""")
curdoc().js_on_event(DocumentReady, _vmpt_drag_init_js)

# `_MODAL_HEADER_STYLES` is defined near the top of the file
# (alongside SIDEBAR_W) so it's available before any modal card is
# built. The drag JS above only relies on the css_class for
# targeting; the visual styling lives there.
# Attach the GlobalInlineStyleSheet to the main figure's stylesheets
# list — Bokeh emits it into <head> for the whole document. (Adding
# it as a root via curdoc().add_root() isn't the documented API for
# stylesheets; `stylesheets=` on any LayoutDOM is.)
# Header bar visuals are now inline on each header row via
# `_MODAL_HEADER_STYLES` — see the modal definitions above. No
# document-level CSS to attach here.

# SlickGrid intercepts Cmd-C / Ctrl-C at the grid container level and
# copies the entire selected column / row from its internal data
# model, which defeats the in-cell drag-select + copy workflow. We
# install a CAPTURE-phase keydown listener at the document level that
# runs BEFORE SlickGrid: if an input/textarea is focused and a
# selection exists, we copy that selection ourselves via
# `navigator.clipboard.writeText`, then `preventDefault` + `stop-
# ImmediatePropagation` so SlickGrid never sees the event.
_cat_edit_clipboard_js = CustomJS(code="""
if (!window.__vmpt_copy_installed) {
    window.__vmpt_copy_installed = true;
    document.addEventListener('keydown', function(e) {
        if (!(e.metaKey || e.ctrlKey)) return;
        const k = (e.key || '').toLowerCase();
        if (k !== 'c' && k !== 'x') return;
        const a = document.activeElement;
        if (!a) return;
        if (a.tagName !== 'INPUT' && a.tagName !== 'TEXTAREA') return;
        const s = a.selectionStart, t = a.selectionEnd;
        if (s == null || t == null || s === t) return;
        const text = (a.value || '').substring(s, t);
        try {
            navigator.clipboard.writeText(text);
        } catch (err) {
            // Fallback for non-secure contexts.
            const ta = document.createElement('textarea');
            ta.value = text;
            document.body.appendChild(ta);
            ta.select();
            try { document.execCommand('copy'); } catch (_) {}
            document.body.removeChild(ta);
        }
        if (k === 'x') {
            // Implement cut: clear the selected slice from the input.
            a.value = (a.value || '').substring(0, s) +
                      (a.value || '').substring(t);
            a.selectionStart = a.selectionEnd = s;
            a.dispatchEvent(new Event('input', { bubbles: true }));
        }
        e.preventDefault();
        e.stopImmediatePropagation();
    }, true);
}
""")
curdoc().js_on_event(DocumentReady, _cat_edit_clipboard_js)
# Re-install on every catalog populate too — covers the case where
# the page never fired DocumentReady (rare) or the function got wiped.
_cat_edit_source.js_on_change("data", _cat_edit_clipboard_js)

# ── Keyboard panning of the canvas ───────────────────────────────────
# WASD and the arrow keys pan the figure view (shift the x/y ranges),
# mirroring a left-drag with the pan tool: a fixed fraction of the
# visible span per press (proportional to zoom), and holding a key
# auto-repeats for a glide. Shift = coarser step. The shift is by the
# *signed* span so it stays correct on the flipped RA axis (← always
# moves the view left on screen, etc.). Skipped while a text field is
# focused or any modal is open, so it never steals the arrow keys from
# inputs or the catalog editor's cell navigation. Pure client-side, so
# panning is instant (no server round-trip per keypress).
# Keyboard-pan → re-render bridge. A JS `setv` on the ranges (below) pans
# the view but — unlike the drag pan tool — does NOT emit a RangesUpdate
# event, so `on_ranges_update` never fires and the shutter overlays don't
# re-cull to the newly exposed region. The keypan handler bumps this
# source once the pan settles (debounced, mirroring drag-pan's
# render-on-release); Python then re-runs refresh_overlays.
keypan_trigger = ColumnDataSource(data=dict(n=[0]))
keypan_trigger.on_change("data", lambda attr, old, new: refresh_overlays())
# Park the trigger source on an invisible, zero-size renderer so the CDS
# is part of the figure's (hence the document's) model graph. Without a
# renderer the CDS is only referenced by the CustomJS args and never gets
# serialized into the session, so a JS `.data` write would never sync back
# to Python. The figure's DataRange1d is pinned to `img_glyph` above, so
# this phantom point can never affect auto-ranging or appear on screen.
fig.scatter(x="n", y="n", source=keypan_trigger,
            size=0, fill_alpha=0, line_alpha=0)

# Hover + spacebar → toggle the single hovered shutter. The keydown handler
# (below) reads the last hover position (tracked in JS via MouseMove) and
# writes it here; Python then snaps to the nearest shutter and toggles it.
# Same phantom-renderer trick as keypan_trigger so the JS `.data` write
# syncs back to the server.
spacetoggle_trigger = ColumnDataSource(data=dict(x=[0.0], y=[0.0], n=[0]))
fig.scatter(x="x", y="y", source=spacetoggle_trigger,
            size=0, fill_alpha=0, line_alpha=0)


def _on_spacebar_toggle(attr, old, new) -> None:
    """Toggle the single shutter under the cursor (from the hover position
    the keydown handler stashed). Mirrors on_tap's snap-to-nearest, but
    always acts on ONE shutter via `_toggle_single_shutter_at`."""
    img = state["image"]
    fiducial = _pointing_skycoord()
    if img is None or fiducial is None:
        return
    try:
        xs = spacetoggle_trigger.data.get("x") or []
        ys = spacetoggle_trigger.data.get("y") or []
        x_data, y_data = float(xs[0]), float(ys[0])
    except (TypeError, ValueError, IndexError):
        return
    try:
        sky = img.wcs.pixel_to_world(x_data, y_data)
        v2, v3 = _sky_to_v2v3(sky, fiducial, state["pa_v3"])
        nearest = _nearest_shutter(v2, v3, require_operable=False,
                                   max_dist_arcsec=0.5)
    except Exception:  # noqa: BLE001
        return
    if nearest is None:
        return  # cursor wasn't over a shutter
    q, s, d = nearest
    _toggle_single_shutter_at(q, s, d)


spacetoggle_trigger.on_change("data", _on_spacebar_toggle)

_canvas_keypan_js = CustomJS(
    args=dict(xr=fig.x_range, yr=fig.y_range, trig=keypan_trigger,
              trig2=spacetoggle_trigger), code="""
if (window.__vmptKeypanInstalled) return;
window.__vmptKeypanInstalled = true;
const KEYS = {
  ArrowLeft:'L', a:'L', A:'L', ArrowRight:'R', d:'R', D:'R',
  ArrowUp:'U',   w:'U', W:'U', ArrowDown:'D',  s:'D', S:'D',
};
function modalOpen() {
  const out = [];
  (function walk(r){
    r.querySelectorAll('.vmpt-modal-card').forEach(c => out.push(c));
    r.querySelectorAll('*').forEach(el => { if (el.shadowRoot) walk(el.shadowRoot); });
  })(document);
  return out.some(c => c.offsetParent !== null);
}
document.addEventListener('keydown', function(e) {
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  const isSpace = (e.key === ' ' || e.code === 'Space');
  const dir = KEYS[e.key];
  if (!isSpace && !dir) return;
  // Bokeh widgets render their <input>/<select> inside shadow roots, so
  // document.activeElement is the shadow HOST — pierce into the focused
  // shadow tree to find the real focused element, else arrows/WASD/space
  // would fire while the user is typing in a field (filters, V3 PA, …).
  let ae = document.activeElement;
  while (ae && ae.shadowRoot && ae.shadowRoot.activeElement) ae = ae.shadowRoot.activeElement;
  if (ae) {
    const tag = (ae.tagName || '').toUpperCase();
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || ae.isContentEditable) return;
  }
  if (modalOpen()) return;
  // Spacebar → toggle the single shutter currently under the cursor. Only
  // fires while the mouse is over the canvas (hover position known); off
  // the canvas we leave the key alone (don't swallow it / scroll the page).
  if (isSpace) {
    const h = window.__vmptHoverXY;
    if (!h) return;
    e.preventDefault();
    trig2.data = {x: [h.x], y: [h.y],
                  n: [(((trig2.data.n || [0])[0]) | 0) + 1]};
    return;
  }
  if (xr.start == null || xr.end == null || yr.start == null || yr.end == null) return;
  e.preventDefault();
  const f = e.shiftKey ? 0.25 : 0.10;
  const dx = xr.end - xr.start, dy = yr.end - yr.start;
  if      (dir === 'R') xr.setv({start: xr.start + f*dx, end: xr.end + f*dx});
  else if (dir === 'L') xr.setv({start: xr.start - f*dx, end: xr.end - f*dx});
  else if (dir === 'U') yr.setv({start: yr.start + f*dy, end: yr.end + f*dy});
  else if (dir === 'D') yr.setv({start: yr.start - f*dy, end: yr.end - f*dy});
  // Re-cull the shutter overlays to the new view once the pan settles.
  // Debounced so holding a key (key-repeat) glides smoothly and only
  // refreshes after the last step, like the drag pan tool on release.
  if (window.__vmptKeypanTimer) clearTimeout(window.__vmptKeypanTimer);
  window.__vmptKeypanTimer = setTimeout(function() {
    trig.data = {n: [(((trig.data.n || [0])[0]) | 0) + 1]};
  }, 200);
}, true);
""")
curdoc().js_on_event(DocumentReady, _canvas_keypan_js)

# Track the cursor's data-space position over the canvas so the spacebar
# shortcut (above) knows which shutter is hovered. MouseMove fires only over
# the plot frame and `cb_obj.x/y` are already in image-pixel data coords —
# pure JS, no server round-trip. MouseLeave clears it so spacebar off the
# canvas is a no-op.
fig.js_on_event(MouseMove, CustomJS(code="window.__vmptHoverXY = {x: cb_obj.x, y: cb_obj.y};"))
fig.js_on_event(MouseLeave, CustomJS(code="window.__vmptHoverXY = null;"))

# When the user clicks a SlickGrid column header to sort, the table
# stays at whatever row was on screen — so a sort hidden 200 rows
# down looks like nothing happened. Install a document-level click
# delegate that scrolls the affected table's viewport back to row 0
# after every header click.
_cat_edit_sort_scroll_js = CustomJS(code="""
if (!window.__vmpt_sort_scroll_installed) {
    window.__vmpt_sort_scroll_installed = true;
    document.addEventListener('click', function(e) {
        const header = e.target.closest('.slick-header-column');
        if (!header) return;
        // Walk up to the SlickGrid container, then find its viewport.
        // Bokeh wraps the DataTable in a `.bk-DataTable` (Bokeh 3.x)
        // — but the SlickGrid viewport class is stable.
        const grid = header.closest('.slick-pane') ||
                     header.closest('.slick-container') ||
                     header.closest('.bk-DataTable') ||
                     header.parentElement;
        if (!grid) return;
        const viewport = grid.querySelector('.slick-viewport') ||
                         grid.parentElement.querySelector('.slick-viewport');
        if (!viewport) return;
        // SlickGrid commits the sort synchronously on click; the
        // re-render runs immediately after. A small delay ensures
        // we scroll the NEW layout, not the pre-sort one.
        setTimeout(function() { viewport.scrollTop = 0; }, 80);
    });
}
""")
curdoc().js_on_event(DocumentReady, _cat_edit_sort_scroll_js)
_cat_edit_source.js_on_change("data", _cat_edit_sort_scroll_js)

# Keep the user's scroll position when they edit a cell. When a column
# sort is active, Bokeh re-sorts the table on every data change and the
# view follows the moved row — so committing an edit otherwise flings the
# viewport to wherever the edited row lands, which is jarring when editing
# a long catalog. The fix is two parts:
#   • SETUP: attach a capture-phase `mousedown` listener to the editor's
#     SlickGrid viewport that snapshots the scroll position the instant
#     the user presses on a cell — i.e. BEFORE the edit commits and the
#     re-sort jumps the view. (Capturing later, on the commit, is too
#     late: the jump has already happened.) Attached with a short retry
#     since the viewport only renders once the editor modal is open.
#   • RESTORE: on every `patching` event (the client cell-commit and the
#     server-side numeric coercion both fire one), pin the viewport back
#     to the snapshot after the re-render settles.
# Keyed to `patching`/`mousedown` only, so it never touches header-click
# sorts (those reorder the view without patching source.data) — the
# deliberate scroll-to-top-after-sort behaviour above is preserved.
_VMPT_FIND_CAT_VP = """
function _vmptCatVP() {
    const out = [];
    (function walk(root) {
        root.querySelectorAll('.slick-viewport').forEach(v => out.push(v));
        root.querySelectorAll('*').forEach(el => { if (el.shadowRoot) walk(el.shadowRoot); });
    })(document);
    return out.find(v => v.scrollHeight > 1000 && v.offsetParent !== null)
        || out.find(v => v.scrollHeight > 1000) || null;
}
"""
_cat_edit_scroll_setup_js = CustomJS(code=_VMPT_FIND_CAT_VP + """
let tries = 0;
(function attach() {
    const vp = _vmptCatVP();
    if (!vp) { if (tries++ < 25) setTimeout(attach, 150); return; }
    window.__vmptCatVP = vp;
    if (!vp.__vmptKeep) {
        vp.__vmptKeep = true;
        vp.__vmptSaved = null;
        // Snapshot the scroll the instant the user presses on a cell —
        // before the edit commits and the table re-sorts.
        vp.addEventListener('mousedown', () => { vp.__vmptSaved = vp.scrollTop; }, true);
    }
})();
""")
_cat_edit_scroll_restore_js = CustomJS(code=_VMPT_FIND_CAT_VP + """
let vp = window.__vmptCatVP;
if (!vp || vp.offsetParent === null || vp.scrollHeight < 1000) {
    vp = _vmptCatVP();
    window.__vmptCatVP = vp;
}
if (!vp || vp.__vmptSaved == null) return;
const saved = vp.__vmptSaved;
const restore = () => {
    const max = vp.scrollHeight - vp.clientHeight;
    vp.scrollTop = Math.max(0, Math.min(saved, max));
};
// Re-pin across the synchronous re-sort and any async follow-up render.
setTimeout(restore, 0);
setTimeout(restore, 90);
setTimeout(restore, 240);
""")
curdoc().js_on_event(DocumentReady, _cat_edit_scroll_setup_js)
_cat_edit_source.js_on_change("data", _cat_edit_scroll_setup_js)
catalog_edit_btn.js_on_click(_cat_edit_scroll_setup_js)
_cat_edit_source.js_on_change("patching", _cat_edit_scroll_restore_js)


def _apply_optimizer_result(
    ra_p: float, dec_p: float, pa_v3: float,
    *, clear_existing: bool = True, config_idx: int | None = None,
) -> None:
    """Apply one of the optimizer's pointings.

    1. Push current open_shutters to history so the Undo button
       reverts the whole apply in one step.
    2. If `clear_existing` is True, drop every previously open shutter
       before placing the optimizer's slitlets. (This is the default —
       the user is warned via a confirm dialog before the Apply
       button fires.)
    3. Set RA / Dec / V3 PA via the widgets.
    4. Use the in-flight `_opt_run["evaluator"]` to find, for this
       pointing, which sources are observable and where their shutter
       centres are; open an N-shutter slitlet (N from Setting tab)
       at each, auto-tagged with the catalog source ID.

    The optimizer's `axy_to_shutter` returns 0-based fractional
    indices (`s_frac ∈ [0,170]`, `d_frac ∈ [0,364]`); vMPT's
    `_add_slitlet` takes 1-based `s, d`. We convert here so the
    opened slitlets centre exactly on the target.

    `config_idx` (v1.4.0): when given, the result is applied to that
    MPT config — the config list is grown and the active config switched
    first, so the slitlets land on the right config. The per-config
    observation budget (config 1 → none, config 2 → the pass-2 mask) is
    applied to the evaluator so Config 2 never reopens Config 1's targets.
    """
    if config_idx is not None:
        if config_idx + 1 > int(state.get("n_configs", 1)):
            _ensure_n_configs(config_idx + 1)
            _prefs_save_suppress["flag"] = True
            try:
                mpt_num_configs_spinner.value = int(state["n_configs"])
            finally:
                _prefs_save_suppress["flag"] = False
        if config_idx != state["active_config"]:
            _switch_active_config(config_idx)

    _push_history()  # snapshot for Undo

    if clear_existing:
        _set_open_shutters({})

    ra_input.value = f"{ra_p:.6f}"
    dec_input.value = f"{dec_p:.6f}"
    _sync_pa_widgets(float(pa_v3))

    n_targets = 0
    n_opened = 0
    ev = _opt_run.get("evaluator") if _opt_run else None
    ids = _opt_run.get("source_ids", []) if _opt_run else []
    if ev is not None:
        # Apply the per-config evaluator budget (v1.4.0). config 0 / single
        # config → no budget; config 1 → the pass-2 budget mask so it
        # excludes Config 1's already-observed sources.
        try:
            _eff_cfg = (config_idx if config_idx is not None
                        else int(state.get("active_config", 0)))
            _budget = (_opt_run.get("budgets", {}) or {}).get(_eff_cfg)
            if _budget is not None:
                ev._budget = np.asarray(_budget, dtype=bool)
                ev._budget_enabled = not bool(ev._budget.all())
            else:
                ev._budget = np.ones(len(ev.ra), dtype=bool)
                ev._budget_enabled = False
        except Exception:  # noqa: BLE001
            pass
        try:
            detected, _tp, shutters = ev.evaluate(ra_p, dec_p, pa_v3)
            quad, s_frac, d_frac = shutters
            for i in range(len(detected)):
                if not bool(detected[i]):
                    continue
                n_targets += 1
                q = int(quad[i])
                # Optimizer indices are 0-based; vMPT shutters use
                # 1-based `s, d`. Without the +1 the slitlet opens
                # one row up and one column left of the target.
                s = int(round(float(s_frac[i]))) + 1
                d = int(round(float(d_frac[i]))) + 1
                if not (1 <= q <= 4 and 1 <= s <= 171 and 1 <= d <= 365):
                    continue
                target_id = str(ids[i]) if i < len(ids) else None
                try:
                    if _add_slitlet(q, s, d, target_id=target_id) > 0:
                        n_opened += 1
                except Exception:  # noqa: BLE001
                    # Edge cases (e.g. slitlet partially off the MSA)
                    # are silently skipped; we still try the rest.
                    pass
        except Exception as e:  # noqa: BLE001
            _set_status(f"Auto-open after Apply failed: {e}", "warn")
            traceback.print_exc()

    n_shutters = int(state.get("slitlet_height", 3))
    cleared_note = " (cleared previous picks)" if clear_existing else ""
    _set_status(
        f"Applied: pointing → RA={ra_p:.5f}, Dec={dec_p:.5f}, "
        f"V3 PA={pa_v3:.2f}° · opened {n_opened}/{n_targets} "
        f"{n_shutters}-shutter slitlets{cleared_note}.",
        "ok", clear_after=14,
    )


def _on_opt_apply_trigger(attr, old, new):
    """Fires when the Apply button's JS confirm dialog returns OK.
    The trigger value is `"<ra>,<dec>,<pa>,<config_idx>,<stamp>"`."""
    if not new:
        return
    cfg = None
    try:
        parts = new.split(",", 4)
        ra_p = float(parts[0]); dec_p = float(parts[1]); pa_v3 = float(parts[2])
        if len(parts) >= 5:
            cfg = int(parts[3])
    except (ValueError, TypeError, IndexError):
        opt_apply_trigger.value = ""
        return
    # Reset the trigger so the next click on the same row still fires.
    # (Suppress the recursive on_change by checking `new` at the top.)
    opt_apply_trigger.value = ""
    _apply_optimizer_result(ra_p, dec_p, pa_v3, clear_existing=True,
                            config_idx=cfg)
    _opt_hide_modal()


opt_apply_trigger.on_change("value", _on_opt_apply_trigger)


def _on_opt_apply_both_trigger(attr, old, new):
    """Apply an N-config plan in one click: Config 1 #k + … + Config N #k.

    The trigger payload is ``"<rank>,<stamp>"`` (rank is 0-based); each
    config is filled from its own k-th solution (clamped to range)."""
    if not new:
        return
    opt_apply_both_trigger.value = ""
    try:
        k = int(str(new).split(",", 1)[0])
    except (ValueError, IndexError):
        k = 0
    # Leading run of configs that produced solutions (matches the table).
    results: list = []
    for r in (_opt_run.get("pass_results") or []):
        if len(r.get("score", [])) == 0:
            break
        results.append(r)
    if not results:
        return
    for ci, r in enumerate(results):
        kk = max(0, min(k, len(r["score"]) - 1))
        _apply_optimizer_result(
            float(r["ra"][kk]), float(r["dec"][kk]), float(r["pa"][kk]),
            clear_existing=True, config_idx=ci,
        )
    # Leave the user on Config 1 to review.
    _switch_active_config(0)
    _opt_hide_modal()
    _set_status(
        f"Applied {len(results)}-config plan #{k + 1} "
        f"(Config 1 … Config {len(results)}).", "ok", clear_after=12)


opt_apply_both_trigger.on_change("value", _on_opt_apply_both_trigger)


# ── Pop-up modal helpers ─────────────────────────────────────────────────


def _opt_show_modal() -> None:
    """Reveal the optimizer modal in its 'progress' state."""
    opt_modal_backdrop.visible = True
    opt_modal_card.visible = True
    opt_modal_progress_box.visible = True
    opt_modal_results_box.visible = False
    opt_modal_results_summary.text = ""
    opt_modal_results_rows.children = []


def _opt_hide_modal() -> None:
    opt_modal_backdrop.visible = False
    opt_modal_card.visible = False


def _opt_update_progress(text: str, frac: float) -> None:
    """Update the progress text + CSS progress bar. `frac` ∈ [0, 1].

    The bar's `text` is set ONCE at construction (so the stripe + pulse
    CSS animations keep running uninterrupted). Here we only swap the
    `--vmpt-pct` custom property on the wrapper Div's inline style;
    Bokeh applies that without replacing innerHTML, so animations
    don't restart.
    """
    pct = max(0.0, min(1.0, float(frac))) * 100.0
    opt_modal_progress_text.text = f"<i>{text}</i>"
    opt_modal_progress_bar.styles = {"--vmpt-pct": f"{pct:.1f}%"}


# Human-readable labels for the per-pointing drop reasons surfaced in the
# results Score-cell tooltips. Keys mirror `optimizer.DROP_REASONS`.
# `budget` (v1.4.0) only fires in the 2-config pass-2 search: a source is
# dropped because the other config already observed it up to its
# max-configs cap. Ordered budget-first since that dominates Config 2.
_DROP_REASON_LABELS = {
    "budget":       "already observed in another config",
    "collision":    "spectral collision",
    "required_lam": "required λ-range missing",
    "no_gap":       "detector gap inside spectrum",
    "extend_blue":  "blue-edge truncated",
    "extend_red":   "red-edge truncated",
}


def _opt_drop_breakdown_lines(reasons_i: dict, n_drop_i: int) -> list:
    """Tooltip lines breaking `n_drop_i` down by reason (best-effort).

    Falls back to a generic line when the per-reason dict is empty.
    """
    lines: list = []
    if reasons_i:
        lines.append(f"−{n_drop_i} dropped:")
        for key, label in _DROP_REASON_LABELS.items():
            count = int(reasons_i.get(key, 0) or 0)
            if count > 0:
                lines.append(f"   {count}× {label}")
    else:
        lines.append(f"(−{n_drop_i} source(s) dropped at this pointing.)")
    return lines


def _opt_build_result_rows(
    results: dict, ra_ref: float, dec_ref: float, pa_ref: float,
    *, method: str, config_idx: int, confirm_what: str,
) -> list:
    """Build the Bokeh rows (header + one per solution) for one results
    table. Shared by the single-config and 2-config (pair) renderers.

    `config_idx` is baked into each Apply button's trigger payload so the
    Python handler knows which config to fill. `confirm_what` is the body
    of the JS confirm() dialog.
    """
    n = len(results["score"])
    if n == 0:
        return []
    n_show = min(10, n)
    cos_dec = np.cos(np.deg2rad(dec_ref))
    breakdowns = results.get("tier_breakdown")
    totals = results.get("total_count")
    sum_weights = results.get("sum_weight")
    top_targets = results.get("top_targets")
    n_dropped = results.get("n_dropped")
    drop_reasons = results.get("drop_reasons")  # per-pointing dict (v1.3.0+)
    protect_enabled = bool(results.get("protect_enabled"))
    is_hierarchy = breakdowns is not None and len(breakdowns) >= n_show
    is_meritocracy = method == "Meritocracy"

    # Header row (Bokeh row of Divs aligning with the data rows).
    HEADER_BG = "#eef2f8"
    ROW_H = 26
    # `overflow: hidden + white-space: nowrap` keeps a too-wide label
    # from wrapping below the row (e.g. "Σw 2009.0 (46)" was wrapping
    # into a second line and rendering outside the modal). Long
    # strings get a hover ellipsis-truncate instead.
    CELL_STYLES = {
        "padding": "0 6px",
        "line-height": f"{ROW_H}px",
        "border-bottom": "1px solid #d8dee8",
        "box-sizing": "border-box",
        "height": f"{ROW_H}px",
        "overflow": "hidden",
        "white-space": "nowrap",
        "text-overflow": "ellipsis",
    }
    def _cell(text, width, *, header=False, bold=False, mono=False):
        styles = dict(CELL_STYLES)
        if header:
            styles["background"] = HEADER_BG
            styles["font-weight"] = "600"
            styles["text-align"] = "center"
        if bold:
            styles["font-weight"] = "600"
        if mono:
            styles["font-family"] = "monospace"
            styles["text-align"] = "right"
        else:
            styles["text-align"] = "center"
        return Div(text=text, width=width, height=ROW_H, styles=styles)

    # Score column width is method-dependent — Democracy shows a
    # single integer ("46"), Meritocracy a weight sum + count
    # ("Σw 287.0 (46)"), Hierarchy a per-tier breakdown +
    # count ("P0:4 · P1:12 · P2:30 (46)"). Narrow widths caused
    # the label to wrap below the row. When collision protection is
    # on, the label gains a " −K" suffix so the column needs to be
    # wider by ~30 px to keep the ellipsis from kicking in.
    if is_hierarchy:
        score_w = 200
    elif is_meritocracy:
        score_w = 140
    else:
        score_w = 80
    if protect_enabled:
        score_w += 36
    COL_WIDTHS = (32, score_w, 90, 90, 100, 168)
    total_w = sum(COL_WIDTHS)
    header_row = row(
        _cell("#",     COL_WIDTHS[0], header=True),
        _cell("Score", COL_WIDTHS[1], header=True),
        _cell("ΔRA",   COL_WIDTHS[2], header=True),
        _cell("ΔDec",  COL_WIDTHS[3], header=True),
        _cell("ΔPA",   COL_WIDTHS[4], header=True),
        _cell("",      COL_WIDTHS[5], header=True),  # spacer over the Apply col
        spacing=0, width=total_w,
    )
    data_rows: list = []
    for i in range(n_show):
        s = float(results["score"][i])
        ra_i = float(results["ra"][i])
        dec_i = float(results["dec"][i])
        pa_i = float(results["pa"][i])
        d_ra = (ra_i - ra_ref) * 3600.0 * cos_dec
        d_dec = (dec_i - dec_ref) * 3600.0
        d_pa = (pa_i - pa_ref + 180.0) % 360.0 - 180.0

        # Per-row headline score. Method-dependent:
        #   Democracy:   "<count>"
        #   Meritocracy: "<Σw> (<count>)"
        #   Hierarchy:   "P0:4 · P1:12 · …  (<count>)"
        total_i = (int(totals[i]) if totals is not None and i < len(totals)
                   else None)
        if is_hierarchy:
            # Only show tiers that actually placed a source — dropping the
            # "P5:0" empties keeps the label short and readable.
            breakdown_label = " · ".join(
                f"P{int(t)}:{c}" for t, c in breakdowns[i] if c > 0
            )
            if total_i is not None:
                score_label = f"{breakdown_label}  ({total_i})"
            else:
                score_label = breakdown_label or f"{s:.1f}"
        elif is_meritocracy:
            sw = (float(sum_weights[i])
                  if sum_weights is not None and i < len(sum_weights)
                  else s)
            score_label = (f"Σw {sw:.1f}  ({total_i})"
                           if total_i is not None else f"Σw {sw:.1f}")
        else:
            score_label = (str(total_i) if total_i is not None
                           else f"{s:.1f}")
        # Append "−K" when the collision-protection rules dropped any
        # detected sources at this pointing. Hover tooltip explains.
        n_drop_i = (int(n_dropped[i]) if n_dropped is not None
                    and i < len(n_dropped) else 0)
        if n_drop_i > 0:
            score_label = f"{score_label} −{n_drop_i}"

        # Hover tooltip — top-10 placed sources in priority order.
        tooltip = ""
        if top_targets is not None and i < len(top_targets):
            tt = top_targets[i]
            if tt:
                tip_lines = [f"Top {len(tt)} placed sources at this pointing:"]
                for r, t in enumerate(tt, 1):
                    pid = "—" if t["p"] is None else f"{int(t['p'])}"
                    wid = "—" if t["w"] is None else f"{int(t['w'])}"
                    # 🛡 prefix marks sources in the protected set
                    # (collision-protection mode). Plain row otherwise.
                    prot_mark = "🛡 " if t.get("prot") else "   "
                    tip_lines.append(
                        f"{r}. {prot_mark}ID={t['id']}  P={pid}  W={wid}"
                    )
                if n_drop_i > 0:
                    tip_lines.append("")
                    # Break the drop count down by reason if the
                    # per-pointing dict was populated (v1.3.0+).
                    reasons_i = (drop_reasons[i]
                                 if drop_reasons is not None
                                 and i < len(drop_reasons) else None)
                    if reasons_i:
                        tip_lines.append(f"−{n_drop_i} dropped:")
                        REASON_LABELS = {
                            "collision":    "spectral collision",
                            "required_lam": "required λ-range missing",
                            "no_gap":       "detector gap inside spectrum",
                            "extend_blue":  "blue-edge truncated",
                            "extend_red":   "red-edge truncated",
                        }
                        for key, label in REASON_LABELS.items():
                            count = int(reasons_i.get(key, 0) or 0)
                            if count > 0:
                                tip_lines.append(f"   {count}× {label}")
                    else:
                        tip_lines.append(
                            f"(−{n_drop_i} non-protected source(s) dropped "
                            f"to protect the 🛡 spectra above.)"
                        )
                tooltip = "\n".join(tip_lines)
        # `title=` renders as a native browser tooltip on hover.
        # HTML-escape minimally so quotes don't break the attribute.
        title_attr = (tooltip
                      .replace("&", "&amp;")
                      .replace('"', "&quot;")
                      .replace("<", "&lt;"))
        score_html = (
            f'<span title="{title_attr}" '
            f'style="cursor:help; border-bottom:1px dotted #888">'
            f'{score_label}</span>'
        )

        btn = Button(
            label=f"Apply #{i+1}",
            button_type="primary" if i == 0 else "default",
            width=COL_WIDTHS[5] - 12, height=ROW_H,
        )
        # Apply via JS confirm → trigger TextInput → Python handler.
        # The browser's `window.confirm` blocks the JS event loop so
        # the user has to OK before we tell Python to apply (which
        # CLEARS all previously open shutters). Cancel = no-op.
        #
        # NOTE: CustomJS.args only accepts Bokeh Model instances, not
        # plain floats. We embed the per-button scalars via Python
        # f-string interpolation into the JS body, and pass only the
        # trigger TextInput as a real model arg.
        _confirm_js = confirm_what.replace('"', '\\"')
        btn.js_on_click(CustomJS(
            args=dict(trig=opt_apply_trigger),
            code=f"""
            const ra = {float(ra_i)};
            const dec = {float(dec_i)};
            const pa = {float(pa_i)};
            const rank = {i + 1};
            const cfg = {int(config_idx)};
            const msg = "Apply solution #" + rank +
                        " ?\\n\\n{_confirm_js}";
            if (!window.confirm(msg)) {{
                return;
            }}
            // Payload: ra,dec,pa,config_idx,stamp. The stamp makes Bokeh
            // see a fresh value even when the same row is re-applied.
            trig.value = ra.toFixed(8) + "," + dec.toFixed(8) + "," +
                         pa.toFixed(6) + "," + cfg + "," +
                         Date.now() + "_" + Math.random();
            """,
        ))

        zebra = "#f7f9fc" if i % 2 else "#ffffff"
        row_styles_bg = lambda c=zebra: {"background": c}
        data_rows.append(row(
            _cell(str(i + 1), COL_WIDTHS[0]),
            _cell(score_html, COL_WIDTHS[1], bold=True),
            _cell(f"{d_ra:+.2f}″", COL_WIDTHS[2], mono=True),
            _cell(f"{d_dec:+.2f}″", COL_WIDTHS[3], mono=True),
            _cell(f"{d_pa:+.3f}°", COL_WIDTHS[4], mono=True),
            btn,
            spacing=0, width=total_w, styles=row_styles_bg(),
        ))

    return [header_row, *data_rows]


def _opt_render_results_in_modal(
    results: dict, ra_ref: float, dec_ref: float, pa_ref: float,
    n_sources: int,
    *, method: str = "Democracy",
) -> None:
    """Single-config results table (n_pass == 1). Applies to the active
    config (config 0 in the common case)."""
    n = len(results["score"])
    if n == 0:
        opt_modal_progress_box.visible = False
        opt_modal_results_box.visible = True
        opt_modal_results_summary.text = (
            "<i>No solutions found inside the search box. "
            "Try widening ΔRA / ΔDec / ΔPA, lowering the centration "
            "class, or relaxing the priority cutoff.</i>"
        )
        opt_modal_results_rows.children = []
        return
    n_show = min(10, n)
    method_blurb = {
        "Democracy":   "ranked by raw count.",
        "Meritocracy": "ranked by Σ weight; <b>Σw</b> shown, total count in parens.",
        "Hierarchy":   ("ranked lexicographically by priority tier; "
                        "<b>P<sub>i</sub>:n</b> = sources placed at tier i."),
    }.get(method, "")
    protect_blurb = ""
    if bool(results.get("protect_enabled")):
        protect_blurb = (
            " · <b>🛡 collision protection</b> ON — "
            "<code>−K</code> counts targets dropped to protect "
            "high-priority spectra."
        )
    opt_modal_results_summary.text = (
        f"<div style='font-size:12px; margin-bottom:4px'>"
        f"<b>Top {n_show}</b> distinct solutions of {n} total · "
        f"{n_sources} candidate sources · "
        f"<b>Method:</b> {method} — {method_blurb}{protect_blurb} "
        f"Hover a Score cell to see the top-10 placed sources. "
        f"Offsets are from the search centre.</div>"
    )
    rows_built = _opt_build_result_rows(
        results, ra_ref, dec_ref, pa_ref, method=method,
        config_idx=int(state.get("active_config", 0)),
        confirm_what=("This will CLEAR all previously open shutters and "
                      "replace them with the optimizer's slitlets."),
    )
    opt_modal_progress_box.visible = False
    opt_modal_results_box.visible = True
    opt_modal_results_rows.children = rows_built


def _opt_render_multi_results(
    results_list: list, ra_ref: float, dec_ref: float, pa_ref: float,
    n_sources: int, *, method: str = "Democracy",
) -> None:
    """N-config (auto-all) results as ONE combined table.

    Each row is a *plan* = Config 1 #k + Config 2 #k + … + Config N #k,
    rendered as one line per config sharing a combined headline score,
    with each config's own score and ΔRA/ΔDec/ΔPA on its line. The
    per-plan Apply fills ALL configs at once. Config j is optimized
    against the cumulative budget of configs 1..j-1 (each at its own #1),
    so plan #1 is the consistent recommended set; lower rows pair the
    k-th best of each config independently.
    """
    # Keep the leading run of configs that produced solutions. Under the
    # monotonic per-pass budget, once one config finds nothing the later
    # ones are empty too, so this keeps config index == display position.
    results: list = []
    for r in results_list:
        if len(r.get("score", [])) == 0:
            break
        results.append(r)
    if len(results) == 0:
        _opt_render_results_in_modal(
            results_list[0] if results_list else {"score": []},
            ra_ref, dec_ref, pa_ref, n_sources, method=method)
        return
    if len(results) == 1:
        # Only Config 1 found targets — show the single-config table.
        _opt_render_results_in_modal(results[0], ra_ref, dec_ref, pa_ref,
                                     n_sources, method=method)
        if len(results_list) > 1:
            _set_status(
                "Configs beyond #1 found no additional targets under the "
                "current max-configs cap — showing Config 1 only.",
                "warn", clear_after=12,
            )
        return
    n_cfg = len(results)

    cos_dec = np.cos(np.deg2rad(dec_ref))
    is_hier = method == "Hierarchy"
    is_merit = method == "Meritocracy"

    def _parts(res: dict, k: int) -> dict:
        """Per-config stats for solution k: count, Σweight, a display
        label (method-aware, with the −drop suffix), and Δ offsets."""
        totals = res.get("total_count")
        sweights = res.get("sum_weight")
        breakdowns = res.get("tier_breakdown")
        ndrop = res.get("n_dropped")
        top_targets = res.get("top_targets")
        drop_reasons = res.get("drop_reasons")
        cnt = (int(totals[k]) if totals is not None and k < len(totals)
               else int(round(float(res["score"][k]))))
        sw = (float(sweights[k]) if sweights is not None
              and k < len(sweights) else None)
        drop = int(ndrop[k]) if ndrop is not None and k < len(ndrop) else 0
        if is_merit:
            lbl = f"Σw {sw:.0f} ({cnt})" if sw is not None else f"{cnt}"
        elif is_hier and breakdowns is not None and k < len(breakdowns):
            # Visible label shows only tiers that placed a source; the
            # full breakdown (incl. empty tiers) goes in the hover tooltip
            # below, so the cell stays short even with many tiers.
            tb = " · ".join(f"P{int(t)}:{c}"
                            for t, c in breakdowns[k] if c > 0)
            lbl = f"{tb} ({cnt})" if tb else f"{cnt}"
        else:
            lbl = f"{cnt}"
        if drop > 0:
            lbl += f" −{drop}"
        # Hover tooltip: top placed sources + the −drop breakdown by
        # reason (budget / collision / spectral). Mirrors the single-
        # config table so the prominent "−K" on Config 2 is explained.
        tip_lines: list = []
        # Full per-tier breakdown (incl. empty tiers) for Hierarchy, so a
        # clipped Score cell never hides a tier.
        if is_hier and breakdowns is not None and k < len(breakdowns):
            full_tb = " · ".join(f"P{int(t)}:{c}" for t, c in breakdowns[k])
            if full_tb:
                tip_lines.append(f"Tiers: {full_tb}")
                tip_lines.append("")
        tt = (top_targets[k] if top_targets is not None
              and k < len(top_targets) else None)
        if tt:
            tip_lines.append(f"Top {len(tt)} placed sources at this pointing:")
            for r, t in enumerate(tt, 1):
                pid = "—" if t["p"] is None else f"{int(t['p'])}"
                wid = "—" if t["w"] is None else f"{int(t['w'])}"
                tip_lines.append(f"{r}. ID={t['id']}  P={pid}  W={wid}")
        if drop > 0:
            reasons_i = (drop_reasons[k] if drop_reasons is not None
                         and k < len(drop_reasons) else None)
            if tip_lines:
                tip_lines.append("")
            tip_lines.extend(_opt_drop_breakdown_lines(reasons_i, drop))
        tooltip = "\n".join(tip_lines)
        if tooltip:
            title_attr = (tooltip.replace("&", "&amp;")
                          .replace('"', "&quot;").replace("<", "&lt;"))
            score_html = (
                f'<span title="{title_attr}" style="cursor:help; '
                f'border-bottom:1px dotted #888">{lbl}</span>'
            )
        else:
            score_html = lbl
        ra_i = float(res["ra"][k])
        dec_i = float(res["dec"][k])
        pa_i = float(res["pa"][k])
        return dict(
            cnt=cnt, sw=(sw if sw is not None else float(cnt)), lbl=lbl,
            score_html=score_html,
            dra=(ra_i - ra_ref) * 3600.0 * cos_dec,
            ddec=(dec_i - dec_ref) * 3600.0,
            dpa=(pa_i - pa_ref + 180.0) % 360.0 - 180.0,
        )

    def _combined(parts: list) -> str:
        """Combined headline across all configs (disjoint under the
        max-configs cap, so counts/weights simply add)."""
        cnt = sum(p["cnt"] for p in parts)
        if is_merit:
            return f"Σw {sum(p['sw'] for p in parts):.0f} ({cnt})"
        return f"{cnt}"

    p_first = [_parts(r, 0) for r in results]
    # Legend for the "−K" suffix — shown only when some plan actually
    # drops sources (usually later configs, where the max-configs cap
    # skips sources earlier configs already claimed). Hover a Score for
    # the exact per-reason breakdown.
    n_plans_pre = min(min(len(r.get("score", [])) for r in results), 10)
    _any_drop = any(
        any(_parts(r, k)["lbl"].find("−") >= 0 for r in results)
        for k in range(n_plans_pre)
    )
    drop_legend = (
        " <span style='color:#8a5a00'>“<b>−K</b>” = sources that landed in "
        "shutters but were skipped (usually because another config already "
        "observed them under the max-configs cap; hover a Score for the "
        "breakdown).</span>"
        if _any_drop else ""
    )
    opt_modal_results_summary.text = (
        f"<div style='font-size:12px; margin-bottom:4px'>"
        f"<b>{n_cfg}-config plan</b> · {n_sources} candidate sources · "
        f"<b>Method:</b> {method}. Each row is a plan "
        f"(<b>Config 1 #k + … + Config {n_cfg} #k</b>) — combined score on "
        f"the left, each config's own score + offsets on its line. "
        f"<b>Recommended: plan #1</b> = <b>{_combined(p_first)}</b> combined. "
        f"Each config is optimized against the cumulative budget of the "
        f"earlier ones, so plan #1 is the consistent set. Offsets are from "
        f"the search centre.{drop_legend}</div>"
    )

    # Banner: one-click apply of the recommended plan (#1).
    apply_both_btn = Button(
        label=f"✔ Apply recommended plan #1  ({_combined(p_first)} combined)",
        button_type="primary", width=380, height=30,
    )
    apply_both_btn.js_on_click(CustomJS(
        args=dict(trig=opt_apply_both_trigger),
        code=f"""
        if (!window.confirm("Apply the recommended {n_cfg}-config plan (#1)?\\n\\n"
            + "This CLEARS and refills all {n_cfg} configs.")) {{
            return;
        }}
        trig.value = "0," + Date.now() + "_" + Math.random();
        """,
    ))

    HEADER_BG = "#eef2f8"
    ROW_H = 24
    score_w = 250 if is_hier else 150 if is_merit else 72
    W = dict(rank=30, comb=98, cfg=34, dra=80, ddec=80, dpa=74,
             score=score_w, apply=150)
    total_w = sum(W.values())
    _base = {
        "padding": "0 6px", "line-height": f"{ROW_H}px",
        "box-sizing": "border-box", "height": f"{ROW_H}px",
        "overflow": "hidden", "white-space": "nowrap",
        "text-overflow": "ellipsis",
    }

    def _cell(text, w, *, header=False, bold=False, mono=False,
              align="center", bg=None):
        st = dict(_base)
        st["text-align"] = align
        if header:
            st.update({"background": HEADER_BG, "font-weight": "600",
                       "text-align": "center"})
        if bold:
            st["font-weight"] = "600"
        if mono:
            st["font-family"] = "monospace"
        if bg:
            st["background"] = bg
        return Div(text=text, width=w, height=ROW_H, styles=st)

    header_row = row(
        _cell("#", W["rank"], header=True),
        _cell("Combined", W["comb"], header=True),
        _cell("Cfg", W["cfg"], header=True),
        _cell("ΔRA", W["dra"], header=True),
        _cell("ΔDec", W["ddec"], header=True),
        _cell("ΔPA", W["dpa"], header=True),
        _cell("Score", W["score"], header=True),
        _cell("", W["apply"], header=True),
        spacing=0, width=total_w,
    )

    n_plans = min(min(len(r["score"]) for r in results), 10)
    plan_blocks: list = []
    for k in range(n_plans):
        parts_k = [_parts(r, k) for r in results]
        bg = "#f7f9fc" if k % 2 else "#ffffff"
        apply_btn = Button(
            label=f"Apply plan #{k + 1}" if k == 0 else f"Plan #{k + 1}",
            button_type="primary" if k == 0 else "default",
            width=W["apply"] - 12, height=ROW_H,
        )
        apply_btn.js_on_click(CustomJS(
            args=dict(trig=opt_apply_both_trigger),
            code=f"""
            if (!window.confirm("Apply {n_cfg}-config plan #{k + 1}?\\n\\n"
                + "This CLEARS and refills all {n_cfg} configs.")) {{
                return;
            }}
            trig.value = "{k}," + Date.now() + "_" + Math.random();
            """,
        ))
        # One line per config; the first carries the rank, combined score
        # and the per-plan Apply button, the rest just their own offsets.
        lines: list = []
        for ci, pk in enumerate(parts_k):
            cfg_cell = _cell(
                f"<b style='color:{_config_color(ci)}'>C{ci + 1}</b>",
                W["cfg"], bg=bg)
            offsets = (
                _cell(f"{pk['dra']:+.2f}″", W["dra"], mono=True,
                      align="right", bg=bg),
                _cell(f"{pk['ddec']:+.2f}″", W["ddec"], mono=True,
                      align="right", bg=bg),
                _cell(f"{pk['dpa']:+.3f}°", W["dpa"], mono=True,
                      align="right", bg=bg),
                _cell(pk["score_html"], W["score"], bg=bg),
            )
            if ci == 0:
                lines.append(row(
                    _cell(str(k + 1), W["rank"], bold=True, bg=bg),
                    _cell(_combined(parts_k), W["comb"], bold=True, bg=bg),
                    cfg_cell, *offsets, apply_btn,
                    spacing=0, width=total_w,
                ))
            else:
                lines.append(row(
                    _cell("", W["rank"], bg=bg),
                    _cell("", W["comb"], bg=bg),
                    cfg_cell, *offsets,
                    _cell("", W["apply"], bg=bg),
                    spacing=0, width=total_w,
                ))
        plan_blocks.append(column(
            *lines, spacing=0, width=total_w,
            styles={"border-bottom": "2px solid #d8dee8"},
        ))

    opt_modal_progress_box.visible = False
    opt_modal_results_box.visible = True
    opt_modal_results_rows.children = [
        row(apply_both_btn), header_row, *plan_blocks,
    ]


# ── Chunked optimizer runner ─────────────────────────────────────────────


_OPT_GRID_CHUNK = 400      # ~400 ms per chunk → bar updates ~3 ×/s
_OPT_DE_PER_TICK = 1       # one DE refinement per tick


def _opt_drive() -> None:
    """One state-machine step. Re-schedules itself until done."""
    if not _opt_run:
        return  # cancelled / wiped
    try:
        if _opt_run["phase"] == "grid":
            _opt_grid_step()
        elif _opt_run["phase"] == "hierarchy":
            _opt_hierarchy_step()
        elif _opt_run["phase"] == "de":
            _opt_de_step()
        elif _opt_run["phase"] == "done":
            return
    except Exception as e:  # noqa: BLE001
        _set_status(f"Optimizer failed: {e}", "err")
        traceback.print_exc()
        _opt_run.clear()
        _opt_hide_modal()


def _opt_grid_step() -> None:
    """Process the next batch of grid pointings, update progress, and
    either schedule another chunk or transition to DE refinement."""
    ev = _opt_run["evaluator"]
    weights = _opt_run["weights"]
    use_flux = (_opt_run["objective"] == "flux")
    ras = _opt_run["ra_cube"]
    decs = _opt_run["dec_cube"]
    pas = _opt_run["pa_cube"]
    scores = _opt_run["grid_scores"]
    idx = _opt_run["grid_idx"]
    n_total = len(ras)
    end = min(idx + _OPT_GRID_CHUNK, n_total)
    for i in range(idx, end):
        det, tp, _ = ev.evaluate(ras[i], decs[i], pas[i])
        if use_flux:
            scores[i] = float(np.sum(tp * ev.flux * weights))
        else:
            scores[i] = float(np.sum(det * weights))
    _opt_run["grid_idx"] = end

    elapsed = _now() - _opt_run["started"]
    rate = end / max(0.01, elapsed)
    eta = (n_total - end) / rate if rate > 0 else 0
    _opt_update_progress(
        f"Grid: {end:,} / {n_total:,} pointings evaluated · "
        f"{elapsed:.1f}s elapsed · ~{eta:.1f}s left",
        end / max(1, n_total) * 0.85,   # leave 15 % for DE
    )

    if end < n_total:
        curdoc().add_next_tick_callback(_opt_drive)
        return

    # Grid phase done → sort.
    order = np.argsort(-scores)
    grid_res = {
        "score": scores[order],
        "ra": ras[order], "dec": decs[order], "pa": pas[order],
    }
    n_top = _opt_run["n_top"]
    _opt_run["grid_result"] = grid_res
    # Hierarchy: prepare the multi-stage filter. We keep a pool of up
    # to 100 grid candidates and rank them lexicographically by
    # tier-by-tier score from highest priority (smallest p) downward.
    if _opt_run.get("method") == "Hierarchy":
        K = min(100, len(scores))
        _opt_run["hier_pool"] = list(range(K))   # indices into grid_res
        _opt_run["hier_scores"] = {}             # per-(cand, tier) score
        # Tier list: ascending priority value = descending importance.
        pri = _opt_run["priorities"]
        tiers = sorted(set(float(p) for p in pri if np.isfinite(p)))
        _opt_run["hier_tiers"] = tiers
        _opt_run["hier_tier_idx"] = 0
        _opt_run["phase"] = "hierarchy"
    else:
        _opt_run["phase"] = "de"
        _opt_run["de_idx"] = 0
        _opt_run["de_total"] = min(n_top, len(scores))
        _opt_run["de_scores"] = []
        _opt_run["de_params"] = []
    curdoc().add_next_tick_callback(_opt_drive)


def _opt_hierarchy_step() -> None:
    """One priority tier of the multi-stage hierarchy filter.

    Score each surviving candidate by `1` per source at the current
    tier; keep only those tied with the per-stage max. After the last
    tier, the survivors are passed to DE refinement.
    """
    ev = _opt_run["evaluator"]
    pri = _opt_run["priorities"]
    grid_res = _opt_run["grid_result"]
    pool = _opt_run["hier_pool"]
    tiers = _opt_run["hier_tiers"]
    t_idx = _opt_run["hier_tier_idx"]

    if t_idx >= len(tiers) or not pool:
        # Filter done → DE refinement on the surviving pool.
        n_top = min(_opt_run["n_top"], len(pool))
        # Build a "fake" grid_result from the pool ordering so the
        # existing DE machinery can reuse its array layout.
        pool_arr = np.asarray(pool, dtype=int)
        _opt_run["grid_result"] = {
            "score": np.asarray(
                [_opt_run["hier_scores"].get((i, 0), 0.0) for i in pool],
                dtype=float,
            ),
            "ra": grid_res["ra"][pool_arr],
            "dec": grid_res["dec"][pool_arr],
            "pa": grid_res["pa"][pool_arr],
        }
        _opt_run["phase"] = "de"
        _opt_run["de_idx"] = 0
        _opt_run["de_total"] = n_top
        _opt_run["de_scores"] = []
        _opt_run["de_params"] = []
        # DE refinement weights must PRESERVE the multi-tier lex
        # ordering established by the filter — otherwise DE, given
        # only tier-0 weights, would happily slide to a pointing
        # with same tier-0 count but fewer tier-1 / tier-2 placements.
        # Use the auto-weight formula (smallest int weights s.t. each
        # tier strictly outweighs the sum of all lower tiers); their
        # sum is then a lex-equivalent scalar that DE can maximise
        # without violating the priority ordering.
        if tiers:
            pri_strs = [
                str(int(p)) if np.isfinite(p) else "" for p in pri
            ]
            w_strs = _compute_weights_from_priorities(pri_strs)
            if w_strs is None:
                lex_weights = np.ones_like(pri, dtype=float)
            else:
                lex_weights = np.array(
                    [float(s) if s else 0.0 for s in w_strs],
                    dtype=float,
                )
            _opt_run["weights"] = lex_weights
            # Store the tiers so post-DE we can compute per-tier
            # breakdowns for the results display.
            _opt_run["lex_tiers"] = tiers
        curdoc().add_next_tick_callback(_opt_drive)
        return

    tier = tiers[t_idx]
    tier_mask = (pri == tier)
    new_pool: list[int] = []
    scores_this_tier: list[float] = []
    for k in pool:
        ra_k = float(grid_res["ra"][k])
        dec_k = float(grid_res["dec"][k])
        pa_k = float(grid_res["pa"][k])
        det, _tp, _ = ev.evaluate(ra_k, dec_k, pa_k)
        s = float(np.sum(det[tier_mask])) if tier_mask.any() else 0.0
        _opt_run["hier_scores"][(k, t_idx)] = s
        scores_this_tier.append(s)
    # Keep only candidates that tie the per-stage max.
    if scores_this_tier:
        best = max(scores_this_tier)
        for k, s in zip(pool, scores_this_tier):
            if s >= best - 1e-9:
                new_pool.append(k)
    _opt_run["hier_pool"] = new_pool
    _opt_run["hier_tier_idx"] = t_idx + 1

    frac = 0.85 + 0.10 * ((t_idx + 1) / max(1, len(tiers)))
    elapsed = _now() - _opt_run["started"]
    _opt_update_progress(
        f"Hierarchy filter: tier {t_idx + 1} / {len(tiers)} "
        f"(p={tier:g}) — survivors: {len(new_pool)} · "
        f"{elapsed:.1f}s elapsed",
        frac,
    )
    curdoc().add_next_tick_callback(_opt_drive)


def _opt_de_step() -> None:
    """Refine one (or a few) top grid candidates via DE."""
    from vmpt.optimizer import refine_top
    de_idx = _opt_run["de_idx"]
    de_total = _opt_run["de_total"]
    if de_idx >= de_total:
        # Already done — final ranking + dedup pass via refine_top's
        # internal sort/dedup. Stitch our per-candidate refinements
        # into a fake grid_results and re-feed (no extra work — DE box
        # is zero so it terminates instantly).
        refined = {
            "score": np.asarray(_opt_run["de_scores"]),
            "ra": np.asarray([p[0] for p in _opt_run["de_params"]]),
            "dec": np.asarray([p[1] for p in _opt_run["de_params"]]),
            "pa": np.asarray([p[2] for p in _opt_run["de_params"]]),
        }
        order = np.argsort(-refined["score"])
        for k in ("score", "ra", "dec", "pa"):
            refined[k] = refined[k][order]

        # Manual dedup using same tolerances as refine_top default.
        ra_tol = 0.3 / 3600.0 / np.cos(np.deg2rad(_opt_run["dec_ref"]))
        dec_tol = 0.3 / 3600.0
        pa_tol = 0.05
        keep: list[int] = []
        for i in range(len(refined["score"])):
            dup = False
            for j in keep:
                if (abs(refined["ra"][i] - refined["ra"][j]) <= ra_tol
                        and abs(refined["dec"][i] - refined["dec"][j]) <= dec_tol
                        and abs(((refined["pa"][i] - refined["pa"][j]
                                 + 180) % 360) - 180) <= pa_tol):
                    dup = True
                    break
            if not dup:
                keep.append(i)
        refined = {k: refined[k][keep] for k in ("score", "ra", "dec", "pa")}

        # For Hierarchy mode, compute a per-tier source-count
        # breakdown for each surviving candidate (e.g. P0:4·P1:12·P2:30).
        # This is what the user actually wants to see — the headline
        # score from the lex-weighted DE objective is hard to
        # interpret on its own.
        # For every refined candidate: compute (a) the total count of
        # observable sources, (b) the top-10 sources at this pointing
        # sorted by priority ascending then weight descending. Both
        # feed the results-modal display — the count appears inline,
        # the top-10 list as a hover tooltip.
        ev = _opt_run["evaluator"]
        pri_arr = _opt_run["priorities"]
        wt_arr = _opt_run.get(
            "weight_arr",
            np.full(len(pri_arr), np.nan, dtype=float),
        )
        ids_arr = _opt_run.get("source_ids", [])
        protect_mask = _opt_run.get("protect_mask")
        totals: list[int] = []
        sum_weights: list[float] = []
        top_targets: list[list[dict]] = []
        n_dropped_list: list[int] = []
        reasons_list: list[dict] = []
        for k in range(len(refined["score"])):
            det, _, _, reasons = ev.evaluate_with_reasons(
                float(refined["ra"][k]),
                float(refined["dec"][k]),
                float(refined["pa"][k]),
            )
            totals.append(int(det.sum()))
            n_dropped_list.append(int(sum(reasons.values())))
            reasons_list.append(dict(reasons))
            # Weight sum across detected sources (Meritocracy headline).
            w_finite = np.where(np.isfinite(wt_arr), wt_arr, 0.0)
            sum_weights.append(float(np.sum(det * w_finite)))
            # Top-10 detected sources, sorted by priority asc (NaN last)
            # then weight desc. The list goes into a `title=` tooltip.
            idxs = np.where(det)[0]
            if idxs.size:
                pri_key = np.where(np.isfinite(pri_arr[idxs]),
                                    pri_arr[idxs], np.inf)
                wt_key = np.where(np.isfinite(wt_arr[idxs]),
                                   wt_arr[idxs], 0.0)
                # Lexicographic sort: priority asc, weight desc.
                order = np.lexsort((-wt_key, pri_key))
                top_idx = idxs[order[:10]]
                top_targets.append([
                    {
                        "id": (str(ids_arr[i]) if i < len(ids_arr)
                               else f"<row {i}>"),
                        "p": (float(pri_arr[i])
                              if np.isfinite(pri_arr[i]) else None),
                        "w": (float(wt_arr[i])
                              if np.isfinite(wt_arr[i]) else None),
                        "prot": (bool(protect_mask[i])
                                 if protect_mask is not None
                                 and i < len(protect_mask) else False),
                    }
                    for i in top_idx
                ])
            else:
                top_targets.append([])
        refined["total_count"] = totals
        refined["sum_weight"] = sum_weights
        refined["top_targets"] = top_targets
        refined["n_dropped"] = n_dropped_list
        # Per-pointing per-reason drop counts. Dict shape:
        # {"collision": int, "required_lam": int, "no_gap": int,
        #  "extend_blue": int, "extend_red": int}.  The results
        # modal renders the breakdown in the Score-cell tooltip.
        refined["drop_reasons"] = reasons_list
        refined["protect_enabled"] = bool(_opt_run.get("protect_enabled"))

        # Hierarchy mode also wants the per-tier breakdown.
        breakdowns: list[list[tuple[float, int]]] | None = None
        if _opt_run.get("method") == "Hierarchy" and _opt_run.get("lex_tiers"):
            tiers_list = _opt_run["lex_tiers"]
            breakdowns = []
            for k in range(len(refined["score"])):
                det, _, _ = ev.evaluate(
                    float(refined["ra"][k]),
                    float(refined["dec"][k]),
                    float(refined["pa"][k]),
                )
                per_tier: list[tuple[float, int]] = []
                for t in tiers_list:
                    cnt = int(np.sum(det & (pri_arr == t)))
                    per_tier.append((float(t), cnt))
                breakdowns.append(per_tier)
            refined["tier_breakdown"] = breakdowns

        # ── Multi-config (v1.4.0; generalised to N passes in v1.5.0) ──
        # Sequential greedy: after each config's DE finishes, charge its
        # best pointing's observed sources against their max-configs cap
        # and, while configs remain, launch a fresh search for the next
        # one on the reduced budget. Each refined result is stashed in
        # `pass_results` (index i = config i).
        n_pass = int(_opt_run.get("n_pass", 1))
        cur = int(_opt_run.get("pass", 1))            # 1-based pass number
        if n_pass >= 2:
            _opt_run.setdefault("pass_results", []).append(dict(refined))

        if n_pass >= 2 and cur < n_pass:
            ev = _opt_run["evaluator"]
            # Sources observed by this config's best pick. `ev` still
            # carries the budget this config ran under, so `det_best` is
            # already disjoint from earlier configs except where a source's
            # max_configs cap permits re-observation.
            try:
                det_best, _, _ = ev.evaluate(
                    float(refined["ra"][0]), float(refined["dec"][0]),
                    float(refined["pa"][0]))
                observed_now = det_best.astype(int)
            except Exception:  # noqa: BLE001
                observed_now = np.zeros(
                    len(_opt_run["effective_max"]), dtype=int)
            observed_total = _opt_run.get("observed_total")
            if observed_total is None:
                observed_total = np.zeros_like(observed_now)
            observed_total = observed_total + observed_now
            _opt_run["observed_total"] = observed_total
            eff = np.asarray(_opt_run["effective_max"], dtype=float)
            budget_next = np.asarray(eff > observed_total, dtype=bool)
            ev._budget = budget_next
            ev._budget_enabled = not bool(budget_next.all())
            # budgets[config_index]: the next config is index `cur` (0-based).
            _opt_run["budgets"][cur] = budget_next.copy()
            # Reset the state machine for a fresh search of the next config.
            _opt_run["pass"] = cur + 1
            _opt_run["phase"] = "grid"
            _opt_run["grid_idx"] = 0
            _opt_run["grid_scores"] = np.zeros(
                int(_opt_run["n_total"]), dtype=float)
            _opt_run["de_idx"] = 0
            _opt_run["de_scores"] = []
            _opt_run["de_params"] = []
            for _k in ("grid_result", "hier_pool", "hier_scores",
                       "hier_tiers", "hier_tier_idx", "lex_tiers"):
                _opt_run.pop(_k, None)
            _opt_run["started"] = _now()
            _opt_update_progress(
                f"Config {cur + 1} of {n_pass}: optimizing on the "
                f"remaining budget…", 0.0)
            curdoc().add_next_tick_callback(_opt_drive)
            return

        if n_pass >= 2:
            _opt_render_multi_results(
                _opt_run.get("pass_results") or [refined],
                _opt_run["ra_ref"], _opt_run["dec_ref"], _opt_run["pa_ref"],
                _opt_run["n_sources"],
                method=_opt_run.get("method", "Democracy"),
            )
            _set_status(
                f"Optimization complete: {n_pass}-config plan ready. Apply "
                f"a row or the recommended plan.", "ok", clear_after=12,
            )
        else:
            _opt_render_results_in_modal(
                refined,
                _opt_run["ra_ref"], _opt_run["dec_ref"], _opt_run["pa_ref"],
                _opt_run["n_sources"],
                method=_opt_run.get("method", "Democracy"),
            )
            if breakdowns and breakdowns[0]:
                best_summary = " · ".join(
                    f"P{int(t)}={c}" for t, c in breakdowns[0]
                )
                _set_status(
                    f"Optimization complete (Hierarchy): "
                    f"{best_summary} of {_opt_run['n_sources']} sources.",
                    "ok", clear_after=12,
                )
            else:
                _set_status(
                    f"Optimization complete: best score "
                    f"{refined['score'][0]:.1f} of "
                    f"{_opt_run['n_sources']} sources.",
                    "ok", clear_after=10,
                )
        _opt_run["phase"] = "done"
        return

    # Refine one candidate. Re-uses optimizer.refine_top with n_top=1
    # so we get the same dedup-aware bound-aware logic for free.
    grid_res = _opt_run["grid_result"]
    single = {
        "score": grid_res["score"][de_idx:de_idx + 1],
        "ra": grid_res["ra"][de_idx:de_idx + 1],
        "dec": grid_res["dec"][de_idx:de_idx + 1],
        "pa": grid_res["pa"][de_idx:de_idx + 1],
    }
    refined = refine_top(
        _opt_run["evaluator"], single,
        n_top=1,
        dra_arcsec=_opt_run["de_dra_arcsec"],
        ddec_arcsec=_opt_run["de_ddec_arcsec"],
        dpa_deg=_opt_run["de_dpa_deg"],
        maxiter=_opt_run["maxiter"],
        weights=_opt_run["weights"], objective=_opt_run["objective"],
    )
    _opt_run["de_scores"].append(float(refined["score"][0]))
    _opt_run["de_params"].append((
        float(refined["ra"][0]),
        float(refined["dec"][0]),
        float(refined["pa"][0]),
    ))
    _opt_run["de_idx"] = de_idx + 1

    elapsed = _now() - _opt_run["started"]
    frac = 0.85 + 0.15 * ((de_idx + 1) / max(1, de_total))
    _opt_update_progress(
        f"Refining top {de_total}: {de_idx + 1} / {de_total} · "
        f"{elapsed:.1f}s elapsed",
        frac,
    )
    curdoc().add_next_tick_callback(_opt_drive)


# ---------------------------------------------------------------------
# Collision protection — helpers
# ---------------------------------------------------------------------


def _protect_mask_for_catalog(
    cat,
    enabled: bool,
    mode_idx: int,
    threshold_text: str,
) -> tuple[np.ndarray | None, str | None]:
    """Build a per-source boolean mask of "protected" sources.

    Parameters mirror the UI controls. ``mode_idx`` is the RadioGroup's
    ``active`` index: 0 → "By priority ≤", 1 → "By weight ≥". Returns
    ``(mask, error_msg)``. When protection is disabled or impossible
    (no catalog, blank threshold, …), returns ``(None, None)``. When
    the user's threshold cannot be parsed or the requested column is
    missing, returns ``(None, error_msg)`` so the caller can surface a
    friendly status line.
    """
    if not enabled:
        return None, None
    if cat is None or len(cat.ra_deg) == 0:
        return None, "Load a catalog first."
    text = (threshold_text or "").strip()
    if not text:
        return None, "Enter a threshold."
    try:
        threshold = float(text)
    except ValueError:
        return None, "Threshold must be numeric."
    if mode_idx == 0:
        pri = np.asarray(cat.priority, dtype=float)
        if not np.isfinite(pri).any():
            return None, "Catalog has no priority values."
        mask = np.isfinite(pri) & (pri <= threshold)
    else:
        wgt = np.asarray(getattr(cat, "weight", []), dtype=float)
        if wgt.size != len(cat.ra_deg) or not np.isfinite(wgt).any():
            return None, "Catalog has no weight values."
        mask = np.isfinite(wgt) & (wgt >= threshold)
    return mask, None


def _update_protect_status_div() -> None:
    """Refresh the protect-status Div based on the current widget state.

    Shows the count of sources that would be marked protected at the
    current threshold, or a small warning when the configuration is
    invalid (missing column, blank threshold, etc.). Wired to the
    relevant widgets' on_change events and to catalog-change events
    via the Tabs container so it updates live."""
    enabled = bool(opt_protect_enable_cb.active)
    if not enabled:
        opt_protect_status_div.text = (
            "<small style='color:#5a6b85'>"
            "Disabled — all detected targets count toward the score."
            "</small>"
        )
        return
    cat = state.get("catalog")
    mode_idx = int(opt_protect_mode_radio.active or 0)
    mask, err = _protect_mask_for_catalog(
        cat, enabled, mode_idx, opt_protect_threshold_input.value,
    )
    if err is not None:
        opt_protect_status_div.text = (
            f"<small style='color:#a05a30'>⚠ {err}</small>"
        )
        return
    n_prot = int(mask.sum())
    n_total = int(mask.size)
    n_other = n_total - n_prot
    if n_prot == 0:
        opt_protect_status_div.text = (
            "<small style='color:#a05a30'>"
            "⚠ Threshold matches no catalog source.</small>"
        )
        return
    # Add a one-line hint about the current Disperser / Filter — H
    # gratings span ~500″ in V2 and therefore drop nearly every
    # co-observable source; PRISM's ~35″ is benign.
    disp = (state.get("disperser") or "?").upper()
    filt = (state.get("filter") or "?").upper()
    try:
        from vmpt.wavelengths import v2_overlap_distance as _v2od
        overlap = _v2od(disp, filt)
    except Exception:
        overlap = None
    overlap_txt = (f" · V2 overlap ≈ {overlap:.0f}″"
                   if overlap is not None else "")
    opt_protect_status_div.text = (
        f"<small style='color:#1a3b66'>"
        f"<b>{n_prot}</b> protected · {n_other} other "
        f"<span style='color:#5a6b85'>"
        f"({disp} / {filt}{overlap_txt})</span></small>"
    )


# Wire the live-update.
opt_protect_enable_cb.on_change(
    "active", lambda attr, old, new: _update_protect_status_div(),
)
opt_protect_mode_radio.on_change(
    "active", lambda attr, old, new: _update_protect_status_div(),
)
opt_protect_threshold_input.on_change(
    "value", lambda attr, old, new: _update_protect_status_div(),
)


def _now() -> float:
    """Wall-clock seconds since epoch (relative wall-time is fine)."""
    import time as _time
    return _time.time()


def on_optimize():
    """Validate inputs and start the chunked optimization run.

    All work happens via `_opt_drive` ticks so the IO loop keeps
    rendering — the user sees the progress bar advance instead of
    a frozen UI."""
    cat = state.get("catalog")
    if cat is None or len(cat.ra_deg) == 0:
        _set_status("Load a catalog first.", "warn")
        return

    try:
        d_ra = float(opt_dra_input.value)
        d_dec = float(opt_ddec_input.value)
        d_pa = float(opt_dpa_input.value)
        n_top = int(float(opt_n_top_input.value))
    except (TypeError, ValueError):
        _set_status("Optimizer: ΔRA/ΔDec/ΔPA/N must be numeric.", "err")
        return
    try:
        n_ra = int(float(opt_grid_n_ra_input.value))
        n_dec = int(float(opt_grid_n_dec_input.value))
        n_pa = int(float(opt_grid_n_pa_input.value))
        maxiter = int(float(opt_de_maxiter_input.value))
        sigma = float(opt_sigma_input.value)
    except (TypeError, ValueError):
        _set_status("Optimizer: advanced numeric inputs invalid.", "err")
        return

    centration = opt_centration_select.value
    objective = opt_objective_select.value
    method = opt_method_select.value or "Democracy"

    ra0 = float(state.get("ra_deg") or 0.0)
    dec0 = float(state.get("dec_deg") or 0.0)
    pa0 = float(state.get("pa_v3") or 0.0)

    # Filter catalog by priority cutoff if set.
    ra_arr = np.asarray(cat.ra_deg, dtype=float)
    dec_arr = np.asarray(cat.dec_deg, dtype=float)
    pri = np.asarray(cat.priority, dtype=float)
    # Source IDs as strings, parallel to ra_arr. Used at Apply time
    # to tag each opened slitlet with the right catalog-source ID.
    ids_arr = np.asarray([str(v) for v in np.asarray(cat.ids).tolist()],
                         dtype=object)
    weight_arr_full = np.asarray(getattr(cat, "weight", []), dtype=float)
    if weight_arr_full.size != len(ra_arr):
        weight_arr_full = np.full(len(ra_arr), np.nan, dtype=float)
    pri_cut_text = (opt_priority_input.value or "").strip()
    keep = None
    if pri_cut_text:
        try:
            cutoff = float(pri_cut_text)
            keep = np.where(np.isnan(pri), False, pri <= cutoff)
            if int(keep.sum()) == 0:
                _set_status(
                    f"Priority cutoff ≤ {cutoff} excludes every source.",
                    "err",
                )
                return
            ra_arr = ra_arr[keep]
            dec_arr = dec_arr[keep]
            pri = pri[keep]
            ids_arr = ids_arr[keep]
            weight_arr_full = weight_arr_full[keep]
        except ValueError:
            _set_status("Priority cutoff must be numeric.", "err")
            return

    # --- Method-specific weight array + validation ---
    # Democracy: every source counts 1, NaN-priority sources still in.
    # Meritocracy: weight column required; NaN → 0 (silent skip).
    # Hierarchy: priority column required; multi-stage filter built later.
    if method == "Meritocracy":
        if not np.isfinite(weight_arr_full).any():
            _set_status(
                "Meritocracy needs a Weight column. Use the catalog "
                "editor to add weights (or `Compute w from p`).", "err",
            )
            return
        weights = np.where(np.isfinite(weight_arr_full),
                           weight_arr_full, 0.0)
    elif method == "Hierarchy":
        if not np.isfinite(pri).any():
            _set_status(
                "Hierarchy needs a Priority column. Use the catalog "
                "editor to fill priorities.", "err",
            )
            return
        # Grid phase uses uniform weights so we cover candidates that
        # serve any tier; the multi-stage filter applies the lex order.
        weights = np.ones_like(ra_arr)
    else:  # Democracy
        weights = np.ones_like(ra_arr)

    flux_arr = None
    if objective == "flux":
        mag = np.asarray(cat.mag, dtype=float)
        if keep is not None:
            mag = mag[keep]
        flux_arr = np.where(np.isfinite(mag),
                            10.0 ** (-0.4 * mag), 1.0)

    # Honour ΔX=0 as "freeze that axis" — same convention as the
    # underlying optimizer module.
    if d_ra <= 0:
        n_ra = 1
    if d_dec <= 0:
        n_dec = 1
    if d_pa <= 0:
        n_pa = 1
    n_total = n_ra * n_dec * n_pa
    n_sources = len(ra_arr)

    # Build the (ra, dec, pa) cube.
    cos_dec = max(np.cos(np.deg2rad(dec0)), 1e-3)
    dra_deg = max(d_ra, 0.0) / 3600.0 / cos_dec
    ddec_deg = max(d_dec, 0.0) / 3600.0
    ras = ra0 + (np.array([0.0]) if n_ra == 1
                 else np.linspace(-dra_deg, dra_deg, n_ra))
    decs = dec0 + (np.array([0.0]) if n_dec == 1
                   else np.linspace(-ddec_deg, ddec_deg, n_dec))
    pas_grid = pa0 + (np.array([0.0]) if n_pa == 1
                      else np.linspace(-d_pa, d_pa, n_pa))
    R, D, P = np.meshgrid(ras, decs, pas_grid, indexing="ij")
    ra_cube = R.ravel()
    dec_cube = D.ravel()
    pa_cube = P.ravel()

    # ---- Collision-protection mask (post-priority-cutoff slice) ----
    protect_enabled = bool(opt_protect_enable_cb.active)
    protect_mask_evaluator = None
    protect_mode_idx = int(opt_protect_mode_radio.active or 0)
    if protect_enabled:
        # Build over the FULL catalog first, then slice with `keep` so
        # indices align with `ra_arr / dec_arr / weight_arr_full`.
        full_mask, err = _protect_mask_for_catalog(
            cat, True, protect_mode_idx,
            opt_protect_threshold_input.value,
        )
        if err is not None or full_mask is None:
            _set_status(
                f"Collision protection: {err or 'invalid configuration'}.",
                "err",
            )
            return
        if keep is not None:
            full_mask = full_mask[keep]
        if not full_mask.any():
            _set_status(
                "Collision protection: threshold matches no participating "
                "source (after priority cutoff).", "warn",
            )
            return
        protect_mask_evaluator = full_mask

    # ---- Per-target spectral constraints (v1.3.0+) ----
    # Pull the five per-target arrays off the merged catalog, slice
    # them by the priority-cutoff `keep` mask if set, and pass through
    # to the evaluator. All optional — defaults preserve v1.2.x
    # behaviour when the user hasn't set any constraints.
    n_full = len(cat.ra_deg)
    def _slice_or_default(arr, default, dtype):
        a = np.asarray(arr if arr is not None and len(arr) == n_full
                       else [default] * n_full, dtype=dtype)
        return a[keep] if keep is not None else a

    cat_required_lam = getattr(cat, "required_lam", None)
    if cat_required_lam is not None and len(cat_required_lam) == n_full:
        full_req = np.asarray(cat_required_lam, dtype=object)
    else:
        # Build a true 1D object array of [] (not the 2D
        # shape=(n,0) array that np.array([[],[],…]) yields).
        full_req = np.empty(n_full, dtype=object)
        for _i in range(n_full):
            full_req[_i] = []
    required_lam_arr = (full_req[keep] if keep is not None else full_req)
    no_gap_arr = _slice_or_default(
        getattr(cat, "no_gap", None), False, bool)
    extend_blue_arr = _slice_or_default(
        getattr(cat, "extend_blue", None), False, bool)
    extend_red_arr = _slice_or_default(
        getattr(cat, "extend_red", None), False, bool)
    cat_protect_arr = _slice_or_default(
        getattr(cat, "protect", None), False, bool)
    # Per-target centration override (v1.3.1+). Strings (not bools),
    # so we can't reuse `_slice_or_default`. Treat missing / wrong-
    # length as all-empty (= use the global setting for every row).
    #
    # NB: `getattr(...) or []` triggers `bool(arr)` on a numpy array,
    # which raises `ValueError: The truth value of an array with more
    # than one element is ambiguous` whenever the catalog has 2+
    # sources. Check `is None` explicitly instead.
    _cent_raw = getattr(cat, "centration", None)
    if _cent_raw is None:
        _cent_full = np.array([], dtype=object)
    else:
        _cent_full = np.asarray(_cent_raw, dtype=object)
    if _cent_full.size != len(cat.ra_deg):
        _cent_full = np.array([""] * len(cat.ra_deg), dtype=object)
    centration_per_target_arr = (
        _cent_full[keep] if keep is not None else _cent_full
    )

    # ---- Multi-config observation budget (v1.4.0) ----
    # n_pass = how many configs the optimizer should plan (1 or 2). When
    # 2, it runs sequential-greedy: optimize config 1 (full budget), then
    # optimize config 2 with the config-1 best's observed sources charged
    # against their cap. effective_max[i] = per-source override (Catalog
    # max_configs) if set, else the optimizer-modal global default (blank
    # = unlimited → +inf).
    n_pass = max(1, min(_MAX_CONFIGS, int(mpt_num_configs_spinner.value or 1)))
    _gmax_text = (opt_global_max_configs_input.value or "").strip()
    try:
        global_max = float(_gmax_text) if _gmax_text else np.inf
    except ValueError:
        global_max = np.inf
    per_src_max = _slice_or_default(
        getattr(cat, "max_configs", None), np.nan, float)
    effective_max = np.where(np.isfinite(per_src_max), per_src_max,
                             global_max)

    # Build the evaluator now so the heavy CloughTocher Delaunay step
    # happens before we show the modal — the user sees the bar start
    # at 0 % and advance, instead of staring at a frozen 0 % for 2 s.
    from vmpt.optimizer import PointingEvaluator
    _set_status(
        "Building MSA inverse map (first time only)…",
        "info", clear_after=0,
    )
    # Load the reason mask alongside operability so the evaluator can
    # apply the "protected target on stuck-open row" rule. Loading is
    # cached by `app/msa.load_operability` after the first call.
    # Also load when per-target `protect` is set (per-row collision
    # protection needs the same stuck-open data).
    need_reason = (protect_mask_evaluator is not None
                   or bool(cat_protect_arr.any()))
    if need_reason:
        _op_mask, _op_reason = load_operability()
    else:
        _op_mask, _op_reason = None, None
    ev = PointingEvaluator(
        ra_arr, dec_arr, flux_sources=flux_arr,
        sigma_arcsec=sigma, centration=centration,
        slit_length=int(state.get("slitlet_height", 3)),
        # Operability is loaded fresh by PointingEvaluator from
        # `app/msa.load_operability()` — the same CRDS file the main
        # canvas uses. No state-side caching to drift between them.
        operable=_op_mask,
        protect_mask=protect_mask_evaluator,
        priorities=pri,
        weights=weight_arr_full,
        disperser=state.get("disperser"),
        filt=state.get("filter"),
        reason=_op_reason,
        # Per-target spectral constraints (v1.3.0+). The evaluator
        # ORs `protect` with `protect_mask` so either source flagging
        # a target as protected switches on the v1.2.x collision rules.
        required_lam=required_lam_arr,
        no_gap=no_gap_arr,
        extend_blue=extend_blue_arr,
        extend_red=extend_red_arr,
        protect=cat_protect_arr,
        # Per-target centration override (v1.3.1+). Wins unconditionally
        # over `centration=...` for any row with a non-empty entry.
        centration_per_target=centration_per_target_arr,
    )

    _opt_run.clear()
    _opt_run.update({
        "phase": "grid",
        "method": method,
        "evaluator": ev,
        # Catalog IDs parallel to evaluator.ra/dec — used by
        # `_apply_optimizer_result` to tag each opened slitlet with
        # its catalog source ID.
        "source_ids": ids_arr,
        "ra_cube": ra_cube, "dec_cube": dec_cube, "pa_cube": pa_cube,
        "grid_scores": np.zeros(n_total, dtype=float),
        "grid_idx": 0,
        "n_top": max(1, n_top),
        "weights": weights, "objective": objective,
        # Hierarchy needs the per-source priority array to compute
        # tier indicator weights during the filter phase.
        "priorities": pri,
        "weight_arr": weight_arr_full,
        "ra_ref": ra0, "dec_ref": dec0, "pa_ref": pa0,
        "n_sources": n_sources,
        "maxiter": maxiter,
        "de_dra_arcsec": max(2.0, d_ra / 10.0) if d_ra > 0 else 0.0,
        "de_ddec_arcsec": max(2.0, d_dec / 10.0) if d_dec > 0 else 0.0,
        "de_dpa_deg": max(0.5, d_pa / 10.0) if d_pa > 0 else 0.0,
        "started": _now(),
        # Collision-protection bookkeeping for the results modal +
        # _apply_optimizer_result. None when protection is disabled.
        "protect_enabled": protect_mask_evaluator is not None,
        "protect_mask": protect_mask_evaluator,
        "protect_mode_idx": protect_mode_idx,
        "protect_threshold": opt_protect_threshold_input.value,
        # ── Multi-config (v1.4.0; up to _MAX_CONFIGS passes in v1.5.0) ──
        "n_pass": n_pass,            # requested config count (1.._MAX_CONFIGS)
        "pass": 1,                   # current pass (1-based)
        "n_total": n_total,          # grid points (for per-pass grid reset)
        "effective_max": effective_max,
        "pass_results": [],          # one refined dict per config, in order
        "observed_total": None,      # cumulative observed count across passes
        # Per-config evaluator budgets used by Apply: config 0 → None (no
        # budget); config k → the cumulative budget after configs 0..k-1.
        "budgets": {0: None},
    })

    # Close the config modal now that we've kicked off the run; the
    # progress + results modal takes over the screen.
    _close_opt_config_modal()
    _opt_show_modal()
    _opt_update_progress(
        f"Grid: 0 / {n_total:,} pointings over {n_sources} sources…", 0.0,
    )
    _set_status(
        f"Optimization started: {n_total:,} grid pointings × "
        f"{n_sources} sources.", "info", clear_after=3,
    )
    curdoc().add_next_tick_callback(_opt_drive)


def _on_modal_close() -> None:
    """Close the modal mid-run or post-run. If a run is in progress
    we clear state too so chunks stop firing."""
    _opt_hide_modal()
    if _opt_run.get("phase") in ("grid", "hierarchy", "de"):
        _opt_run.clear()
        _set_status("Optimization cancelled.", "warn", clear_after=6)


opt_modal_close_btn.on_click(_on_modal_close)
opt_modal_top_close_btn.on_click(_on_modal_close)
opt_run_btn.on_click(on_optimize)


# ---------------------------------------------------------------------------
# Wiring
# ---------------------------------------------------------------------------

fits_path_input.on_change("value", on_fits_path)
jpg_path_input.on_change("value", on_jpg_path)
sidecar_path_input.on_change("value", on_sidecar_path)
catalog_path_input.on_change("value", on_catalog_path)
catalog_add_btn.on_click(on_catalog_add)
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


def _on_stats_bar_choice(attr, old, new):
    """Translate the picker's value list (display labels) back into
    the cell-key list that ``refresh_overlays`` reads from
    ``state["stats_bar_order"]`` and trigger a redraw.

    Note: must call the FULL :func:`refresh_overlays`, not the
    "light" variant — the stats bar HTML is built at the tail end
    of the full refresh (after all the per-shutter computation).
    The light variant only redraws the cheap glyphs (MSA outline,
    fixed slits, pointing handle) and does not touch ``stats_div``.
    The user only retoggles cells occasionally, so the extra cost
    is invisible in practice.
    """
    keys = [
        _STATS_BAR_LABEL_TO_KEY[label]
        for label in (new or [])
        if label in _STATS_BAR_LABEL_TO_KEY
    ]
    state["stats_bar_order"] = keys
    refresh_overlays()


stats_bar_choice.on_change("value", _on_stats_bar_choice)


def _reset_stats_bar_order() -> None:
    """Restore the canonical stats-bar order. Triggered by the
    "Reset to default" button in the Settings tab."""
    default_labels = [STATS_BAR_CELL_LABELS[k]
                      for k in STATS_BAR_DEFAULT_ORDER]
    stats_bar_choice.value = default_labels
    # `on_change` fires automatically — state + stats_div both update.


stats_bar_reset_btn = Button(
    label="Reset to default order",
    button_type="default",
    width=SIDEBAR_W - 20,
)
stats_bar_reset_btn.on_click(_reset_stats_bar_order)


def _refresh_catalog_hover_tooltip() -> None:
    """Rebuild ``catalog_hover.tooltips`` from the picker's value list.

    Each selected field contributes one of the template fragments
    from :data:`CATALOG_HOVER_FIELDS`, joined by a thin grey middle-
    dot separator. The result replaces the HoverTool's `.tooltips`
    in place (Bokeh propagates the change to the browser without
    re-creating the tool, so the user can keep their hover state).
    """
    keys = [
        _CATALOG_HOVER_LABEL_TO_KEY[label]
        for label in (catalog_hover_choice.value or [])
        if label in _CATALOG_HOVER_LABEL_TO_KEY
    ]
    if not keys:
        catalog_hover.tooltips = ""
        return
    fragments = [CATALOG_HOVER_FIELDS[k][1] for k in keys]
    sep = '<span style="color:#888;"> · </span>'
    catalog_hover.tooltips = (
        f'<div style="{_TIP_BASE_STYLE} color:#1a3b66;">'
        + sep.join(fragments)
        + '</div>'
    )


def _on_catalog_hover_choice(attr, old, new):
    _refresh_catalog_hover_tooltip()


catalog_hover_choice.on_change("value", _on_catalog_hover_choice)


def _reset_catalog_hover_order() -> None:
    catalog_hover_choice.value = [
        CATALOG_HOVER_FIELDS[k][0] for k in CATALOG_HOVER_DEFAULT_ORDER
    ]


catalog_hover_reset_btn = Button(
    label="Reset hover to default",
    button_type="default",
    width=SIDEBAR_W - 20,
)
catalog_hover_reset_btn.on_click(_reset_catalog_hover_order)

# Apply the initial tooltip so the hover matches the picker on first
# load (the picker's default value is CATALOG_HOVER_DEFAULT_ORDER,
# which maps to "id · ra/dec · priority" — same content as v1.0).
_refresh_catalog_hover_tooltip()


# ── Build the two "Customise…" modal cards (v1.3.0+) ─────────────────
# The cards reference the picker + reset-button widgets, so they have
# to be built AFTER those widgets are defined. Same layout pattern as
# `opt_config_modal_card` — title Div, content, Done button, with a
# top-right × dismiss and a fixed-position centered card.

_CUSTOMISE_MODAL_STYLES = {
    "position": "fixed",
    "top": "50%", "left": "50%",
    "transform": "translate(-50%, -50%)",
    "background": "white",
    "border": "1px solid #c0c8d6",
    "border-radius": "6px",
    "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
    "padding": "16px 18px",
    "z-index": "1000",
    "max-height": "92vh",
    "overflow-y": "auto",
}

# Bokeh injects ``body > div > .bk-Column { height: 100% !important; }``
# at the page root, which makes any top-level Column stretch to fill
# the viewport. That's correct for the main layout but wrong for our
# little customise modals — they have ~4 children and end up with a
# huge empty white area. Override it with a more-specific selector
# (matches both classes; specificity beats Bokeh's `body > div > class`)
# and re-assert it with !important so the picker modal hugs its content.
_CUSTOMISE_MODAL_CSS = GlobalInlineStyleSheet(css="""
.bk-Column.vmpt-customise-modal {
  height: auto !important;
}
""")

stats_bar_modal_card = column(
    row(
        Div(text="<h3>Customise top stats bar</h3>",
            sizing_mode="stretch_width"),
        stats_bar_modal_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    Div(text="<div style='font-size:12px; color:#5a6b85'>"
             "Pick which cells appear in the bar above the figure, "
             "and in what order. Drop a chip with its × to hide that "
             "cell; click in the dropdown to re-add. The order of "
             "the chips is the on-screen order.</div>",
        width=SIDEBAR_W + 100),
    stats_bar_choice,
    stats_bar_reset_btn,
    row(stats_bar_modal_close_btn, spacing=10),
    spacing=10,
    width=SIDEBAR_W + 130,
    visible=False,
    styles=_CUSTOMISE_MODAL_STYLES,
    css_classes=["vmpt-customise-modal", "vmpt-modal-card"],
    stylesheets=[_CUSTOMISE_MODAL_CSS],
)

catalog_hover_modal_card = column(
    row(
        Div(text="<h3>Customise catalog hover</h3>",
            sizing_mode="stretch_width"),
        catalog_hover_modal_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    Div(text="<div style='font-size:12px; color:#5a6b85'>"
             "Pick which fields show up in the tooltip when you "
             "hover a catalog target marker on the canvas, and in "
             "what order. The Constraints field renders a compact "
             "summary (e.g. <code>λ:1·G·🛡</code>) when the source "
             "has any per-target spectral constraint set.</div>",
        width=SIDEBAR_W + 100),
    catalog_hover_choice,
    catalog_hover_reset_btn,
    row(catalog_hover_modal_close_btn, spacing=10),
    spacing=10,
    width=SIDEBAR_W + 130,
    visible=False,
    styles=_CUSTOMISE_MODAL_STYLES,
    css_classes=["vmpt-customise-modal", "vmpt-modal-card"],
    stylesheets=[_CUSTOMISE_MODAL_CSS],
)


# ── Overwrite-confirmation modal (v1.3.3+) ───────────────────────────────
# Shown when a Save-as-CSV or Save-session would clobber an existing
# file. Lists the path(s) being overwritten, defaults to Cancel
# (safer), and only proceeds with the queued callback when the user
# explicitly clicks "Overwrite". For new paths the modal is bypassed
# entirely — same UX as before.
overwrite_modal_backdrop = Div(
    text="", width=0, height=0, visible=False,
    styles={
        "position": "fixed", "top": "0", "left": "0",
        "right": "0", "bottom": "0",
        "background": "rgba(20, 30, 50, 0.45)",
        # Stacks above the catalog-editor backdrop (z-index 999) so a
        # save from inside the editor still shows the confirmation
        # on top of the table.
        "z-index": "1020",
    },
)
overwrite_modal_top_close_btn = Button(
    label="×", button_type="default",
    width=32, height=28,
    css_classes=["vmpt-modal-x"],
)
overwrite_modal_body = Div(text="", width=520)
overwrite_modal_yes_btn = Button(
    label="Overwrite", button_type="danger", width=140,
)
overwrite_modal_no_btn = Button(
    label="Cancel", button_type="default", width=80,
)
overwrite_modal_card = column(
    row(
        Div(text="<h3 style='color:#a04030'>"
                 "Overwrite existing file?</h3>",
            sizing_mode="stretch_width"),
        overwrite_modal_top_close_btn,
        css_classes=["vmpt-modal-header"],
        styles=_MODAL_HEADER_STYLES,
        sizing_mode="stretch_width",
    ),
    overwrite_modal_body,
    row(overwrite_modal_yes_btn, overwrite_modal_no_btn, spacing=10),
    spacing=10,
    width=560,
    # As a fixed-position Bokeh root this column otherwise inherits a
    # full-viewport height (a near-empty 100vh box). "min" asks Bokeh to
    # shrink-wrap; the real enforcement is `_overwrite_modal_fit_js` below
    # (Bokeh's layout sets the height with !important, so only an inline
    # !important can shrink it). `max-height` + scroll is a safety net.
    height_policy="min",
    visible=False,
    # `vmpt-confirm-card` is a unique hook for the fit-to-content JS below
    # (the title lives in a nested shadow root, so we can't find this card
    # by its text — a dedicated class is the reliable handle).
    css_classes=["vmpt-modal-card", "vmpt-confirm-card"],
    styles={
        "position": "fixed",
        "top": "50%", "left": "50%",
        "transform": "translate(-50%, -50%)",
        "background": "white",
        "border": "1px solid #c0c8d6",
        "border-radius": "6px",
        "box-shadow": "0 10px 32px rgba(0, 30, 80, 0.3)",
        "padding": "16px 18px",
        "z-index": "1021",
        "max-height": "85vh",
        "overflow-y": "auto",
    },
)
# Bokeh's layout CSS forces the card to a full-viewport height with
# `!important`; a plain inline style or a document rule can't beat it, so
# pin the height to its content with an inline `!important` whenever the
# dialog opens (the only declaration that wins).
_overwrite_modal_fit_js = CustomJS(code="""
if (!cb_obj.visible) return;
const pin = () => {
    const out = [];
    (function w(r){
        r.querySelectorAll('.vmpt-confirm-card').forEach(e => out.push(e));
        r.querySelectorAll('*').forEach(el => { if (el.shadowRoot) w(el.shadowRoot); });
    })(document);
    if (out[0]) out[0].style.setProperty('height', 'fit-content', 'important');
};
pin();
setTimeout(pin, 30);
""")
overwrite_modal_card.js_on_change("visible", _overwrite_modal_fit_js)

# Stash the pending callback here while the user decides. Keying by a
# dict (not a closure variable) so the JS-free Python handlers stay
# composable.
_overwrite_pending: dict = {"callback": None}


def _confirm_overwrite_if_exists(paths, callback, what: str = "file") -> None:
    """If any of `paths` already exists, show the overwrite modal and
    only call `callback()` when the user confirms. Otherwise call
    `callback()` immediately.

    Parameters
    ----------
    paths : str | iterable[str]
        Path(s) that the save operation is about to write. Strings
        get expanded with ~ resolution.
    callback : callable
        Zero-argument function that performs the actual write.
    what : str
        Short label used in the modal (e.g. "CSV", "session JSON").
        Defaults to "file".
    """
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    expanded = [Path(p).expanduser() for p in paths]
    existing = [p for p in expanded if p.exists()]
    if not existing:
        callback()
        return
    if len(existing) == 1:
        body = (
            f"<p>The {what} already exists at:</p>"
            f"<p><code>{existing[0]}</code></p>"
            f"<p>Click <b>Overwrite</b> to replace it, or "
            f"<b>Cancel</b> to keep the existing file and pick a "
            f"different path.</p>"
        )
    else:
        items = "".join(f"<li><code>{p}</code></li>" for p in existing)
        body = (
            f"<p>{len(existing)} {what}s already exist:</p>"
            f"<ul>{items}</ul>"
            f"<p>Click <b>Overwrite</b> to replace them all, or "
            f"<b>Cancel</b>.</p>"
        )
    overwrite_modal_body.text = body
    _overwrite_pending["callback"] = callback
    overwrite_modal_backdrop.visible = True
    overwrite_modal_card.visible = True


def _on_overwrite_yes() -> None:
    cb = _overwrite_pending["callback"]
    _overwrite_pending["callback"] = None
    overwrite_modal_backdrop.visible = False
    overwrite_modal_card.visible = False
    if cb is not None:
        cb()


def _on_overwrite_no() -> None:
    _overwrite_pending["callback"] = None
    overwrite_modal_backdrop.visible = False
    overwrite_modal_card.visible = False
    _set_status("Save cancelled — existing file preserved.", "info",
                clear_after=10)


overwrite_modal_yes_btn.on_click(_on_overwrite_yes)
overwrite_modal_no_btn.on_click(_on_overwrite_no)
overwrite_modal_top_close_btn.on_click(_on_overwrite_no)


def _open_stats_bar_modal(_e=None) -> None:
    stats_bar_modal_backdrop.visible = True
    stats_bar_modal_card.visible = True


def _close_stats_bar_modal(_e=None) -> None:
    stats_bar_modal_backdrop.visible = False
    stats_bar_modal_card.visible = False


def _open_catalog_hover_modal(_e=None) -> None:
    catalog_hover_modal_backdrop.visible = True
    catalog_hover_modal_card.visible = True


def _close_catalog_hover_modal(_e=None) -> None:
    catalog_hover_modal_backdrop.visible = False
    catalog_hover_modal_card.visible = False


stats_bar_open_btn.on_click(_open_stats_bar_modal)
stats_bar_modal_close_btn.on_click(_close_stats_bar_modal)
stats_bar_modal_top_close_btn.on_click(_close_stats_bar_modal)
catalog_hover_open_btn.on_click(_open_catalog_hover_modal)
catalog_hover_modal_close_btn.on_click(_close_catalog_hover_modal)
catalog_hover_modal_top_close_btn.on_click(_close_catalog_hover_modal)

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

# SIDEBAR_W / HELPPANEL_W are defined above (near the figure init) and
# referenced by the figure-sizing comment. Don't redeclare here.

# Compact sidebar organised as tabs around the natural workflow:
# Load (image + catalog) → Aim (pointing + PA + visibility) →
# Pick (instrument, layers, slitlet, filters, undo/clear) →
# Save (session save/load, APT export).

# A little left/right breathing room inside every tab, created within the
# tab's own column so the Tabs/sidebar widths are unchanged.
_TAB_PAD = {"padding": "2px 4px 0 8px"}

image_tab = TabPanel(title="Input", child=column(
    Div(text="<div style='font-size:12px; color:#5a6b85; margin:2px 0 6px'>"
             "Load a background image and target catalogs.</div>",
        width=SIDEBAR_W - 20),
    load_image_open_btn,
    _tab_caption("An example field, a local FITS, or a JPG/PNG + WCS sidecar."),
    load_catalog_open_btn,
    _tab_caption("CSV / ASCII / FITS with ID, RA, Dec — add several to layer."),
    _section_header("Loaded catalogs",
                    "Target catalogs layered on the image. Tick to show/hide, "
                    "▲▼ to reorder draw order, ✕ to remove."),
    catalog_list_column,
    catalog_edit_btn,
    _section_header("Display filters",
                    "Show only catalog sources at or below these limits; "
                    "leave blank to show all."),
    catalog_priority_input,
    catalog_mag_input,
    width=SIDEBAR_W - 6, styles=_TAB_PAD,
))

aim_tab = TabPanel(title="Pointing", child=column(
    _section_header("Disperser / filter",
                    "NIRSpec disperser + filter — sets the wavelength range and "
                    "the dispersed-spectrum length used for collision checks."),
    disperser_filter_select,
    _section_header("Pointing center",
                    "Sky RA/Dec the MSA aperture is centred on. Shift-click the "
                    "image to move it."),
    row(ra_input, dec_input),
    _section_header("Rotation",
                    "Aperture roll on sky. Type the V3 PA or NIRSpec APA "
                    "directly (APA = V3 PA + 138.575°)."),
    v3pa_slider,
    row(v3pa_input, apa_input),
    pa_help_div,
    _section_header("Visibility window",
                    "Compute the V3 PA range JWST allows for a given date "
                    "(from the JWST visibility tool)."),
    row(visibility_date_input, visibility_btn),
    visibility_div,
    _section_header("Optimize MSA pointing",
                    "Search (RA, Dec, V3 PA) for the placement that captures the "
                    "most targets in operable, well-centred shutters — by count, "
                    "weight, or priority tier. Opens the optimizer dialog."),
    opt_open_btn,
    _section_header("MPT configurations",
                    "Plan up to 5 separate exposures (APT-style); choose how many "
                    "and which one you're editing. 'View MPT catalog…' lists the "
                    "sources placed across all of them."),
    row(mpt_num_configs_spinner, mpt_config_select, spacing=10),
    mpt_active_config_div,
    mpt_view_btn,
    width=SIDEBAR_W - 6, styles=_TAB_PAD,
))

pick_tab = TabPanel(title="Settings", child=column(
    _section_header("Layers", "Show or hide overlay layers on the canvas."),
    layers_box,
    _section_header("Slitlet",
                    "Number of shutters opened per click, and whether clicks "
                    "snap to the nearest operable shutter."),
    slitlet_select,
    snap_box,
    _section_header("Overlay appearance",
                    "Per-layer transparency (alpha) and outline width."),
    overlay_layer_select,
    overlay_alpha_slider,
    overlay_stroke_slider,
    _section_header("Canvas (pixels)", "On-screen size of the image canvas."),
    row(canvas_x_spinner, canvas_y_spinner, spacing=10),
    _section_header("Customise display",
                    "Choose which fields appear in the top status bar and the "
                    "catalog hover tooltip."),
    stats_bar_open_btn,
    catalog_hover_open_btn,
    _section_header("Actions",
                    "Undo the last pick, clear all open shutters, or reset "
                    "display settings to defaults."),
    row(undo_btn, clear_btn),
    reset_prefs_btn,
    width=SIDEBAR_W - 6, styles=_TAB_PAD,
))

# MPT tab — everything to do with APT / MPT plans: import (from JSON,
# shutter CSV, or .aptx archive / program ID), the session save/load
# round-trip for collaboration, and the eMPT export bundle.
def _mpt_tab_caption(text):
    return Div(text=f"<div style='font-size:11px; color:#7a8699; "
                    f"margin:0 0 12px 2px'>{text}</div>", width=SIDEBAR_W - 20)


mpt_tab = TabPanel(title="MPT", child=column(
    Div(text="<div style='font-size:12px; color:#5a6b85; margin:2px 0 6px'>"
             "Move plans in and out of vMPT — each opens a dialog.</div>",
        width=SIDEBAR_W - 20),
    mpt_open_import_btn,
    _mpt_tab_caption("An APT/MPT plan, shutter mask, APT program, "
                     "or a saved vMPT session."),
    mpt_open_save_btn,
    _mpt_tab_caption("A shareable vMPT session bundle (your full picking state)."),
    mpt_open_export_btn,
    _mpt_tab_caption(f"The eMPT bundle + {MPT_PLAN_FILENAME} + an "
                     "APT-importable .cat."),
    width=SIDEBAR_W - 6, styles=_TAB_PAD,
))

# Tab strip styling. Document-level CSS can't reach the Tabs' shadow
# root (see the _MODAL_HEADER_STYLES note), so the tab buttons are
# styled via a per-widget InlineStyleSheet. Each tab gets a visible
# border + light fill so it reads as a clickable button even at rest;
# the active tab is white with a blue top accent, and hovering an
# inactive tab highlights it.
_SIDEBAR_TABS_CSS = """
.bk-header { border-bottom: 1px solid #c2d2e6; padding-bottom: 0; }
.bk-tab {
  padding: 4px 13px;
  margin-right: 3px;
  border: 1px solid #cdd8e8;
  border-top: 2px solid transparent;
  border-radius: 6px 6px 0 0;
  background: #f6f8fc;
  color: #41557a;
  font-weight: 600;
  cursor: pointer;
  transition: background .12s, color .12s, border-color .12s;
}
.bk-tab:hover { background: #e6eef9; color: #163a63; border-color: #9db8da; }
.bk-tab.bk-active {
  background: #ffffff;
  border-color: #9db8da;
  border-top: 2px solid #2f6fb3;
  color: #163a63;
}
"""

sidebar_tabs = Tabs(
    tabs=[image_tab, aim_tab, pick_tab, mpt_tab],
    width=SIDEBAR_W,
    stylesheets=[InlineStyleSheet(css=_SIDEBAR_TABS_CSS)],
)

# Sidebar (no stats — stats now live above the figure as a wide bar).
# The MPT tab content (import + save/load + export) is taller than the
# viewport on most laptops, so the sidebar scrolls internally on overflow.
sidebar = column(
    loading_banner,
    sidebar_tabs,
    width=SIDEBAR_W,
    height_policy="max",
    styles={
        "overflow-y": "auto",
        # Leave room at the bottom for the position-fixed status bar
        # (42 px tall) so the last widget in any tab isn't covered.
        "max-height": "calc(100vh - 46px)",
        "padding-bottom": "8px",
        # A shallow tint — distinct from the white central canvas but light
        # enough that the (darkened) section titles read clearly; a crisp
        # divider keeps the panel edge defined.
        "background": "#eef1f7",
        "border-right": "1px solid #d3dae6",
    },
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
    # `stretch_width` so the column itself fills horizontal space in
    # the root row (the figure inside is fixed-pixel; the column lets
    # the stats_div span the full available width).
    sizing_mode="stretch_width",
)

curdoc().add_root(row(
    sidebar, figure_column, help_panel,
    sizing_mode="stretch_both",
))
# Optimizer modal — added as separate roots so the position:fixed
# CSS isn't clipped by the parent row's flex layout.
curdoc().add_root(opt_modal_backdrop)
curdoc().add_root(opt_modal_card)
curdoc().add_root(opt_advanced_modal_backdrop)
curdoc().add_root(opt_advanced_modal_card)
curdoc().add_root(opt_config_modal_backdrop)
curdoc().add_root(opt_config_modal_card)
curdoc().add_root(mpt_view_modal_backdrop)
curdoc().add_root(mpt_view_modal_card)
curdoc().add_root(cat_edit_modal_backdrop)
curdoc().add_root(cat_edit_modal_card)
curdoc().add_root(cat_constraints_modal_backdrop)
curdoc().add_root(cat_constraints_modal_card)
curdoc().add_root(stats_bar_modal_backdrop)
curdoc().add_root(stats_bar_modal_card)
curdoc().add_root(catalog_hover_modal_backdrop)
curdoc().add_root(catalog_hover_modal_card)
# Overwrite-confirmation modal — gates Save-as-CSV and Save-session
# from clobbering existing files. Lives at z-index 1020/1021 so it
# stacks above every other modal.
curdoc().add_root(overwrite_modal_backdrop)
curdoc().add_root(overwrite_modal_card)
# MPT-tab dialogs (Import / Save / Export).
curdoc().add_root(import_modal_backdrop)
curdoc().add_root(import_modal_card)
curdoc().add_root(save_modal_backdrop)
curdoc().add_root(save_modal_card)
curdoc().add_root(export_modal_backdrop)
curdoc().add_root(export_modal_card)
# Input-tab dialogs (Load image / Load catalog).
curdoc().add_root(load_image_modal_backdrop)
curdoc().add_root(load_image_modal_card)
curdoc().add_root(load_catalog_modal_backdrop)
curdoc().add_root(load_catalog_modal_card)
# Status bar — separate root so its position:fixed style escapes the
# sidebar's scrollable container. Lives at the bottom-left of the
# viewport, under the sidebar.
curdoc().add_root(status)
curdoc().title = "vMPT — visual MSA Planning Tool"


# ── User preferences (v1.3.4+) ──────────────────────────────────────────
# Settings the user touches in the Settings tab persist across
# sessions via ~/.vmpt/preferences.json. Loaded once at startup and
# auto-saved on every widget change. See `vmpt/preferences.py` for
# the IO layer.

from .preferences import (
    load_preferences as _load_prefs,
    save_preferences as _save_prefs,
    reset_preferences as _reset_prefs,
)

# When True, prefs-mutating callbacks return early — used while we
# push saved prefs into widgets at startup (otherwise each setattr
# would trigger a write of the partially-applied snapshot).
_prefs_save_suppress = {"flag": False}


def _collect_prefs() -> dict:
    """Snapshot the current widget state into a JSON-serialisable dict."""
    overlay_alphas: dict = {}
    overlay_strokes: dict = {}
    for name, cfg in _OVERLAY_LAYER_CONFIG.items():
        glyph = cfg["glyph"].glyph
        # Field-referenced alpha (the MPT-style overlap layers) lives
        # on `state` — `alpha_attr` is None there. Use the state key
        # instead of trying to read the glyph attribute (which would
        # be the field-name string, not a number).
        if cfg.get("alpha_attr") is None:
            key = cfg.get("alpha_state_key")
            if key and key in state:
                try:
                    overlay_alphas[name] = float(state[key])
                except (TypeError, ValueError):
                    pass
        else:
            try:
                overlay_alphas[name] = float(
                    getattr(glyph, cfg["alpha_attr"])
                )
            except (AttributeError, TypeError, ValueError):
                pass
        try:
            overlay_strokes[name] = float(
                getattr(glyph, cfg["stroke_attr"])
            )
        except (AttributeError, TypeError, ValueError):
            pass
    return {
        "frame_x": int(canvas_x_spinner.value or 800),
        "frame_y": int(canvas_y_spinner.value or 600),
        "slitlet_size": int(slitlet_select.value or 3),
        "snap_to_operable": 0 in (snap_box.active or []),
        "layers_visible": list(layers_box.active or []),
        "overlay_alphas": overlay_alphas,
        "overlay_strokes": overlay_strokes,
        "stats_bar_order": list(stats_bar_choice.value or []),
        "catalog_hover_order": list(catalog_hover_choice.value or []),
        "help_visible": bool(help_div.visible),
        "default_num_configs": int(mpt_num_configs_spinner.value or 1),
    }


def _apply_prefs(prefs: dict) -> None:
    """Push a saved prefs dict back into the widgets. Guarded by
    `_prefs_save_suppress` so the cascade of setattr's doesn't
    trigger another save round-trip per widget change.
    """
    if not prefs:
        return
    _prefs_save_suppress["flag"] = True
    try:
        if "frame_x" in prefs:
            try:
                v = int(prefs["frame_x"])
                canvas_x_spinner.value = v
                state["frame_x"] = v
            except (ValueError, TypeError):
                pass
        if "frame_y" in prefs:
            try:
                v = int(prefs["frame_y"])
                canvas_y_spinner.value = v
                state["frame_y"] = v
            except (ValueError, TypeError):
                pass
        if "slitlet_size" in prefs:
            try:
                slitlet_select.value = str(int(prefs["slitlet_size"]))
            except (ValueError, TypeError):
                pass
        if "snap_to_operable" in prefs:
            snap_box.active = [0] if bool(prefs["snap_to_operable"]) else []
        if "layers_visible" in prefs:
            try:
                layers_box.active = [int(i) for i in prefs["layers_visible"]]
            except (ValueError, TypeError):
                pass
        # Per-layer overlay properties — set glyph attrs directly so
        # `_on_overlay_layer` then re-syncs the sliders' visible
        # values to whichever layer is currently selected.
        # v1.3.0→v1.3.1 migration: the old "Overlapping shutters"
        # key (one orange layer with scalar fill_alpha) became three
        # field-referenced layers driven by state. Map the legacy
        # key onto all three new state buckets so existing prefs
        # files don't lose their setting on first launch.
        prefs_alphas = dict(prefs.get("overlay_alphas") or {})
        prefs_strokes = dict(prefs.get("overlay_strokes") or {})

        def _migrate_overlay_keys(d: dict) -> dict:
            # v1.3.0→v1.3.1: the single "Overlapping shutters" layer became
            # three field-referenced spec-overlap layers — fan its value
            # out so old prefs don't lose the setting.
            if "Overlapping shutters" in d:
                legacy = d.pop("Overlapping shutters")
                for nm in ("Mask Stuck (pink)", "Masked (overlapping warning)",
                           "Mask Conflict (purple)"):
                    d.setdefault(nm, legacy)
            # v1.6.x: "Masked (orange)" → "Masked (overlapping warning)".
            if "Masked (orange)" in d:
                d.setdefault("Masked (overlapping warning)",
                             d.pop("Masked (orange)"))
            return d

        _migrate_overlay_keys(prefs_alphas)
        _migrate_overlay_keys(prefs_strokes)
        for name, alpha in prefs_alphas.items():
            cfg = _OVERLAY_LAYER_CONFIG.get(name)
            if cfg is None:
                continue
            try:
                if cfg.get("alpha_attr") is not None:
                    setattr(cfg["glyph"].glyph, cfg["alpha_attr"],
                            float(alpha))
                elif cfg.get("alpha_state_key"):
                    state[cfg["alpha_state_key"]] = float(alpha)
            except (AttributeError, ValueError, TypeError):
                pass
        for name, stroke in prefs_strokes.items():
            cfg = _OVERLAY_LAYER_CONFIG.get(name)
            if cfg is None:
                continue
            try:
                setattr(cfg["glyph"].glyph, cfg["stroke_attr"],
                        float(stroke))
                extra = cfg.get("stroke_extra")
                if extra is not None:
                    extra(float(stroke))
            except (AttributeError, ValueError, TypeError):
                pass
        # Resync the visible sliders to the currently-selected layer.
        try:
            _on_overlay_layer("value", None, overlay_layer_select.value)
        except Exception:  # noqa: BLE001
            pass
        if "stats_bar_order" in prefs:
            try:
                stats_bar_choice.value = list(prefs["stats_bar_order"])
                state["stats_bar_order"] = [
                    _STATS_BAR_LABEL_TO_KEY[l]
                    for l in stats_bar_choice.value
                    if l in _STATS_BAR_LABEL_TO_KEY
                ]
            except (KeyError, TypeError):
                pass
        if "catalog_hover_order" in prefs:
            try:
                catalog_hover_choice.value = list(
                    prefs["catalog_hover_order"]
                )
                state["catalog_hover_order"] = [
                    _CATALOG_HOVER_LABEL_TO_KEY[l]
                    for l in catalog_hover_choice.value
                    if l in _CATALOG_HOVER_LABEL_TO_KEY
                ]
                try:
                    _refresh_catalog_hover_tooltip()
                except Exception:  # noqa: BLE001
                    pass
            except (KeyError, TypeError):
                pass
        if "help_visible" in prefs:
            try:
                v = bool(prefs["help_visible"])
                help_div.visible = v
                tip_div.visible = v
                help_toggle_btn.label = "Hide help" if v else "Show help"
                help_panel.width = HELPPANEL_W if v else 130
            except Exception:  # noqa: BLE001
                pass
        if "default_num_configs" in prefs:
            try:
                v = max(1, min(_MAX_CONFIGS, int(prefs["default_num_configs"])))
                mpt_num_configs_spinner.value = v
                _ensure_n_configs(v)
            except (ValueError, TypeError):
                pass
    finally:
        _prefs_save_suppress["flag"] = False


def _save_current_prefs_now() -> None:
    """Snapshot + save the full widget state. Suppressed during
    `_apply_prefs` so startup doesn't write back the same prefs
    once per widget."""
    if _prefs_save_suppress["flag"]:
        return
    try:
        _save_prefs(_collect_prefs())
    except Exception:  # noqa: BLE001
        # Disk problems shouldn't break the app — just lose this save.
        pass


def _save_current_prefs(attr, old, new) -> None:
    """Bokeh on_change wrapper — strict (attr, old, new) signature
    (Bokeh validates parameter names exactly, no defaults allowed).
    Args are discarded; we always re-snapshot the full state."""
    _save_current_prefs_now()


# Apply persisted prefs (if any) before wiring auto-save callbacks.
# Doing it in this order means the apply itself doesn't generate
# spurious save events (the suppress flag is on during apply, and
# the change-callbacks haven't been registered yet anyway).
_apply_prefs(_load_prefs())

# Auto-save: every widget change writes the FULL snapshot back to
# disk. Cheap (a few hundred bytes of JSON) and immune to partial
# state — the file is always whatever the user currently sees.
canvas_x_spinner.on_change("value", _save_current_prefs)
canvas_y_spinner.on_change("value", _save_current_prefs)
slitlet_select.on_change("value", _save_current_prefs)
snap_box.on_change("active", _save_current_prefs)
layers_box.on_change("active", _save_current_prefs)
overlay_alpha_slider.on_change("value", _save_current_prefs)
overlay_stroke_slider.on_change("value", _save_current_prefs)
stats_bar_choice.on_change("value", _save_current_prefs)
catalog_hover_choice.on_change("value", _save_current_prefs)
mpt_num_configs_spinner.on_change("value", _save_current_prefs)
# The help-panel toggle button already has `on_help_toggle` bound;
# adding `_save_current_prefs` as a SECOND on_click means both fire
# on every click — the first toggles the panel state, the second
# saves the new state. This wiring is also post-attach, so re-subscribe
# (a no-op today since `on_help_toggle` was bound pre-attach, but it
# keeps the handler alive if that early binding ever moves).
help_toggle_btn.on_click(_save_current_prefs_now)
_resubscribe_late_event_handlers(help_toggle_btn)


def _on_reset_prefs():
    """User clicked 'Reset to defaults'. Wipe the file, then push
    hard-coded defaults into every widget the prefs system covers."""
    _reset_prefs()
    _prefs_save_suppress["flag"] = True
    try:
        canvas_x_spinner.value = 800
        canvas_y_spinner.value = 600
        state["frame_x"] = 800
        state["frame_y"] = 600
        slitlet_select.value = "3"
        snap_box.active = [0]
        layers_box.active = [0, 1, 2]
        # Overlay defaults — restore each layer's OWN default alpha /
        # stroke (`default_alpha` / `default_stroke` in the config), not
        # a uniform value. A blanket 0.20/1.0 would mangle most layers
        # (e.g. shrink the catalog marker size to 1 px, dim stuck-open's
        # outline). Masked / Mask-Stuck / Mask-Conflict default to alpha
        # 0.20, stroke 0.5; the silver operable edge to 0.20 / 1.0.
        for name, cfg in _OVERLAY_LAYER_CONFIG.items():
            try:
                da = cfg.get("default_alpha")
                ds = cfg.get("default_stroke")
                if cfg.get("alpha_attr") is not None:
                    setattr(cfg["glyph"].glyph, cfg["alpha_attr"], da)
                elif cfg.get("alpha_state_key"):
                    state[cfg["alpha_state_key"]] = float(da)
                if ds is not None:
                    setattr(cfg["glyph"].glyph, cfg["stroke_attr"], ds)
                    extra = cfg.get("stroke_extra")
                    if extra is not None:
                        extra(ds)
            except (AttributeError, ValueError, TypeError):
                pass
        try:
            _on_overlay_layer("value", None, overlay_layer_select.value)
        except Exception:  # noqa: BLE001
            pass
        # Customise pickers — full default order.
        stats_bar_choice.value = [
            STATS_BAR_CELL_LABELS[k] for k in STATS_BAR_DEFAULT_ORDER
        ]
        state["stats_bar_order"] = list(STATS_BAR_DEFAULT_ORDER)
        catalog_hover_choice.value = [
            CATALOG_HOVER_FIELDS[k][0] for k in CATALOG_HOVER_DEFAULT_ORDER
        ]
        state["catalog_hover_order"] = list(CATALOG_HOVER_DEFAULT_ORDER)
        try:
            _refresh_catalog_hover_tooltip()
        except Exception:  # noqa: BLE001
            pass
        # Help panel back on (matches the v1.3.3 default-on change).
        help_div.visible = True
        tip_div.visible = True
        help_toggle_btn.label = "Hide help"
        help_panel.width = HELPPANEL_W
    finally:
        _prefs_save_suppress["flag"] = False
    # Apply the new sizes to the canvas if an image is loaded.
    if state.get("image") is not None:
        refresh_image_glyph()
    refresh_overlays()
    _set_status("Display settings reset to defaults.", "ok",
                clear_after=10)


# `reset_prefs_btn` was forward-declared near the other action
# buttons so the Settings tab layout can include it. Wire the
# handler now that `_on_reset_prefs` exists.
reset_prefs_btn.on_click(_on_reset_prefs)
# Bokeh gotcha: `on_event`/`on_click` only enrols a model in the
# document's event dispatcher (`Document.callbacks._subscribed_models`)
# at *attach* time, via `_update_event_callbacks()` in
# `_attach_document`. This button is forward-declared and only wired
# here — long after its Settings-tab layout root was added to
# `curdoc()` above — so the subscription never happened and the server
# silently dropped its ButtonClick events (the "Reset display to
# defaults does nothing" bug). Re-run the subscription explicitly now
# that the handler is attached. Idempotent (it `set.add`s the model).
_resubscribe_late_event_handlers(reset_prefs_btn)


# ── CLI auto-load ─────────────────────────────────────────────────────
# Args forwarded by `run.sh --args ...` arrive in sys.argv. We pre-fill
# the relevant path inputs; the existing on_change handlers do the
# actual loading. Loads are deferred to the next IO tick so the
# document is fully wired before we trigger heavy work.
def _parse_startup_args(argv: list[str]) -> dict:
    """Tolerant parser for `--fits`, `--jpg`, `--wcs`, `--catalog`
    (repeatable), and a default roll via `--v3pa` or `--apa` (degrees).

    Unknown args are silently ignored — Bokeh prefixes some of its own
    flags before --args and we don't want to throw on them."""
    out: dict = {"fits": None, "jpg": None, "wcs": None, "catalogs": [],
                 "v3pa": None, "apa": None}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--fits", "--jpg", "--wcs", "--v3pa", "--apa") and i + 1 < len(argv):
            out[a[2:]] = argv[i + 1]
            i += 2
        elif a == "--catalog" and i + 1 < len(argv):
            out["catalogs"].append(argv[i + 1])
            i += 2
        else:
            i += 1
    return out


def _autoload_from_args() -> None:
    """Sequence the path-loads the user asked for via run.sh.

    Loads must run one at a time: a catalog overlay needs the image's
    WCS to project source positions to pixels. If they fire in the
    same tick the canvas aspect ratio comes out wrong because the
    catalog refresh races the image-set-and-recenter handler. Each
    loader now accepts an `on_complete` callback we chain through.
    """
    args = _parse_startup_args(list(sys.argv[1:]))

    # Build an explicit step queue. Each step receives the next-step
    # callback so it can chain after its load finishes.
    steps: list = []

    if args["fits"]:
        fits_path = args["fits"]
        def step_fits(cb, p=fits_path):
            if not Path(p).exists():
                _set_status(f"--fits path not found: {p}", "err")
                cb()
                return
            # Set the input AND mirror the value through the loader
            # directly with `on_complete` so we can chain. Setting the
            # TextInput also fires `on_fits_path` which kicks off its
            # own _deferred load; that's redundant but harmless (the
            # second load is a no-op since state["image"] is already
            # set, but we suppress that race by NOT calling on_fits_path
            # — set the value silently is fine because we own this codepath).
            fits_path_input.value = p
            _show_loading(f"Loading FITS: {Path(p).name}…")
            _deferred(_load_fits_from_path, p, on_complete=cb)
        steps.append(step_fits)
    elif args["jpg"] and args["wcs"]:
        jpg_path = args["jpg"]
        wcs_path = args["wcs"]
        def step_jpg(cb, j=jpg_path, w=wcs_path):
            if not Path(j).exists() or not Path(w).exists():
                missing = j if not Path(j).exists() else w
                _set_status(f"--jpg/--wcs path not found: {missing}", "err")
                cb()
                return
            # Surface the paths in the inputs (mostly cosmetic; saves
            # the user from re-typing if they want to reload).
            sidecar_path_input.value = w
            jpg_path_input.value = j
            _show_loading(f"Loading JPG + WCS sidecar…")
            _deferred(_load_jpg_pair_from_paths, j, w, on_complete=cb)
        steps.append(step_jpg)

    for cat_path in args["catalogs"]:
        def step_cat(cb, c=cat_path):
            if not Path(c).exists():
                _set_status(f"--catalog path not found: {c}", "err")
                cb()
                return
            # Surface the FIRST catalog's path in the input box too —
            # later ones don't, since the input only holds one path.
            if not catalog_path_input.value:
                catalog_path_input.value = c
            _show_loading(f"Loading catalog: {Path(c).name}…")
            _deferred(_load_catalog_from_path, c, on_complete=cb)
        steps.append(step_cat)

    # Default pointing roll from --v3pa / --apa. Applied LAST (after the
    # image + catalogs load) so the MSA overlay renders at the requested
    # angle. --v3pa wins if both are given. APA = V3 PA + V3IdlYAngle.
    pa_target = None
    try:
        if args.get("v3pa") is not None:
            pa_target = float(args["v3pa"]) % 360.0
        elif args.get("apa") is not None:
            pa_target = (float(args["apa"]) - V3_IDL_Y_ANGLE) % 360.0
    except (TypeError, ValueError):
        _set_status("--v3pa / --apa must be a number in degrees; ignoring.",
                    "err")

    if pa_target is not None:
        def step_pa(cb, v=pa_target):
            _sync_pa_widgets(v)
            # Re-render the overlay at the new roll only if an image is
            # loaded; otherwise just leaving the PA widgets set is enough.
            if state.get("image") is not None:
                _rebuild_shutter_catalog_index()
                refresh_overlays()
            cb()
        steps.append(step_pa)

    # Chain runner — each step calls its `cb` (= run_next) in `finally`
    # so the next step starts only after the previous one releases.
    # `state["_autoload_active"]` muzzles the path-input on_change
    # handlers while we're driving — otherwise a value-set would
    # trigger its own _deferred load and race with our sequenced one.
    def run_next():
        if not steps:
            state["_autoload_active"] = False
            _set_status("Autoload complete.", "ok", clear_after=4)
            return
        steps.pop(0)(run_next)

    if steps:
        state["_autoload_active"] = True
        run_next()


curdoc().add_next_tick_callback(_autoload_from_args)


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


# ── PyPI version-check (for `pip install jwst-vmpt` users) ─────────────────
# A git checkout is nagged against GitHub `main` (above); an installed
# package is nagged against the latest release on PyPI instead.
_PYPI_JSON_URL = "https://pypi.org/pypi/jwst-vmpt/json"


def _installed_vmpt_version() -> str | None:
    """Installed jwst-vmpt version, or None if it can't be determined
    (e.g. running from a source tree that was never pip-installed)."""
    try:
        from importlib.metadata import PackageNotFoundError, version
        return version("jwst-vmpt")
    except Exception:  # noqa: BLE001 — PackageNotFoundError or anything else
        return None


def _pypi_is_newer(installed: str, latest: str) -> bool:
    """True iff `latest` (from PyPI) is a newer release than `installed`.
    Uses PEP 440 ordering (so 1.10.0 > 1.9.0); falls back to a string
    compare if `packaging` is unavailable or a version is unparseable."""
    if not installed or not latest:
        return False
    try:
        from packaging.version import parse as _parse_version
        return _parse_version(latest) > _parse_version(installed)
    except Exception:  # noqa: BLE001
        return latest != installed


def _on_pypi_update_available(installed: str, latest: str) -> None:
    """Runs on the Bokeh document thread. Reveals the PyPI update notice
    (reuses the floating update banner)."""
    inst = (installed or "").replace("<", "&lt;").replace(">", "&gt;")
    new = (latest or "").replace("<", "&lt;").replace(">", "&gt;")
    update_div.text = (
        '<div style="display:flex; align-items:flex-start; gap:10px;">'
        '  <div style="font-size:22px; line-height:1; flex-shrink:0;">📦</div>'
        '  <div style="min-width:0;">'
        '    <div style="font-weight:700; color:#7a5d00; font-size:13px; '
        '                margin-bottom:4px; letter-spacing:0.2px;">'
        '      vMPT update available'
        '    </div>'
        '    <div style="font-size:11.5px; color:#6f5b1f; margin-bottom:8px;">'
        f'      You have <code>{inst}</code>; the latest on PyPI is '
        f'      <code>{new}</code>.'
        '    </div>'
        '    <div style="font-size:11.5px; color:#6f5b1f; margin-bottom:5px;">'
        '      Update with:'
        '    </div>'
        '    <div style="font-family:monospace; font-size:11.5px; '
        '                background:#fff7db; border:1px solid #e0c873; '
        '                border-radius:4px; padding:4px 7px; color:#5a4500;">'
        '      pip install -U jwst-vmpt'
        '    </div>'
        '    <div style="font-size:10.5px; color:#8a7530; margin-top:6px;">'
        '      Then restart vMPT.'
        '    </div>'
        '  </div>'
        '</div>'
    )
    update_box.visible = True


def _check_pypi_updates_blocking(session_doc) -> None:
    """Background thread: compare the installed version against the latest
    PyPI release and nag (once) if PyPI is newer. Offline / parse / not-
    installed failures are silent — it's a courtesy, not a dependency."""
    installed = _installed_vmpt_version()
    if not installed:
        return
    try:
        req = _urllib_request.Request(
            _PYPI_JSON_URL,
            headers={"Accept": "application/json",
                     "User-Agent": "vMPT-version-check"},
        )
        with _urllib_request.urlopen(req, timeout=6) as resp:
            data = _json.loads(resp.read().decode("utf-8"))
    except (_urllib_error.URLError, _urllib_error.HTTPError, TimeoutError,
            OSError, _json.JSONDecodeError, UnicodeDecodeError):
        return  # offline / rate-limited / unparseable — silently skip
    latest = ((data.get("info") or {}).get("version") or "").strip()
    if not _pypi_is_newer(installed, latest):
        return
    try:
        session_doc.add_next_tick_callback(
            lambda: _on_pypi_update_available(installed, latest)
        )
    except Exception:  # noqa: BLE001
        pass


def _check_for_updates_blocking(session_doc) -> None:
    """Run in a background thread. `session_doc` is the user-session doc
    captured BEFORE the thread started — curdoc() is thread-local, so
    re-querying it from this thread returns a different (irrelevant)
    document and the next-tick callback never fires on the user's UI."""
    local = _local_git_head()
    if not local:
        # Not a git checkout (pip install) → compare against PyPI instead.
        _check_pypi_updates_blocking(session_doc)
        return
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
if (_session_doc.session_context is not None
        and os.environ.get("VMPT_UPDATE_CHECK", "1") != "0"):
    _threading.Thread(
        target=_check_for_updates_blocking,
        args=(_session_doc,),
        daemon=True,
    ).start()
