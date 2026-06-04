"""vMPT — visual MSA Planning Tool. Bokeh server entry point.

Run:  bokeh serve vmpt/ --show   (or `vmpt` after `pip install jwst-vmpt`)
"""
from __future__ import annotations

import base64
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

from bokeh.events import DoubleTap, RangesUpdate, Tap
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
    HTMLTemplateFormatter,
    MultiChoice,
    NumberEditor,
    NumberFormatter,
    Range1d,
    RadioGroup,
    Select,
    Slider,
    StringEditor,
    TableColumn,
    TabPanel,
    Tabs,
    TextInput,
    Toggle,
    WheelZoomTool,
)
from bokeh.plotting import figure

from vmpt.catalog import Catalog, catalog_in_view, load_catalog
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
from vmpt.msa import load_msa_grid, load_operability
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
    export_session_json,
    import_session_json,
)
from vmpt.wavelengths import (
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
    "catalog": None,         # merged-active cache; rebuilt from "catalogs"
    "catalogs": [],          # list of {"name", "catalog", "enabled"} entries
                             # — the source of truth. `state["catalog"]` is
                             # recomputed on every add/remove/toggle.
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
SIDEBAR_W = 340      # left tab panel (Input / Pointing / Setting / MPT)
HELPPANEL_W = 340    # right help panel (Quick guide + rotating tip)

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


def _wrap_path_picker(text_input, browse_btn, *,
                      header_html: str = "", edit_label: str = "Edit path"):
    """Wrap a (path TextInput + Browse button) pair so Browse is the
    primary affordance and the path TextInput is hidden behind an
    "Edit path" toggle.

    The TextInput stays hidden when empty but auto-reveals as soon as
    it has a value (so paths populated by Browse, by autoload, or by
    direct typing all surface to the user). Clicking the toggle
    re-hides it. This keeps the Input + MPT tabs tidy by default
    without forcing power-users to lose visibility of their loaded
    paths.
    """
    visible_now = bool((text_input.value or "").strip())
    text_input.visible = visible_now
    edit_btn = Button(
        label=("Hide path" if visible_now else edit_label),
        button_type="default", width=80, height=26,
        css_classes=["vmpt-help-toggle"],
    )

    def _toggle(_e=None):
        text_input.visible = not text_input.visible
        edit_btn.label = "Hide path" if text_input.visible else edit_label

    edit_btn.on_click(_toggle)

    def _on_value(attr, old, new):
        # Auto-reveal when a path is filled (Browse pick or autoload),
        # but keep the user's manual toggle state otherwise.
        if (new or "").strip() and not text_input.visible:
            text_input.visible = True
            edit_btn.label = "Hide path"
    text_input.on_change("value", _on_value)

    parts = []
    if header_html:
        parts.append(Div(text=header_html, width=SIDEBAR_W - 20))
    parts.append(row(browse_btn, edit_btn, spacing=6))
    parts.append(text_input)
    return column(*parts, spacing=2)

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
    "<small><b>V3 PA</b>: JWST V3 axis PA on sky (drives the overlay). "
    f"<b>NIRSpec APA</b>: aperture PA of NRS_FULL_MSA = V3PA + {V3_IDL_Y_ANGLE:.2f}° "
    "(mod 360). APT/MPT calls this NIRSpec's 'Aperture PA'. "
    "<a href='https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/"
    "jwst-position-angles-ranges-and-offsets' target='_blank'>JDox reference</a>.</small>"
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
opt_n_top_input = TextInput(title="Refine top N", value="10", width=_HALF_W)
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
    width=SIDEBAR_W - 20,
    height=24,
    css_classes=["vmpt-help-toggle"],
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
opt_priority_input = TextInput(
    title="Priority cutoff ≤ (blank = all)", value="", placeholder="e.g. 1",
    width=SIDEBAR_W - 20,
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

_ADV_INPUT_W = 240
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
    Div(text="<h3 style='margin:0 0 6px 0; color:#1a3b66'>"
             "Advanced optimizer settings</h3>"
             "<div style='font-size:12px; color:#5a6b85'>"
             "Tune only if the defaults don't fit. Values stick after Done.</div>",
        width=520),
    row(opt_grid_n_ra_input, opt_grid_n_dec_input, spacing=12),
    row(opt_grid_n_pa_input, opt_de_maxiter_input, spacing=12),
    opt_objective_select,
    opt_sigma_input,
    opt_theta_input,
    opt_advanced_modal_close_btn,
    spacing=10,
    width=540,
    visible=False,
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
  .vmpt-help-toggle button {
    background: transparent;
    border: 0;
    color: #5a6b85;
    text-align: left;
    padding: 2px 0;
    font-size: 12px;
    cursor: pointer;
  }
  .vmpt-help-toggle button:hover { color: #1a3b66; }
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
    # Inline styles work even when the modal's stylesheet doesn't
    # penetrate Bokeh's nested wrappers; pinning the wrapper div to
    # the corner of the modal card via absolute positioning.
    styles={
        "position": "absolute",
        "top": "6px", "right": "8px",
        "z-index": "5",
    },
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
    # Top-right × dismiss. CSS in `_cat_edit_css` floats it into the
    # corner of the modal card — Bokeh layouts don't support absolute
    # positioning natively, so the button sits at the top of the
    # column and CSS does the rest.
    cat_edit_top_close_btn,
    Div(text="<h3 style='margin:0 0 4px 0; color:#1a3b66'>"
             "Edit catalog</h3>"
             "<div style='font-size:12px; color:#5a6b85'>"
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
    styles={
        "position": "absolute",
        "top": "6px", "right": "8px",
        "z-index": "5",
    },
)
cat_constraints_title_div = Div(
    text="<h3 style='margin:0 0 4px 0; color:#1a3b66'>"
         "Per-target spectral constraints</h3>"
         "<div style='font-size:12px; color:#5a6b85'>"
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
    cat_constraints_top_close_btn,
    cat_constraints_title_div,
    cat_constraints_row_label,
    cat_constraints_lam_input,
    cat_constraints_lam_warn,
    cat_constraints_checks,
    cat_constraints_centration_select,
    cat_constraints_centration_hint,
    row(cat_constraints_apply_btn, cat_constraints_cancel_btn, spacing=10),
    spacing=10,
    width=460,
    visible=False,
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


opt_run_btn = Button(label="Run optimization",
                    button_type="primary", width=SIDEBAR_W - 20)
# The Pointing-tab CTA that opens the optimizer's config modal. The
# config widgets used to live inline in the Pointing tab; on a 913 px
# laptop screen they pushed the actual Run button below the fold.
# Wrapping them in a modal keeps the Pointing tab compact.
opt_open_btn = Button(
    label="Open optimizer…",
    button_type="primary",
    width=SIDEBAR_W - 20,
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
    styles={
        "position": "absolute",
        "top": "6px", "right": "8px",
        "z-index": "5",
    },
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
    opt_config_top_close_btn,
    Div(text="<h3 style='margin:0 0 4px 0; color:#1a3b66'>"
             "MSA pointing optimizer</h3>"
             "<div style='font-size:12px; color:#5a6b85'>"
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
    opt_priority_input,
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
opt_modal_results_summary = Div(text="", width=560)
opt_modal_results_rows = column(spacing=0, width=560)
opt_modal_results_box = column(
    opt_modal_results_summary,
    opt_modal_results_rows,
    spacing=4,
    width=560,
    visible=False,
)

opt_modal_close_btn = Button(label="Close", button_type="default", width=80)

opt_modal_card = column(
    opt_modal_title,
    opt_modal_progress_box,
    opt_modal_results_box,
    opt_modal_close_btn,
    visible=False,
    spacing=10,
    # Wider than before to accommodate the Hierarchy mode's per-tier
    # breakdown in the Score column (e.g. "P0:4·P1:12·P2:30 (46)" —
    # ~200 px) plus the rest of the row (~480 px) and the modal's
    # inner padding (~36 px).
    width=740,
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
    ("✋", "Move the pointing", "<b>Shift + click</b> anywhere on the image to recentre the pointing on that spot. The <span style='color:#2e9b3f;font-weight:600'>lime cross</span> marks the current pointing."),
    ("🔁", "Toggle a slitlet", "Click an already-open shutter to close it. Its slitlet siblings come down with it."),
    ("🎨", "Cyan flag", "Double-click a shutter to toggle a <span style='color:#0aa;font-weight:600'>cyan highlight</span> — a visual flag for your own review. It's not exported."),
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
  <li>Pick the <b>N-shutter slitlet</b> size (1/2/3/5) in <b>Setting</b>.</li>
  <li><b>Click</b> → opens N-shutter slitlet at the nearest operable shutter. Click an open shutter to close the slitlet.</li>
  <li><b>Double-click</b> → toggles <span style='color:#0aa;font-weight:600;background:#222;padding:0 4px'>cyan highlight</span> (visual flag, not exported).</li>
  <li>Layers (Setting tab → <b>Layers</b>):
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
  <li><b>Reset</b>: toolbar · <b>Undo</b>: Setting → <b>Undo last</b></li>
</ul>
</div>
<p style='margin:4px 0'>Full reference in <code>README.md</code> · file roles in <code>CONTEXT.md</code>.</p>
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


# Collapsed by default — the Quick guide + rotating tip are useful
# for first-run users but eat ~340 px of horizontal real estate every
# session. Returning users get a wider workspace; one click on
# "Show help" brings the panel back.
help_div.visible = False
tip_div.visible = False
help_toggle_btn.on_click(on_help_toggle)
help_panel = column(
    help_toggle_btn, tip_div, help_div,
    width=130,
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
    start=0.0, end=3.0, step=0.05, value=1.0,
    title="Stroke (px)",
    width=SIDEBAR_W - 40,
)

# Canvas size — two independent sliders for the X (width) and Y
# (height) of the figure's drawing frame in pixels.
# `refresh_image_glyph()` reads state["frame_x"] / state["frame_y"]
# and sets fig.frame_width / fig.frame_height directly.
# `match_aspect=True` on the figure keeps image pixels square — when
# the canvas aspect ≠ image aspect, Bokeh letterboxes the image so
# both the science image AND the NIRSpec FoV stay at their correct
# aspect ratios at any (frame_x, frame_y). Default 800×800.
canvas_x_slider = Slider(
    start=400, end=1600, step=50, value=800,
    title="Canvas width (X, px)",
    width=SIDEBAR_W - 40,
)
canvas_y_slider = Slider(
    start=400, end=1600, step=50, value=800,
    title="Canvas height (Y, px)",
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
    styles={
        "position": "absolute",
        "top": "6px", "right": "8px",
        "z-index": "5",
    },
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
    styles={
        "position": "absolute",
        "top": "6px", "right": "8px",
        "z-index": "5",
    },
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
# open shutters (user picks), highlighted shutters, MSA outline, fixed slits,
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
# Spectral-overlap shutters: any operable shutter in the same s row of the
# same quadrant as a currently-open shutter. If opened, their dispersed
# spectra would overlap on the detector (MPT-style spectral conflict).
spec_overlap_glyph = fig.multi_polygons(
    xs="xs", ys="ys", source=src_spec_overlap,
    # No edge by default — overlap shutters are fill-only, alpha 0.2
    # per conflict so the colour intensifies where multiple open
    # shutters' spectra contribute (alpha compositing stacks). The
    # stroke slider in the appearance picker reveals the edge by
    # bumping line_alpha; line_color is set to orange here so the
    # outline matches the fill (default would be Bokeh's blue-grey).
    line_color="#d97a00", line_alpha=0.0, line_width=0,
    fill_color="orange", fill_alpha=0.20,
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
    stats_div.text = (
        '<div style="display:flex; flex-wrap:wrap; align-items:center; gap:0;">'
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
    # Re-fit figure ranges to the new image's pixel extent. We use a
    # single atomic update per axis (via `.update()`) so match_aspect's
    # aspect-locking pass sees a coherent change and re-evaluates the
    # constraint. Setting `.start` and `.end` separately on DataRange1d
    # can leave the ranges in an intermediate state where match_aspect
    # doesn't fire — that was the "switch from A370 to a JPG and the
    # aspect comes out wrong" bug.
    fig.x_range.update(start=0, end=W)
    fig.y_range.update(start=0, end=H)
    # Lock the canvas frame pixel dimensions to the image's W:H so 1
    # image pixel = (FRAME_SCALE × W / max(W, H)) screen pixels,
    # consistently in X and Y. Window resizes leave these alone —
    # the canvas is a fixed-pixel block; the layout column letterboxes
    # around it. Trade-off: the canvas doesn't grow on big monitors,
    # but image pixels are guaranteed square.
    if W > 0 and H > 0:
        # User-adjustable canvas dimensions (state-backed so values
        # survive image reloads). Each axis is set directly here;
        # `match_aspect=True` on the figure then enforces 1:1 pixel
        # aspect by adjusting the data range — image pixels stay
        # square and the MSA overlay's data range stays consistent.
        # When the canvas aspect ≠ image aspect, the image is
        # letterboxed inside the frame so both stay at their
        # correct ratios at any (frame_x, frame_y).
        # Default 800×800; the legacy "frame_max" key (pre-split)
        # still drives both axes if present, for session-reload
        # backwards compatibility.
        legacy = state.get("frame_max")
        frame_x = int(state.get("frame_x", legacy or 800))
        frame_y = int(state.get("frame_y", legacy or 800))
        fig.frame_width = max(100, frame_x)
        fig.frame_height = max(100, frame_y)
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
                if 0 <= i < len(state["catalogs"]):
                    state["catalogs"][i]["enabled"] = (0 in new)
                    _rebuild_merged_catalog()
                    _rebuild_shutter_catalog_index()
                    refresh_overlays()
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
    if not Path(p).exists():
        _set_status(f"Catalog path not found: {p}", "err")
        return
    _show_loading(f"Loading catalog: {Path(p).name}…")
    _deferred(_load_catalog_from_path, p)


def on_catalog_add():
    """Explicit Add button — re-uses whatever's in catalog_path_input."""
    p = catalog_path_input.value.strip()
    if not p:
        _set_status("Set a catalog path first.", "warn")
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

    def _cat_lookup(tid: str) -> tuple[float, float, int, str] | None:
        """Look up a catalog row by id. Returns (ra, dec, weight, label).
        `label` is the catalog's `label`/`name` column value when
        available, else the original raw ID (so the trace from MPT-side
        integer ID back to the user's source-list name is preserved)."""
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
            label_val = tid  # original string ID as the Label
        return float(cat.ra_deg[k]), float(cat.dec_deg[k]), pr, label_val

    def _push_target_row(
        raw_tid: str | None,
        ra_d: float,
        dec_d: float,
        pr: int,
        label: str,
    ) -> int:
        """Append a target row. Always returns the assigned integer ID.

        If `raw_tid` is None, a fresh sequential integer is generated.
        Otherwise the integer is derived from `raw_tid` via _to_int_id
        (e.g. "RJ0600-10274-P0" → 10274).
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
            "Pr": pr,
            "ra_deg": ra_d,
            "dec_deg": dec_d,
            "label": label,
        })
        return target_no

    # Step 1: open shutters with a known catalog source. The output
    # `.cat` row's integer ID is derived from the original catalog
    # token (digit-run extraction). The original token survives in the
    # Label column for downstream traceability.
    for (q, s, d), sh in state["open_shutters"].items():
        tid = sh.target_id or _shutter_source_id(q, s, d)
        if not tid or tid in real_ids_seen:
            continue
        info = _cat_lookup(str(tid))
        if info is None:
            continue  # synthesise later from geometry
        ra_d, dec_d, pr, label_val = info
        real_ids_seen.add(str(tid))
        _push_target_row(str(tid), ra_d, dec_d, pr, label=label_val)

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
                pr,
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
                None, ra_d, dec_d, pr=5, label="vMPT_synth",
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
        f"Session loaded: {len(state['open_shutters'])} open shutters, "
        f"{len(state['highlighted'])} highlighted.{image_note}{catalog_note}",
        "warn" if (image_note or catalog_note) else "ok", clear_after=14,
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


def _on_canvas_x(attr, old, new):
    """Resize the figure frame's WIDTH in response to the Settings
    slider. Stashes the value on `state["frame_x"]` so the change
    survives image reloads; calls :func:`refresh_image_glyph` to
    apply immediately when an image is loaded. `match_aspect=True`
    on the figure handles pixel-aspect locking, so a non-square
    canvas letterboxes the image to keep the MSA FoV correct.
    """
    try:
        v = int(new)
    except (TypeError, ValueError):
        return
    v = max(400, min(1600, v))
    state["frame_x"] = v
    if state.get("image") is not None:
        refresh_image_glyph()


def _on_canvas_y(attr, old, new):
    """Same as :func:`_on_canvas_x` for the Y (height) axis."""
    try:
        v = int(new)
    except (TypeError, ValueError):
        return
    v = max(400, min(1600, v))
    state["frame_y"] = v
    if state.get("image") is not None:
        refresh_image_glyph()


canvas_x_slider.on_change("value", _on_canvas_x)
canvas_y_slider.on_change("value", _on_canvas_y)


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
            protect=[], _has_constraint=[],
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
        has_constraint_list = [
            int(bool(lam_req_list[i]) or no_gap_list[i] or
                extend_blue_list[i] or extend_red_list[i]
                or protect_list[i] or bool(centration_list[i]))
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
        "centration",
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
    needs_fix = False
    new_dict = {k: list(v) for k, v in new.items()}
    for col in ("priority", "weight", "z"):
        if col not in new_dict:
            continue
        vals = new_dict[col]
        if not any(isinstance(v, str) for v in vals):
            continue
        coerced = []
        for v in vals:
            if isinstance(v, (int, float)) and not (
                isinstance(v, float) and not np.isfinite(v)
            ):
                coerced.append(float(v))
                continue
            s = ("" if v is None else str(v)).strip()
            if not s or s.lower() in ("nan", "none", "null", "--"):
                coerced.append(float("nan"))
            else:
                try:
                    coerced.append(float(s))
                except ValueError:
                    coerced.append(float("nan"))
        new_dict[col] = coerced
        needs_fix = True
    if needs_fix:
        # Suppress recursion: setting `data` would re-trigger us.
        _cat_edit_set_data_silently(new_dict)


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
    """Write the working copy to CSV at the user-supplied path."""
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
                     "centration"}
        extras_cols = [k for k in data.keys()
                       if k not in _CAT_STD_COLS and k not in _INTERNAL]
        # Are any per-target constraints set anywhere in the catalog?
        # If not, omit the constraint columns so v1.2.x users get the
        # same CSV format they had before.
        has_constraints = any(
            int(v) for v in (data.get("_has_constraint") or [])
        )
        constraint_cols = (["lam_req", "no_gap", "extend_blue",
                            "extend_red", "protect", "centration"]
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
                    elif k == "centration":
                        # String label or empty — write verbatim
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
    has_any = int(
        bool(lam_val) or new_no_gap or new_blue or new_red or new_protect
        or bool(cent_val)
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


def _apply_optimizer_result(
    ra_p: float, dec_p: float, pa_v3: float,
    *, clear_existing: bool = True,
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
    """
    _push_history()  # snapshot for Undo

    if clear_existing:
        state["open_shutters"] = {}

    ra_input.value = f"{ra_p:.6f}"
    dec_input.value = f"{dec_p:.6f}"
    _sync_pa_widgets(float(pa_v3))

    n_targets = 0
    n_opened = 0
    ev = _opt_run.get("evaluator") if _opt_run else None
    ids = _opt_run.get("source_ids", []) if _opt_run else []
    if ev is not None:
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
    The trigger value is `"<ra>,<dec>,<pa>,<stamp>"`."""
    if not new:
        return
    try:
        ra_s, dec_s, pa_s, _stamp = new.split(",", 3)
        ra_p = float(ra_s); dec_p = float(dec_s); pa_v3 = float(pa_s)
    except (ValueError, TypeError):
        opt_apply_trigger.value = ""
        return
    # Reset the trigger so the next click on the same row still fires.
    # (Suppress the recursive on_change by checking `new` at the top.)
    opt_apply_trigger.value = ""
    _apply_optimizer_result(ra_p, dec_p, pa_v3, clear_existing=True)
    _opt_hide_modal()


opt_apply_trigger.on_change("value", _on_opt_apply_trigger)


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


def _opt_render_results_in_modal(
    results: dict, ra_ref: float, dec_ref: float, pa_ref: float,
    n_sources: int,
    *, method: str = "Democracy",
) -> None:
    """Build the post-run results table: one Bokeh row per solution
    so the Apply button sits exactly next to its cells (no Bokeh
    column-spacing drift between buttons and table rows).

    Per-row enrichments built upstream and read here:
      • `total_count[i]` — count of observable sources at this pointing.
      • `sum_weight[i]`  — Σ weight of observable sources (Meritocracy
        headline).
      • `tier_breakdown[i]` — per-priority-tier counts (Hierarchy
        headline).
      • `top_targets[i]`  — top-10 sources at this pointing sorted by
        priority asc / weight desc. Rendered as a hover tooltip on
        the row's Score cell.
    """
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

    # Summary line.
    method_blurb = {
        "Democracy":   "ranked by raw count.",
        "Meritocracy": "ranked by Σ weight; <b>Σw</b> shown, total count in parens.",
        "Hierarchy":   ("ranked lexicographically by priority tier; "
                        "<b>P<sub>i</sub>:n</b> = sources placed at tier i."),
    }.get(method, "")
    protect_blurb = ""
    if protect_enabled:
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
            breakdown_label = " · ".join(
                f"P{int(t)}:{c}" for t, c in breakdowns[i]
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
        btn.js_on_click(CustomJS(
            args=dict(trig=opt_apply_trigger),
            code=f"""
            const ra = {float(ra_i)};
            const dec = {float(dec_i)};
            const pa = {float(pa_i)};
            const rank = {i + 1};
            const msg = "Apply solution #" + rank +
                        " ?\\n\\nThis will CLEAR all previously open " +
                        "shutters and replace them with the optimizer's " +
                        "slitlets.";
            if (!window.confirm(msg)) {{
                return;
            }}
            // Include a stamp so Bokeh sees a new value even for
            // repeat applies (same row clicked twice in a session).
            trig.value = ra.toFixed(8) + "," + dec.toFixed(8) + "," +
                         pa.toFixed(6) + "," + Date.now() + "_" + Math.random();
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

    opt_modal_progress_box.visible = False
    opt_modal_results_box.visible = True
    opt_modal_results_rows.children = [header_row, *data_rows]


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
                f"{refined['score'][0]:.1f} of {_opt_run['n_sources']} sources.",
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
    _cent_full = np.asarray(getattr(cat, "centration", None) or [],
                            dtype=object)
    if _cent_full.size != len(cat.ra_deg):
        _cent_full = np.array([""] * len(cat.ra_deg), dtype=object)
    centration_per_target_arr = (
        _cent_full[keep] if keep is not None else _cent_full
    )

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
    stats_bar_modal_top_close_btn,
    Div(text="<h3 style='margin:0 0 4px 0; color:#1a3b66'>"
             "Customise top stats bar</h3>"
             "<div style='font-size:12px; color:#5a6b85'>"
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
    css_classes=["vmpt-customise-modal"],
    stylesheets=[_CUSTOMISE_MODAL_CSS],
)

catalog_hover_modal_card = column(
    catalog_hover_modal_top_close_btn,
    Div(text="<h3 style='margin:0 0 4px 0; color:#1a3b66'>"
             "Customise catalog hover</h3>"
             "<div style='font-size:12px; color:#5a6b85'>"
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
    css_classes=["vmpt-customise-modal"],
    stylesheets=[_CUSTOMISE_MODAL_CSS],
)


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

image_tab = TabPanel(title="Input", child=column(
    Div(text="<b>Image</b> — try an example:"),
    row(example_a370_btn, example_r0600_btn),
    _wrap_path_picker(
        fits_path_input, fits_browse_btn,
        header_html="<small><b>or</b> a local FITS:</small>",
    ),
    Div(text="<small><b>or</b> JPG + sidecar FITS:</small>"),
    _wrap_path_picker(
        sidecar_path_input, sidecar_browse_btn,
        header_html="<small>Sidecar FITS (WCS)</small>",
    ),
    _wrap_path_picker(
        jpg_path_input, jpg_browse_btn,
        header_html="<small>JPG / PNG</small>",
    ),
    Div(text="<b>Catalogs</b> <small>(CSV / ASCII / FITS with ID, RA, DEC; "
             "you can load multiple — each can be toggled on/off or removed)</small>"),
    _wrap_path_picker(
        catalog_path_input, catalog_browse_btn,
        header_html="",
    ),
    row(catalog_add_btn),
    catalog_list_column,
    catalog_edit_btn,
    catalog_priority_input,
    catalog_mag_input,
    width=SIDEBAR_W - 20,
))

aim_tab = TabPanel(title="Pointing", child=column(
    Div(text="<b>Disperser / Filter</b>"),
    disperser_filter_select,
    Div(text="<b>Pointing center</b>"),
    row(ra_input, dec_input),
    Div(text="<b>Rotation</b>"),
    v3pa_slider,
    row(v3pa_input, apa_input),
    pa_help_div,
    Div(text="<b>Visibility window</b>"),
    row(visibility_date_input, visibility_btn),
    visibility_div,
    Div(text="<b>Optimize MSA pointing</b> "
             "<small>(grid search + refine, hMPT-derived)</small>"),
    Div(text=("<small style='color:#5a6b85; line-height:1.4'>"
              "Configure the search (method, ΔRA / ΔDec / ΔPA, "
              "collision protection, …) and run it from a single "
              "dialog. Results land in the same modal as before.</small>"),
        width=SIDEBAR_W - 20),
    opt_open_btn,
    width=SIDEBAR_W - 20,
))

pick_tab = TabPanel(title="Settings", child=column(
    Div(text="<b>Layers</b>"),
    layers_box,
    Div(text="<b>Slitlet</b>"),
    slitlet_select,
    snap_box,
    Div(text="<b>Overlay appearance</b>"),
    overlay_layer_select,
    overlay_alpha_slider,
    overlay_stroke_slider,
    Div(text="<b>Canvas</b>"),
    canvas_x_slider,
    canvas_y_slider,
    Div(text="<b>Customise display</b>"),
    stats_bar_open_btn,
    catalog_hover_open_btn,
    Div(text="<b>Actions</b>"),
    row(undo_btn, clear_btn),
    width=SIDEBAR_W - 20,
))

# MPT tab — everything to do with APT / MPT plans: import (from JSON,
# shutter CSV, or .aptx archive / program ID), the session save/load
# round-trip for collaboration, and the eMPT export bundle.
mpt_tab = TabPanel(title="MPT", child=column(
    Div(text="<b>Import a plan</b>"),
    _wrap_path_picker(
        mpt_json_path_input, mpt_json_browse_btn,
        header_html="<small>From a single MPT JSON</small>",
    ),
    mpt_plan_select,
    mpt_load_btn,
    Div(text="<hr style='border:none; border-top:1px dashed #d8dee8; "
             "margin:8px 0'/>"),
    _wrap_path_picker(
        mpt_csv_path_input, mpt_csv_browse_btn,
        header_html="<small>Or a shutter CSV (open-mask only)</small>",
    ),
    mpt_csv_load_btn,
    Div(text="<hr style='border:none; border-top:1px dashed #d8dee8; "
             "margin:8px 0'/>"),
    _wrap_path_picker(
        apt_path_input, apt_path_browse_btn,
        header_html="<small>Or straight from an APT <code>.aptx</code> file or program ID</small>",
    ),
    apt_program_input,
    apt_fetch_btn,
    apt_plan_select,
    apt_load_btn,
    Div(text="<hr style='border:none; border-top:2px solid #d8dee8; "
             "margin:12px 0 6px 0'/>"),
    Div(text="<b>Save / share session</b>"),
    _wrap_path_picker(
        session_save_path_input, session_save_browse_btn,
        header_html="<small>Save destination</small>",
    ),
    session_save_btn,
    _wrap_path_picker(
        session_load_path_input, session_load_browse_btn,
        header_html="<small>Load from session file</small>",
    ),
    session_load_btn,
    Div(text="<hr style='border:none; border-top:2px solid #d8dee8; "
             "margin:12px 0 6px 0'/>"),
    Div(text=f"<b>Export to APT</b> "
             f"<small>(eMPT bundle + <code>{MPT_PLAN_FILENAME}</code>)</small>"),
    _wrap_path_picker(
        export_dir_input, export_dir_browse_btn,
        header_html="<small>Output directory</small>",
    ),
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
    width=SIDEBAR_W,
    height_policy="max",
    styles={
        "overflow-y": "auto",
        # Leave room at the bottom for the position-fixed status bar
        # (42 px tall) so the last widget in any tab isn't covered.
        "max-height": "calc(100vh - 46px)",
        "padding-bottom": "8px",
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
curdoc().add_root(cat_edit_modal_backdrop)
curdoc().add_root(cat_edit_modal_card)
curdoc().add_root(cat_constraints_modal_backdrop)
curdoc().add_root(cat_constraints_modal_card)
curdoc().add_root(stats_bar_modal_backdrop)
curdoc().add_root(stats_bar_modal_card)
curdoc().add_root(catalog_hover_modal_backdrop)
curdoc().add_root(catalog_hover_modal_card)
# Status bar — separate root so its position:fixed style escapes the
# sidebar's scrollable container. Lives at the bottom-left of the
# viewport, under the sidebar.
curdoc().add_root(status)
curdoc().title = "vMPT — visual MSA Planning Tool"


# ── CLI auto-load ─────────────────────────────────────────────────────
# Args forwarded by `run.sh --args ...` arrive in sys.argv. We pre-fill
# the relevant path inputs; the existing on_change handlers do the
# actual loading. Loads are deferred to the next IO tick so the
# document is fully wired before we trigger heavy work.
def _parse_startup_args(argv: list[str]) -> dict:
    """Tolerant parser for `--fits`, `--jpg`, `--wcs`, `--catalog` (repeatable).

    Unknown args are silently ignored — Bokeh prefixes some of its own
    flags before --args and we don't want to throw on them."""
    out: dict = {"fits": None, "jpg": None, "wcs": None, "catalogs": []}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--fits", "--jpg", "--wcs") and i + 1 < len(argv):
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
