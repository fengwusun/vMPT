"""Every Settings → Overlay appearance layer must be adjustable.

Each entry in the `Adjust layer` Select has to be a real key in
`_OVERLAY_LAYER_CONFIG`, and moving the alpha / stroke sliders has to reach
the right place — the glyph attribute for scalar layers, or
`state["overlap_base_alpha_*"]` for the field-referenced MPT spec-overlap
layers (pink / orange / purple). Regression for the stale "Overlapping
shutters" Select entry, which mapped to no config so its sliders were inert.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402


def test_every_select_option_is_a_config_key():
    cfg = set(m._OVERLAY_LAYER_CONFIG)
    opts = set(m.overlay_layer_select.options)
    assert opts <= cfg, f"Select options with no config: {sorted(opts - cfg)}"
    # And every adjustable layer is reachable from the Select.
    assert cfg <= opts, f"Config layers missing from Select: {sorted(cfg - opts)}"


def test_alpha_and_stroke_sliders_drive_every_layer():
    for layer in m.overlay_layer_select.options:
        cfg = m._OVERLAY_LAYER_CONFIG[layer]
        m.overlay_layer_select.value = layer
        m._on_overlay_layer("value", None, layer)  # ensure slider ranges set

        a_lo, a_hi, _ = cfg["alpha_range"]
        s_lo, s_hi, _ = cfg["stroke_range"]
        a_val = round((a_lo + a_hi) / 2.0, 4)
        s_val = round((s_lo + s_hi) / 2.0, 4)

        # Alpha → glyph attr (scalar layers) or state (spec-overlap layers).
        m.overlay_alpha_slider.value = a_val
        if cfg.get("alpha_attr") is not None:
            got = float(getattr(cfg["glyph"].glyph, cfg["alpha_attr"]))
        else:
            got = float(m.state[cfg["alpha_state_key"]])
        assert abs(got - a_val) < 1e-6, (
            f"{layer!r}: alpha slider did not propagate "
            f"(got {got}, set {a_val})")

        # Stroke → glyph attr (line_width / marker size).
        m.overlay_stroke_slider.value = s_val
        got_s = float(getattr(cfg["glyph"].glyph, cfg["stroke_attr"]))
        assert abs(got_s - s_val) < 1e-6, (
            f"{layer!r}: stroke slider did not propagate "
            f"(got {got_s}, set {s_val})")


def test_every_layer_has_per_layer_defaults():
    """Reset restores each layer's OWN default, so every config entry must
    carry `default_alpha` + `default_stroke`. A blanket value (the old
    0.20 / 1.0) mangled most layers — e.g. the catalog marker size
    collapsing to 1 px and the stuck-open outline dimming to 0.20.
    """
    for name, cfg in m._OVERLAY_LAYER_CONFIG.items():
        assert "default_alpha" in cfg, f"{name}: missing default_alpha"
        assert "default_stroke" in cfg, f"{name}: missing default_stroke"
        # Defaults must sit within the layer's own slider range (except the
        # catalog alpha, which restores the field reference, not a scalar).
        da = cfg["default_alpha"]
        if isinstance(da, (int, float)):
            lo, hi, _ = cfg["alpha_range"]
            assert lo <= da <= hi, f"{name}: default_alpha {da} out of range"
        lo, hi, _ = cfg["stroke_range"]
        assert lo <= cfg["default_stroke"] <= hi, (
            f"{name}: default_stroke {cfg['default_stroke']} out of range")


def test_masked_layers_default_to_alpha_020_stroke_05():
    """User-specified defaults: masked / stuck / conflict = 0.20 alpha,
    0.5 stroke; the silver operable edge keeps 0.20 / 1.0."""
    for name in ("Mask Stuck (pink)", "Masked (overlapping warning)",
                 "Mask Conflict (purple)"):
        cfg = m._OVERLAY_LAYER_CONFIG[name]
        assert cfg["default_alpha"] == 0.20, name
        assert cfg["default_stroke"] == 0.5, name
    op = m._OVERLAY_LAYER_CONFIG["Operable shutters"]
    assert (op["default_alpha"], op["default_stroke"]) == (0.20, 1.0)


# ---------------------------------------------------------------------------
# v1.8.0 — colour + fill alpha are PER ITEM (catalog list / Load Add-on
# dialog), NOT per layer. The Layers dialog keeps only visibility + outline
# alpha + stroke; region/contour glyphs read colour + fill from CDS columns.
# ---------------------------------------------------------------------------
from bokeh.models import CheckboxGroup, ColorPicker, Slider  # noqa: E402


def test_overlay_glyphs_use_per_segment_colour_fields():
    """Region/contour glyphs read colour + fill alpha from per-segment CDS
    columns, so each loaded file can be coloured / shaded independently."""
    assert m.regions_glyph.glyph.line_color == "line_color"
    assert m.contours_glyph.glyph.line_color == "line_color"
    assert m.region_points_glyph.glyph.line_color == "line_color"
    for g in (m.regions_fill_glyph, m.contours_fill_glyph):
        assert g.glyph.fill_color == "fill_color"
        assert g.glyph.fill_alpha == "fill_alpha"
    for src in (m.src_regions, m.src_contours):
        assert {"line_color", "fill_color", "fill_alpha"} <= set(src.data)
    assert "line_color" in m.src_region_points.data
    # Fill overlay still must not drive the figure auto-range.
    assert m.regions_fill_glyph not in (m.fig.x_range.renderers or [])
    assert m.contours_fill_glyph not in (m.fig.y_range.renderers or [])


def test_layers_dialog_has_no_colour_or_fill_widgets():
    """Colour + fill moved out of Settings → Layers into per-item controls."""
    assert not hasattr(m, "overlay_color_picker")
    assert not hasattr(m, "overlay_fill_slider")
    # The 3 colour-capable layers no longer carry layer-wide colour/fill hooks.
    for name in ("Catalog sources", "DS9 regions", "Contours"):
        cfg = m._OVERLAY_LAYER_CONFIG[name]
        assert "color_set" not in cfg and "color_get" not in cfg
        assert "fill_glyph" not in cfg and "default_color" not in cfg


def _flatten(node):
    """Yield a layout node and all its descendant models (children/child)."""
    yield node
    for attr in ("children", "child"):
        v = getattr(node, attr, None)
        if v is None:
            continue
        for c in (v if isinstance(v, list) else [v]):
            yield from _flatten(c)


def test_sidebar_overlay_list_builds_swatch_checkbox_delete_rows():
    """Each overlay row in the sidebar 'Loaded add-ons' list is a compact
    [colour swatch button][on/off checkbox + name][× delete]; the section shows
    only when overlays exist. Colour + fill live in the popover, not the row."""
    from bokeh.models import Button
    saved = list(m.state.get("overlays", []))
    try:
        m.state["overlays"] = [
            {"path": "/x/a.reg", "kind": "region", "coordsys": "sky",
             "enabled": True, "color": "#00ff88", "fill_alpha": 0.4},
            {"path": "/x/b.ctr", "kind": "contour", "coordsys": "sky",
             "enabled": False, "color": "#ffcc00", "fill_alpha": 0.0},
        ]
        m._update_addon_list()
        assert m.sidebar_overlay_section.visible is True
        rows = m.sidebar_overlay_list_column.children
        assert len(rows) == 2
        r0 = list(rows[0].children)
        btns = [w for w in r0 if isinstance(w, Button)]
        cb = [w for w in r0 if isinstance(w, CheckboxGroup)]
        assert len(btns) == 2                       # swatch + delete
        sw = btns[0]
        assert sw.width == sw.height                # square swatch button
        assert "#00ff88" in sw.stylesheets[0].css   # painted with its colour
        assert cb and cb[0].active == [0]
        # No inline colour picker / fill slider in the row anymore.
        assert not [w for w in _flatten(rows[0]) if isinstance(w, ColorPicker)]
        assert not [w for w in _flatten(rows[0]) if isinstance(w, Slider)]
        cb1 = [w for w in rows[1].children if isinstance(w, CheckboxGroup)]
        assert cb1 and cb1[0].active == []          # disabled → unticked
    finally:
        m.state["overlays"] = saved
        m._update_addon_list()


def test_overlay_style_popover_edits_colour_and_alpha():
    """Clicking a row's swatch opens the shared popover targeting that overlay;
    the popover's colour picker + fill-α slider edit it and repaint the swatch."""
    from bokeh.models import Button
    saved = list(m.state.get("overlays", []))
    try:
        o = {"path": "/x/a.ctr", "kind": "contour", "coordsys": "sky",
             "enabled": True, "color": "#22d3ee", "fill_alpha": 0.0}
        m.state["overlays"] = [o]
        m._update_addon_list()
        sw = [w for w in m.sidebar_overlay_list_column.children[0].children
              if isinstance(w, Button)][0]
        m._open_ovl_style(o, sw)                     # click the swatch
        assert m.ovl_style_modal_card.visible is True
        assert m._ovl_style_target["o"] is o
        assert m._ovl_style_color.color == "#22d3ee"
        assert abs(m._ovl_style_alpha.value - 0.0) < 1e-9
        # Edit colour + alpha in the popover.
        m._on_ovl_style_color("color", "#22d3ee", "#ff8800")
        m._on_ovl_style_alpha("value", 0.0, 0.6)
        assert o["color"] == "#ff8800"
        assert abs(o["fill_alpha"] - 0.6) < 1e-9
        assert "#ff8800" in sw.stylesheets[0].css    # swatch repainted
        m._close_ovl_style()
        assert m.ovl_style_modal_card.visible is False
        assert m._ovl_style_target["o"] is None
    finally:
        m.state["overlays"] = saved
        m._update_addon_list()


def test_sidebar_overlay_delete_removes_entry():
    """The × on a sidebar row removes that overlay and hides the section when
    the last one goes."""
    saved = list(m.state.get("overlays", []))
    try:
        o = {"path": "/x/only.reg", "kind": "region", "coordsys": "sky",
             "enabled": True, "color": "#00ff88", "fill_alpha": 0.0}
        m.state["overlays"] = [o]
        m._update_addon_list()
        assert m.sidebar_overlay_section.visible is True
        m._delete_overlay(o)()                       # click the × handler
        assert o not in m.state["overlays"]
        assert m.sidebar_overlay_section.visible is False
    finally:
        m.state["overlays"] = saved
        m._update_addon_list()


def test_addon_per_file_handlers_update_the_entry():
    """The per-row colour / fill / toggle handlers mutate their own overlay
    dict (they close over it, not an index)."""
    o = {"path": "/x/a.reg", "kind": "region", "coordsys": "sky",
         "enabled": True, "color": "#ff3b30", "fill_alpha": 0.0}
    saved = list(m.state.get("overlays", []))
    try:
        m.state["overlays"] = [o]
        m._make_addon_color(o)("color", None, "#123456")
        assert o["color"] == "#123456"
        m._make_addon_fill(o)("value", 0.0, 0.5)
        assert abs(o["fill_alpha"] - 0.5) < 1e-9
        m._make_addon_toggle(o)("active", [0], [])
        assert o["enabled"] is False
    finally:
        m.state["overlays"] = saved
        m._update_addon_list()
