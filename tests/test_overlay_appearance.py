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
