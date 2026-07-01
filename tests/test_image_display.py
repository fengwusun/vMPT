"""Image-display dialog + re-stretch behaviour (v1.8.0 UX follow-ups).

Covers three user-reported fixes:
  1. Changing brightness / contrast / stretch must NOT reset the zoom —
     `_restretch_image_glyph` swaps only the RGBA, leaving the ranges alone.
  2. The image-display controls live in a dialog (Settings tab), not inline.
  3. The DS9 region / contour loaders live in the "Load Add-on" dialog
     (Input tab), not inline.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402
from vmpt.image_io import load_fits  # noqa: E402


def _wcs_header(n):
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2, n / 2]
    w.wcs.cdelt = [-1 / 3600.0, 1 / 3600.0]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w.to_header()


def _load_small_fits():
    td = tempfile.mkdtemp()
    n = 64
    p = Path(td) / "i.fits"
    fits.writeto(str(p), np.random.RandomState(0).rand(n, n).astype(np.float32),
                 _wcs_header(n), overwrite=True)
    return load_fits(str(p))


def _load_big_fits(n=1600, max_dim=400):
    td = tempfile.mkdtemp()
    p = Path(td) / "big.fits"
    data = (np.arange(n * n).reshape(n, n) % 331).astype(np.float32)
    fits.writeto(str(p), data, _wcs_header(n), overwrite=True)
    return load_fits(str(p), max_dim=max_dim)


def test_restretch_preserves_zoom():
    """Brightness / contrast / stretch changes keep the current pan & zoom."""
    m.state["image"] = _load_small_fits()
    m.refresh_image_glyph()          # fits ranges to the frame
    # Simulate the user zooming in.
    m.fig.x_range.update(start=10.0, end=30.0)
    m.fig.y_range.update(start=12.0, end=28.0)
    before = (m.fig.x_range.start, m.fig.x_range.end,
              m.fig.y_range.start, m.fig.y_range.end)

    m._on_rgb_brightness("value", 0.0, 0.3)
    m._on_rgb_contrast("value", 1.0, 1.6)
    m._on_image_stretch("value", "asinh", "log")
    m._on_image_scale_mode("value", "percentile", "zscale")

    after = (m.fig.x_range.start, m.fig.x_range.end,
             m.fig.y_range.start, m.fig.y_range.end)
    assert before == after, f"re-stretch moved the view: {before} -> {after}"


def test_restretch_updates_the_pixels():
    """A re-stretch must actually change the displayed RGBA (else it's a
    no-op), while leaving the image placement (x/y/dw/dh) intact."""
    m.state["image"] = _load_small_fits()
    m.refresh_image_glyph()
    before = np.array(m.src_image.data["image"][0], copy=True)
    placement = (m.src_image.data["x"][0], m.src_image.data["y"][0],
                 m.src_image.data["dw"][0], m.src_image.data["dh"][0])
    m._on_image_stretch("value", "asinh", "linear")
    after = np.array(m.src_image.data["image"][0], copy=True)
    assert not np.array_equal(before, after), "stretch change had no effect"
    placement2 = (m.src_image.data["x"][0], m.src_image.data["y"][0],
                  m.src_image.data["dw"][0], m.src_image.data["dh"][0])
    assert placement == placement2, "re-stretch moved/resized the image glyph"


def test_image_display_dialog_owns_controls_and_toggles():
    """The FITS + RGB control groups live in the image-display modal card;
    the dialog opens and closes."""
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)

    card_ids = {id(x) for x in _walk(m.image_display_modal_card)}
    assert id(m._img_fits_group) in card_ids
    assert id(m._img_rgb_group) in card_ids

    m._open_image_display_modal()
    assert m.image_display_modal_card.visible
    assert m.image_display_modal_backdrop.visible
    m._close_image_display_modal()
    assert not m.image_display_modal_card.visible


def test_load_addon_dialog_owns_overlay_loaders_and_toggles():
    """The DS9 region + contour loaders live in the Load-Add-on modal card;
    the dialog opens and closes."""
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)

    card_ids = {id(x) for x in _walk(m.load_addon_modal_card)}
    assert id(m.addon_add_btn) in card_ids
    assert id(m.addon_clear_btn) in card_ids
    # Per-file management moved to the sidebar "Loaded add-ons" list, NOT the
    # dialog.
    assert id(m.sidebar_overlay_list_column) not in card_ids
    assert m.sidebar_overlay_section is not None

    m._open_load_addon_modal()
    assert m.load_addon_modal_card.visible
    assert m.load_addon_modal_backdrop.visible
    m._close_load_addon_modal()
    assert not m.load_addon_modal_card.visible


def _fits_with_wcs(tmp_path, n=400):
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2, n / 2]
    w.wcs.cdelt = [-1 / 3600.0, 1 / 3600.0]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    ip = tmp_path / "img.fits"
    fits.writeto(str(ip), np.zeros((n, n), np.float32), w.to_header(),
                 overwrite=True)
    return load_fits(str(ip))


def _write_reg_ctr(tmp_path):
    reg = tmp_path / "r.reg"
    reg.write_text("fk5\ncircle(150.0, 2.0, 5\")\n")
    ctr = tmp_path / "c.ctr"  # DS9 7.5 contour, icrs frame
    ctr.write_text("# Contour file format: DS9 version 7.5\nicrs\nlevel=1\n(\n"
                   " 150.0 2.0\n 150.001 2.0\n 150.001 2.001\n)\n")
    return reg, ctr


def test_unified_add_overlay_files_classifies_and_enables(tmp_path):
    """One list of mixed files becomes per-file overlay entries, with a
    contour's coordinate frame auto-detected, and both layer groups enabled."""
    m.state["image"] = _fits_with_wcs(tmp_path)
    m.state["overlays"] = []
    m.layers_box.active = [0, 1, 2]  # simulate < v1.8.0 prefs
    reg, ctr = _write_reg_ctr(tmp_path)
    m._add_overlay_files([str(reg), str(ctr)])

    ovs = m.state["overlays"]
    assert [o["kind"] for o in ovs] == ["region", "contour"]
    assert ovs[0]["path"] == str(reg) and ovs[0]["enabled"] is True
    assert ovs[1]["path"] == str(ctr) and ovs[1]["coordsys"] == "sky"  # detected
    # Both layer groups auto-enabled and their glyphs shown.
    assert 3 in m.layers_box.active and 4 in m.layers_box.active
    assert m.regions_glyph.visible and m.contours_glyph.visible
    # Clear-all empties everything.
    m._on_addon_clear()
    assert m.state["overlays"] == []


def test_per_file_toggle_hides_one_overlay(tmp_path):
    """Un-ticking a file in the toggle list removes only its geometry."""
    img = _fits_with_wcs(tmp_path)
    m.state["image"] = img
    # Pin the view to the whole image so the overlay LOD keeps every segment
    # (off-view segments are culled — irrelevant here).
    m.fig.x_range.update(start=0, end=img.data.shape[1])
    m.fig.y_range.update(start=0, end=img.data.shape[0])
    m.state["overlays"] = []
    m.layers_box.active = [0, 1, 2, 3, 4]
    reg, ctr = _write_reg_ctr(tmp_path)
    m._add_overlay_files([str(reg), str(ctr)])
    # Both present: 1 region polyline, 1 contour polyline.
    assert len(m.src_regions.data["xs"]) == 1
    assert len(m.src_contours.data["xs"]) == 1
    # The sidebar "Loaded add-ons" list has one row per overlay.
    assert len(m.sidebar_overlay_list_column.children) == 2
    # Un-tick the region (via its own per-row handler) → geometry drops.
    reg_entry = m.state["overlays"][0]
    m._make_addon_toggle(reg_entry)("active", [0], [])
    assert reg_entry["enabled"] is False
    assert len(m.src_regions.data["xs"]) == 0
    assert len(m.src_contours.data["xs"]) == 1


def test_per_file_overlay_colour_and_fill_reach_the_cds(tmp_path):
    """Each overlay's colour + fill alpha are painted into per-segment CDS
    columns, so different files render in different colours."""
    img = _fits_with_wcs(tmp_path)
    m.state["image"] = img
    m.fig.x_range.update(start=0, end=img.data.shape[1])   # whole image in view
    m.fig.y_range.update(start=0, end=img.data.shape[0])
    m.state["overlays"] = []
    m.layers_box.active = [0, 1, 2, 3, 4]
    reg, ctr = _write_reg_ctr(tmp_path)
    m._add_overlay_files([str(reg), str(ctr)])
    reg_entry = next(o for o in m.state["overlays"] if o["kind"] == "region")
    # Recolour + shade the region via its per-row handlers.
    m._make_addon_color(reg_entry)("color", None, "#00ff88")
    m._make_addon_fill(reg_entry)("value", 0.0, 0.6)
    assert m.src_regions.data["line_color"][0] == "#00ff88"
    assert m.src_regions.data["fill_color"][0] == "#00ff88"
    assert abs(m.src_regions.data["fill_alpha"][0] - 0.6) < 1e-9
    # The contour keeps its own (default) colour — independent per file.
    assert m.src_contours.data["line_color"][0] == m.DS9_CONTOUR_COLOR
    m.state["overlays"] = []
    m._rebuild_overlay_layers_impl()


def test_overlay_lod_decimates_out_culls_in_full_detail_when_zoomed():
    """A dense contour is level-of-detailed: zoomed out its vertices are
    decimated; zoomed hard onto one segment the others are culled and the
    visible one keeps full detail. Geometry set directly in scaled-image px."""
    import numpy as np
    xl = np.linspace(100.0, 300.0, 4000)          # dense LEFT segment
    xr = np.linspace(3700.0, 3900.0, 4000)         # dense RIGHT segment
    yl = np.full(4000, 200.0)
    m.state["_ovl_full"] = {
        "r_xs": [], "r_ys": [], "r_lc": [], "r_fa": [],
        "p_x": [], "p_y": [], "p_lc": [],
        "c_xs": [xl.tolist(), xr.tolist()],
        "c_ys": [yl.tolist(), yl.tolist()],
        "c_lc": ["#22d3ee", "#22d3ee"], "c_fa": [0.0, 0.0],
    }
    m.state["_ovl_med_spacing"] = m._median_vertex_spacing([xl, xr], [yl, yl])
    m.state["frame_x"] = 800

    # Zoomed OUT (whole 4000 px field): both segments visible but decimated.
    m.fig.x_range.update(start=0, end=4000)
    m.fig.y_range.update(start=-100, end=500)
    m._apply_overlay_lod()
    assert len(m.src_contours.data["xs"]) == 2
    n_out = sum(len(xs) for xs in m.src_contours.data["xs"])
    assert n_out < 500, f"expected heavy decimation, got {n_out} of 8000"

    # Zoomed HARD onto the LEFT segment: right one culled, left at full detail.
    m.fig.x_range.update(start=199.0, end=201.0)
    m.fig.y_range.update(start=195.0, end=205.0)
    m._apply_overlay_lod()
    assert len(m.src_contours.data["xs"]) == 1            # right segment culled
    assert len(m.src_contours.data["xs"][0]) == 4000      # left kept in full

    m.state["_ovl_full"] = None
    m._apply_overlay_lod()


def test_dense_overlay_hidden_while_moving_restored_on_settle():
    """A DENSE overlay is hidden during a pan/zoom (so the drag redraws only
    the image) and re-shown once motion settles; a SPARSE overlay stays drawn."""
    m.state["_ovl_full"] = {
        "r_xs": [], "r_ys": [], "r_lc": [], "r_fa": [],
        "p_x": [], "p_y": [], "p_lc": [],
        "c_xs": [list(range(10))], "c_ys": [[0.0] * 10],
        "c_lc": ["#22d3ee"], "c_fa": [0.0],
    }
    m.state["_ovl_med_spacing"] = 1.0
    m.layers_box.active = [0, 1, 2, 3, 4]        # contours layer ON
    m.contours_glyph.visible = True
    m._ovl_interaction["hidden"] = False

    # Dense → hidden on the first move.
    m.state["_ovl_nverts"] = 50000
    m._hide_overlays_during_interaction()
    assert m._ovl_interaction["hidden"] is True
    assert m.contours_glyph.visible is False
    # Settle → re-shown.
    m._run_ovl_lod()
    assert m._ovl_interaction["hidden"] is False
    assert m.contours_glyph.visible is True

    # Sparse → never hidden (cheap to keep drawing; no flicker).
    m.state["_ovl_nverts"] = 120
    m.contours_glyph.visible = True
    m._hide_overlays_during_interaction()
    assert m._ovl_interaction["hidden"] is False
    assert m.contours_glyph.visible is True

    m.state["_ovl_full"] = None
    m._apply_overlay_lod()


# ---------------------------------------------------------------------------
# On-demand LOD (zoom-dependent FITS resolution)
# ---------------------------------------------------------------------------


def test_lod_sharpens_on_zoom_and_clears_on_zoomout():
    # Base tier (800 px) matches the canvas target, mirroring the real app
    # where DEFAULT_FITS_MAX_DIM (4000) ≫ the ~800 px canvas — so a full-image
    # view needs no sharpening and only zooming in does.
    m.state["frame_x"] = 800
    img = _load_big_fits(n=3200, max_dim=800)   # factor 4, base 800 px
    assert img.factor == 4
    m.state["image"] = img
    m.refresh_image_glyph()
    Hb, Wb = img.shape
    m.fig.x_range.update(start=0, end=Wb)
    m.fig.y_range.update(start=0, end=Hb)
    m._render_fits_lod()
    assert not m.src_lod.data.get("image"), "no LOD when base tier matches screen"

    # Zoom into a small central region → a sharper crop appears (up to native).
    m.fig.x_range.update(start=Wb / 2 - 40, end=Wb / 2 + 40)
    m.fig.y_range.update(start=Hb / 2 - 40, end=Hb / 2 + 40)
    m._render_fits_lod()
    assert m.src_lod.data.get("image"), "LOD crop should populate when zoomed in"
    assert m.state["_lod_tier"] == 1, "deep zoom should reach native tier"
    lod = m.src_lod.data
    base_ppu = m.src_image.data["image"][0].shape[1] / Wb
    crop_ppu = lod["image"][0].shape[1] / lod["dw"][0]
    assert crop_ppu > base_ppu, "crop not sharper than base tier"

    # Zoom back out → LOD clears (base tier only).
    m.fig.x_range.update(start=0, end=Wb)
    m.fig.y_range.update(start=0, end=Hb)
    m._render_fits_lod()
    assert not m.src_lod.data.get("image"), "LOD should clear when zoomed out"


def test_lod_glyph_not_pinned_to_ranges():
    """The LOD overlay must NOT drive the auto-range, or sharpening the view
    would make the figure zoom to the crop."""
    assert m.lod_glyph not in (m.fig.x_range.renderers or [])
    assert m.lod_glyph not in (m.fig.y_range.renderers or [])


def test_small_fits_never_uses_lod():
    """A FITS below the display cap loads at native res (factor 1) → no LOD."""
    img = _load_small_fits()
    assert img.factor == 1
    m.state["image"] = img
    m.refresh_image_glyph()
    m.fig.x_range.update(start=10, end=30)
    m.fig.y_range.update(start=10, end=30)
    m._render_fits_lod()
    assert not m.src_lod.data.get("image")


# ---------------------------------------------------------------------------
# Colormap wiring + Layers / canvas dialog moves
# ---------------------------------------------------------------------------


def test_colormap_and_invert_handlers_restretch():
    img = _load_small_fits()
    m.state["image"] = img
    m._select_cmap("gray")
    m.state["image_invert"] = False
    m.refresh_image_glyph()
    before = np.array(m.src_image.data["image"][0], copy=True)
    m._select_cmap("viridis")
    assert m.state["image_cmap"] == "viridis"
    after = np.array(m.src_image.data["image"][0], copy=True)
    assert not np.array_equal(before, after), "colormap change didn't re-render"
    m._on_image_invert("active", [], [0])
    assert m.state["image_invert"] is True
    m._select_cmap("gray")
    m._on_image_invert("active", [0], [])   # restore


def test_nan_white_renders_blank_pixels_white_and_persists():
    """Toggling 'Render NaN as white' paints blank FITS pixels white and the
    choice round-trips through prefs (unlike the per-session colormap)."""
    td = tempfile.mkdtemp()
    n = 64
    arr = np.random.RandomState(0).rand(n, n).astype(np.float32)
    arr[:8, :8] = np.nan                       # an 8×8 blank corner
    p = Path(td) / "nan.fits"
    fits.writeto(str(p), arr, _wcs_header(n), overwrite=True)
    img = load_fits(str(p))
    assert np.isnan(img.data).sum() == 64      # NaN survives the load

    m.state["image"] = img
    m._select_cmap("gray")
    m.state["image_invert"] = False

    m.state["image_nan_white"] = False
    m.image_nan_white_checkbox.active = []
    m.refresh_image_glyph()
    off = np.array(m.src_image.data["image"][0], copy=True)
    white_off = int(np.count_nonzero(off == 0xFFFFFFFF))

    m.image_nan_white_checkbox.active = [0]    # tick the box …
    m._on_image_nan_white("active", [], [0])   # … + its on_change handler
    assert m.state["image_nan_white"] is True
    on = np.array(m.src_image.data["image"][0], copy=True)
    white_on = int(np.count_nonzero(on == 0xFFFFFFFF))
    assert white_on >= 64                      # the blank corner is now white
    assert white_on > white_off

    # Persisted, unlike the FITS colormap/invert.
    assert m._collect_prefs()["image_nan_white"] is True
    m.state["image_nan_white"] = False
    m.image_nan_white_checkbox.active = []
    m._apply_prefs({"image_nan_white": True})
    assert m.state["image_nan_white"] is True
    assert list(m.image_nan_white_checkbox.active) == [0]

    # The control lives in the Image display dialog.
    assert id(m.image_nan_white_checkbox) in _image_dialog_ids()

    m._on_image_nan_white("active", [0], [])   # restore off
    m.state["image"] = None


def test_layers_dialog_owns_layers_and_appearance():
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)

    ids = {id(x) for x in _walk(m.layers_modal_card)}
    assert id(m.layers_box) in ids
    assert id(m.overlay_layer_select) in ids
    assert id(m.overlay_alpha_slider) in ids
    assert id(m.overlay_stroke_slider) in ids
    m._open_layers_modal()
    assert m.layers_modal_card.visible
    m._close_layers_modal()
    assert not m.layers_modal_card.visible


def test_rich_colormap_dropdown_swatches():
    """The dropdown has one gradient-backed button per colormap (name + a
    0→1 colour bar), and picking one repaints the current-selection swatch."""
    # 10 option buttons, each carrying its colormap's gradient in a stylesheet.
    assert len(m._cmap_option_btns) == len(m.FITS_COLORMAPS)
    for name, b in zip(m.FITS_COLORMAPS, m._cmap_option_btns):
        assert b.label == name.capitalize()
        assert "linear-gradient" in b.stylesheets[0].css

    # Gradient helper reverses under invert (viridis dark→yellow vs yellow→dark).
    non_inv = m._cmap_gradient_css("viridis", invert=False)
    inv = m._cmap_gradient_css("viridis", invert=True)
    assert non_inv != inv
    assert non_inv.index("rgb(253") > inv.index("rgb(253")

    # Selecting a map repaints the current button and collapses the panel.
    m._select_cmap("magma")
    assert m.state["image_cmap"] == "magma"
    assert m.cmap_current_btn.label.startswith("Magma")
    assert "linear-gradient" in m.cmap_current_btn.stylesheets[0].css
    assert not m.cmap_options_col.visible
    assert id(m.cmap_current_btn) in _image_dialog_ids()
    assert id(m.cmap_options_col) in _image_dialog_ids()
    m._select_cmap("gray")


def test_pixel_histogram_fits_and_restretch(tmp_path):
    """Loading a FITS builds the histogram bars + vmin/vmax box + stretch
    curve; a stretch change updates the curve/box but NOT the bars (cached)."""
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    n = 200
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2, n / 2]
    w.wcs.cdelt = [-1 / 3600.0, 1 / 3600.0]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    data = np.random.RandomState(3).normal(100, 10, (n, n)).astype(np.float32)
    p = tmp_path / "h.fits"
    fits.writeto(str(p), data, w.to_header(), overwrite=True)
    m.state["image"] = load_fits(str(p))
    m.state["image_stretch"] = "asinh"      # known starting stretch
    m.refresh_image_glyph()

    # Full data range on x, padded by ONE BIN each side for easier dragging;
    # log y-axis (bars = raw counts, floor 0.5).
    assert 0 < len(m.src_hist.data["left"]) <= 50
    assert m.hist_fig.y_scale.__class__.__name__ == "LogScale"
    bw = (float(data.max()) - float(data.min())) / 50.0
    assert abs(m.hist_fig.x_range.start - (float(data.min()) - bw)) < 1e-2
    assert abs(m.hist_fig.x_range.end - (float(data.max()) + bw)) < 1e-2
    assert all(b == 0.5 for b in m.src_hist.data["bottom"])
    assert m.hist_range_box.visible
    assert (m.hist_fig.x_range.start <= m.hist_range_box.left
            <= m.hist_range_box.right <= m.hist_fig.x_range.end)
    # Two handles at OPPOSITE ends: vmin △ (up) low at the bottom, vmax ▽
    # (down) high at the top.
    assert len(m.src_vhandles.data["x"]) == 2
    assert list(m.src_vhandles.data["mkr"]) == ["triangle", "inverted_triangle"]
    _yh = list(m.src_vhandles.data["y"])
    assert _yh[0] < _yh[1]
    curve = list(m.src_stretch_curve.data["y"])
    assert curve and curve[0] < 0.05 and abs(curve[-1] - 1.0) < 1e-6

    bars = list(m.src_hist.data["top"])
    data_ref = m.state["_hist_data_ref"]
    m._on_image_stretch("value", "asinh", "linear")    # asinh → linear
    assert m.state["_hist_data_ref"] is data_ref          # bars not recomputed
    assert m.src_hist.data["top"] == bars
    assert m.src_stretch_curve.data["y"] != curve         # curve did change

    # Drag the handles → adopt a manual [vmin, vmax] and re-stretch.
    lo, hi = m.hist_fig.x_range.start, m.hist_fig.x_range.end
    nmin, nmax = lo + 0.2 * (hi - lo), lo + 0.6 * (hi - lo)
    m.src_vhandles.data = dict(x=[nmin, nmax], y=list(m.src_vhandles.data["y"]),
                               mkr=list(m.src_vhandles.data["mkr"]))
    assert m.state["image_scale_mode"] == "manual"
    assert abs(m.state["image_vmin"] - nmin) < 1e-6
    assert abs(m.state["image_vmax"] - nmax) < 1e-6
    assert abs(m.hist_range_box.left - nmin) < 1e-4


def test_pixel_histogram_rgb_hides_vmin_vmax_box(tmp_path):
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    n = 64
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2, n / 2]
    w.wcs.cdelt = [-1 / 3600.0, 1 / 3600.0]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hdr = w.to_header()
    hdr["NAXIS1"] = n
    hdr["NAXIS2"] = n
    sidecar = tmp_path / "wcs.fits"
    fits.PrimaryHDU(header=hdr).writeto(str(sidecar), overwrite=True,
                                       output_verify="ignore")
    from PIL import Image
    jpg = tmp_path / "rgb.jpg"
    Image.fromarray(np.random.RandomState(4).randint(
        0, 256, (n, n, 3), dtype=np.uint8), "RGB").save(str(jpg))
    from vmpt.image_io import load_jpg_with_sidecar
    m.state["image"] = load_jpg_with_sidecar(str(jpg), str(sidecar))
    m.refresh_image_glyph()
    assert 0 < len(m.src_hist.data["left"]) <= 50
    assert not m.hist_range_box.visible          # RGB: no vmin/vmax box
    assert not m.src_vhandles.data["x"]          # …and no drag handles
    assert "RGB luminance" in m.hist_info_div.text


def test_histogram_figure_in_image_dialog():
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs", "renderers"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)
    ids = {id(x) for x in _walk(m.image_display_modal_card)}
    assert id(m.hist_fig) in ids
    assert id(m.hist_info_div) in ids


def test_tip_box_navigation():
    """The tip box ‹ / › arrows skim tips manually, update the counter, wrap,
    and pause the auto-rotation briefly after a click."""
    n = len(m._TIPS)
    assert n > 20  # includes the v1.8.0 tips
    m._show_tip(0)
    assert m.tip_counter_div.text == f"1 / {n}"
    m._tip_next()
    assert m._tip_state["idx"] == 1 and m._tip_state["skip"] == 4
    assert m.tip_counter_div.text == f"2 / {n}"
    m._show_tip(0)
    m._tip_prev()                      # wraps to the last tip
    assert m._tip_state["idx"] == n - 1
    # Auto-advance is skipped while a manual-nav pause is active.
    m._tip_state["skip"] = 1
    idx = m._tip_state["idx"]
    m._advance_tip()
    assert m._tip_state["idx"] == idx and m._tip_state["skip"] == 0
    m._advance_tip()                   # pause elapsed → advances
    assert m._tip_state["idx"] == (idx + 1) % n


def test_tips_include_v18_features():
    labels = [t[1] for t in m._TIPS]
    for lab in ("Image display dialog", "Pixel histogram",
                "DS9 regions & contours", "Per-item colour & fill"):
        assert lab in labels, lab


def test_cli_addon_flag_parsed():
    """`--addon` (repeatable) collects DS9 overlay paths from the CLI."""
    args = m._parse_startup_args(
        ["--fits", "x.fits", "--addon", "a.reg", "--addon", "b.ctr",
         "--catalog", "c.csv"])
    assert args["addons"] == ["a.reg", "b.ctr"]
    assert args["fits"] == "x.fits" and args["catalogs"] == ["c.csv"]


def test_fits_colormap_default_is_gray():
    """The FITS colormap must default to gray (not persisted across launches)."""
    prefs = m._collect_prefs()
    assert "image_cmap" not in prefs and "image_invert" not in prefs
    m._apply_prefs({"image_cmap": "viridis", "image_invert": True})
    assert m.state["image_cmap"] == "gray"
    assert m.state["image_invert"] is False


def _image_dialog_ids():
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)
    return {id(x) for x in _walk(m.image_display_modal_card)}


def test_canvas_spinners_moved_into_image_dialog():
    def _walk(node, seen=None):
        seen = seen or set()
        if id(node) in seen:
            return
        seen.add(id(node))
        yield node
        for attr in ("children", "child", "tabs"):
            v = getattr(node, attr, None)
            if v is None:
                continue
            for c in (v if isinstance(v, list) else [v]):
                yield from _walk(c, seen)

    ids = {id(x) for x in _walk(m.image_display_modal_card)}
    assert id(m.canvas_x_spinner) in ids
    assert id(m.canvas_y_spinner) in ids
    assert id(m.cmap_current_btn) in ids
    assert id(m.image_invert_checkbox) in ids
