"""Tests for vmpt.overlays_io — DS9 .reg + .con display-overlay parsers."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from vmpt.overlays_io import (
    classify_overlay_file,
    detect_contour_coordsys,
    load_ds9_contours,
    load_ds9_regions,
)

regions = pytest.importorskip("regions")  # parser depends on the regions pkg


def _wcs(n=400, scale_arcsec=1.0):
    """A simple TAN WCS centred on (150, 2) deg with `scale_arcsec`/pixel."""
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w, n


_REG_TEXT = """# Region file format: DS9 version 4.1
fk5
circle(150.0, 2.0, 5")
ellipse(150.0, 2.0, 6", 3", 30)
box(150.0, 2.0, 8", 4", 0)
polygon(150.0, 2.0, 150.001, 2.0, 150.001, 2.001)
line(150.0, 2.0, 150.002, 2.002)
point(150.0, 2.0)
"""


def test_load_ds9_regions_shapes(tmp_path):
    w, n = _wcs()
    p = tmp_path / "r.reg"
    p.write_text(_REG_TEXT)
    out = load_ds9_regions(str(p), w)
    # circle + ellipse + box + polygon + line = 5 polylines; point = 1 marker.
    assert len(out["lines"]) == 5
    assert len(out["points"]) == 1
    assert out["n"] == 6
    # The point sits at CRVAL → centre pixel ((n/2)-0.5 in 0-indexed).
    px, py = out["points"][0]
    assert abs(px - (n / 2 - 0.5)) < 0.5
    assert abs(py - (n / 2 - 0.5)) < 0.5


def test_load_ds9_regions_circle_radius_in_pixels(tmp_path):
    """A 5" circle on a 1"/pix WCS should sample a ring ~5 px in radius
    around the centre pixel."""
    w, n = _wcs(scale_arcsec=1.0)
    p = tmp_path / "circ.reg"
    p.write_text("fk5\ncircle(150.0, 2.0, 5\")\n")
    out = load_ds9_regions(str(p), w)
    xs, ys = out["lines"][0]
    cx, cy = n / 2 - 0.5, n / 2 - 0.5
    r = np.hypot(np.asarray(xs) - cx, np.asarray(ys) - cy)
    assert np.allclose(r, 5.0, atol=0.2)


def test_load_ds9_regions_factor_scales_pixel_regions(tmp_path):
    """A native image-coordinate region is divided by the downsample factor."""
    w, _ = _wcs()
    p = tmp_path / "imgreg.reg"
    p.write_text("image\ncircle(100, 200, 20)\n")
    full = load_ds9_regions(str(p), w, factor=1)
    half = load_ds9_regions(str(p), w, factor=2)
    fx, fy = full["lines"][0]
    hx, hy = half["lines"][0]
    # Every vertex halves when factor=2.
    assert np.allclose(np.asarray(hx) * 2, np.asarray(fx), atol=1e-6)
    assert np.allclose(np.asarray(hy) * 2, np.asarray(fy), atol=1e-6)


def test_load_ds9_contours_sky(tmp_path):
    """Sky-coordinate .con: two segments separated by a blank line."""
    w, n = _wcs()
    p = tmp_path / "c.con"
    p.write_text(
        "150.0 2.0\n150.001 2.0\n150.001 2.001\n"
        "\n"
        "150.0 2.0\n149.999 1.999\n"
    )
    out = load_ds9_contours(str(p), w, coordsys="sky")
    assert out["n"] == 2
    assert len(out["lines"]) == 2
    # First vertex of segment 1 is CRVAL → centre pixel.
    xs, ys = out["lines"][0]
    assert abs(xs[0] - (n / 2 - 0.5)) < 0.5
    assert abs(ys[0] - (n / 2 - 0.5)) < 0.5


def test_load_ds9_contours_image_factor(tmp_path):
    """Image-coordinate .con divides vertices by the downsample factor and
    ignores comment lines."""
    w, _ = _wcs()
    p = tmp_path / "ci.con"
    p.write_text("# a comment\n100 200\n140 260\n")
    out = load_ds9_contours(str(p), w, coordsys="image", factor=2)
    assert out["n"] == 1
    xs, ys = out["lines"][0]
    assert xs == [50.0, 70.0]
    assert ys == [100.0, 130.0]


def test_classify_overlay_file_by_extension(tmp_path):
    r = tmp_path / "a.reg"; r.write_text("fk5\ncircle(1,2,3\")\n")
    c = tmp_path / "b.ctr"; c.write_text("icrs\nlevel=1\n(\n1 2\n3 4\n)\n")
    o = tmp_path / "c.con"; o.write_text("1 2\n3 4\n")
    assert classify_overlay_file(str(r)) == "region"
    assert classify_overlay_file(str(c)) == "contour"
    assert classify_overlay_file(str(o)) == "contour"


def test_classify_overlay_file_by_content(tmp_path):
    # Ambiguous .txt extensions → sniff the content.
    reg_txt = tmp_path / "shapes.txt"
    reg_txt.write_text("# Region file format: DS9 version 4.1\nfk5\n"
                       "polygon(1,2,3,4,5,6)\n")
    ctr_txt = tmp_path / "levels.txt"
    ctr_txt.write_text("# Contour file format: DS9 version 7.5\nicrs\nlevel=1\n")
    assert classify_overlay_file(str(reg_txt)) == "region"
    assert classify_overlay_file(str(ctr_txt)) == "contour"


def test_detect_contour_coordsys(tmp_path):
    sky = tmp_path / "s.ctr"
    sky.write_text("# Contour file format: DS9 version 7.5\nicrs\nlevel=1\n(\n")
    img = tmp_path / "i.ctr"
    img.write_text("# Contour file format: DS9 version 7.5\nimage\nlevel=1\n(\n")
    none = tmp_path / "n.ctr"
    none.write_text("1 2\n3 4\n")  # no frame line → default sky
    assert detect_contour_coordsys(str(sky)) == "sky"
    assert detect_contour_coordsys(str(img)) == "image"
    assert detect_contour_coordsys(str(none)) == "sky"


def test_load_ds9_contours_skips_degenerate_segments(tmp_path):
    """A single-vertex segment can't form a polyline and is dropped."""
    w, _ = _wcs()
    p = tmp_path / "deg.con"
    p.write_text("150.0 2.0\n\n150.0 2.0\n150.001 2.001\n")
    out = load_ds9_contours(str(p), w, coordsys="sky")
    assert out["n"] == 1  # the lone-vertex segment is skipped


# DS9 7.5 `.ctr` export: a frame line + `level=N` markers, each contour
# wrapped in `(` … `)`. The parentheses delimit segments — NOT blank lines.
_CTR_TEXT = """# Contour file format: DS9 version 7.5
# levels=( 1 )
global color=green width=1 dash=no dashlist=8 3
icrs
level=1
(
 150.0 2.0
 150.001 2.0
 150.001 2.001
)
level=1
(
 149.999 1.999
 149.998 1.998
)
"""


def test_load_ds9_contours_ctr_parenthesised(tmp_path):
    """The DS9 `.ctr` `( … )` layout must yield one segment per contour, not
    one giant merged polyline."""
    w, n = _wcs()
    p = tmp_path / "x.ctr"
    p.write_text(_CTR_TEXT)
    out = load_ds9_contours(str(p), w, coordsys="sky")
    assert out["n"] == 2, out["n"]            # two separate contours
    assert [len(xs) for xs, ys in out["lines"]] == [3, 2]
    # First vertex of contour 1 is at CRVAL → centre pixel.
    xs, ys = out["lines"][0]
    assert abs(xs[0] - (n / 2 - 0.5)) < 0.5
    assert abs(ys[0] - (n / 2 - 0.5)) < 0.5


def test_parse_contour_segments_handles_both_layouts(tmp_path):
    """`_parse_contour_segments` splits on `( )` for .ctr and on blank lines
    for plain .con, skipping metadata either way."""
    from vmpt.overlays_io import _parse_contour_segments
    ctr = tmp_path / "a.ctr"
    ctr.write_text(_CTR_TEXT)
    segs = _parse_contour_segments(str(ctr))
    assert [len(s) for s in segs] == [3, 2]
    con = tmp_path / "b.con"
    con.write_text("1 2\n3 4\n\n5 6\n7 8\n9 10\n")
    segs2 = _parse_contour_segments(str(con))
    assert [len(s) for s in segs2] == [2, 3]
