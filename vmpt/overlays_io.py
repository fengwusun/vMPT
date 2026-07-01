"""Loaders for DS9 region (.reg) and contour (.con) display overlays.

Both return geometry already projected into the loaded image's **pixel /
figure** coordinates (the same frame catalog targets are plotted in), so the
caller just pushes the result into a Bokeh ColumnDataSource.

- ``load_ds9_regions`` parses a DS9 ``.reg`` with the astropy-affiliated
  ``regions`` package (full shape + coordinate-frame support) and samples each
  shape's outline into polyline vertices.
- ``load_ds9_contours`` parses a DS9 ``.con`` (whitespace ``x y`` per line,
  blank lines separating segments) in sky (RA/Dec degrees) or image-pixel
  coordinates.

A returned overlay is ``{"lines": [(xs, ys), …], "points": [(x, y), …]}`` where
each ``xs``/``ys`` is a list of floats (one polyline) and points are single
markers (DS9 ``point`` regions).
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

_N_SAMPLE = 72  # vertices used to sample a circle / ellipse outline

# DS9 coordinate-frame tokens that appear on their own line in a .reg / .ctr.
_SKY_FRAMES = {"fk5", "fk4", "icrs", "galactic", "ecliptic", "j2000", "b1950",
               "wcs", "wcs0"}
_IMAGE_FRAMES = {"image", "physical", "detector", "amplifier", "logical"}
_REGION_SHAPES = ("circle", "ellipse", "box", "polygon", "annulus", "point",
                  "line", "vector", "text", "rotbox")


def classify_overlay_file(path: str) -> str:
    """Return ``"region"`` or ``"contour"`` for a DS9 add-on file, by extension
    then by content — so one file picker can accept both. ``.reg`` / ``.regions``
    → region; ``.ctr`` / ``.con`` → contour; anything else is sniffed (a DS9
    region header or a ``shape(`` line ⇒ region; a contour header / ``level=`` /
    a leading ``(`` ⇒ contour), defaulting to contour."""
    ext = str(path).rsplit(".", 1)[-1].lower() if "." in str(path) else ""
    if ext in ("reg", "regions"):
        return "region"
    if ext in ("ctr", "con"):
        return "contour"
    try:
        with open(path) as f:
            head = "".join(f.readline() for _ in range(40)).lower()
    except OSError:
        return "contour"
    if "region file format" in head:
        return "region"
    if "contour file format" in head or "level=" in head:
        return "contour"
    for shape in _REGION_SHAPES:
        if f"{shape}(" in head:
            return "region"
    return "contour"


def detect_contour_coordsys(path: str) -> str:
    """Sniff a DS9 contour file's coordinate frame → ``"sky"`` or ``"image"``.

    DS9 writes the frame on its own line (e.g. ``icrs`` or ``image``). Sky
    frames (fk5/icrs/galactic/…) project through the WCS; pixel frames
    (image/physical/…) are used directly. Defaults to ``"sky"``."""
    try:
        with open(path) as f:
            for _ in range(40):
                tok = f.readline().strip().lower()
                if tok in _IMAGE_FRAMES:
                    return "image"
                if tok in _SKY_FRAMES:
                    return "sky"
    except OSError:
        pass
    return "sky"


def _ellipse_xy(cx, cy, a, b, ang_rad):
    t = np.linspace(0.0, 2.0 * np.pi, _N_SAMPLE)
    x0, y0 = a * np.cos(t), b * np.sin(t)
    ca, sa = np.cos(ang_rad), np.sin(ang_rad)
    return (cx + x0 * ca - y0 * sa).tolist(), (cy + x0 * sa + y0 * ca).tolist()


def _rect_xy(cx, cy, a, b, ang_rad):
    loc = np.array([(-a, -b), (a, -b), (a, b), (-a, b), (-a, -b)], dtype=float)
    ca, sa = np.cos(ang_rad), np.sin(ang_rad)
    xs = cx + loc[:, 0] * ca - loc[:, 1] * sa
    ys = cy + loc[:, 0] * sa + loc[:, 1] * ca
    return xs.tolist(), ys.tolist()


def _angle_rad(pr) -> float:
    ang = getattr(pr, "angle", 0.0)
    try:
        return float(ang.to_value("rad"))
    except AttributeError:
        return float(np.deg2rad(float(ang)))


def _pixel_region_paths(pr) -> tuple[list, list]:
    """A regions PixelRegion → (polylines, points), each in pixel coords."""
    lines: list = []
    points: list = []
    name = type(pr).__name__
    if name.startswith("Circle") and hasattr(pr, "radius") and "Annulus" not in name:
        lines.append(_ellipse_xy(pr.center.x, pr.center.y,
                                 float(pr.radius), float(pr.radius), 0.0))
    elif name.startswith("Ellipse"):
        lines.append(_ellipse_xy(pr.center.x, pr.center.y,
                                 float(pr.width) / 2.0, float(pr.height) / 2.0,
                                 _angle_rad(pr)))
    elif name.startswith("Rectangle"):
        lines.append(_rect_xy(pr.center.x, pr.center.y,
                              float(pr.width) / 2.0, float(pr.height) / 2.0,
                              _angle_rad(pr)))
    elif name.startswith("Polygon") and hasattr(pr, "vertices"):
        vx = np.asarray(pr.vertices.x, dtype=float)
        vy = np.asarray(pr.vertices.y, dtype=float)
        if len(vx):
            vx = np.append(vx, vx[0])  # close the loop
            vy = np.append(vy, vy[0])
        lines.append((vx.tolist(), vy.tolist()))
    elif name.startswith("Line") and hasattr(pr, "start"):
        lines.append(([pr.start.x, pr.end.x], [pr.start.y, pr.end.y]))
    elif "Annulus" in name and hasattr(pr, "inner_radius"):
        for r in (float(pr.inner_radius), float(pr.outer_radius)):
            lines.append(_ellipse_xy(pr.center.x, pr.center.y, r, r, 0.0))
    elif name.startswith("Point") and hasattr(pr, "center"):
        points.append((float(pr.center.x), float(pr.center.y)))
    # Unknown region types are silently skipped.
    return lines, points


def load_ds9_regions(path: str, wcs: WCS, factor: int = 1) -> dict:
    """Parse a DS9 ``.reg`` and project every shape into image-pixel coords.

    Sky regions are converted with ``to_pixel(wcs)`` (so a downsample baked
    into ``wcs`` is handled automatically); native pixel regions are scaled by
    ``1/factor`` to match a downsampled display.
    """
    from regions import PixelRegion, Regions

    regs = Regions.read(path, format="ds9")
    all_lines: list = []
    all_points: list = []
    for r in regs:
        if isinstance(r, PixelRegion):
            pr, div = r, float(factor)
        else:
            pr, div = r.to_pixel(wcs), 1.0
        lines, points = _pixel_region_paths(pr)
        if div != 1.0:
            lines = [([x / div for x in xs], [y / div for y in ys])
                     for xs, ys in lines]
            points = [(x / div, y / div) for x, y in points]
        all_lines.extend(lines)
        all_points.extend(points)
    return {"lines": all_lines, "points": all_points, "n": len(regs)}


def _parse_contour_segments(path: str) -> list:
    """Split a DS9 contour file into per-segment vertex lists.

    Handles both layouts:

    - **DS9 ``.ctr`` (v7.5+)** — a frame line + ``level=N`` markers, with each
      contour wrapped in ``(`` … ``)``. The parentheses delimit segments.
    - **Plain ``.con``** — whitespace ``x y`` per line, blank lines separating
      segments.

    Non-coordinate lines (``#`` comments, ``global …``, the coordinate-frame
    name, ``level=…``) are skipped. ``(`` / ``)`` and blank lines flush the
    current segment.
    """
    segments: list = []
    cur: list = []

    def _flush():
        nonlocal cur
        if cur:
            segments.append(cur)
            cur = []

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                _flush()
                continue
            if line[0] == "(":            # DS9 .ctr opens a contour
                _flush()
                line = line[1:].strip()
                if not line:
                    continue
            closed = False
            if line.endswith(")"):        # …and closes one
                line = line[:-1].strip()
                closed = True
            if line.startswith(("#", "//")):
                _flush()
                continue
            parts = line.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    cur.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    # Metadata token (frame name, 'level=1', 'global …') —
                    # skip it without breaking the current segment.
                    pass
            if closed:
                _flush()
    _flush()
    return segments


def load_ds9_contours(
    path: str,
    wcs: WCS,
    coordsys: str = "sky",
    factor: int = 1,
) -> dict:
    """Parse a DS9 contour file (``.ctr`` / ``.con``) into image-pixel polylines.

    ``coordsys='sky'`` interprets the vertices as RA/Dec degrees and projects
    them through ``wcs``; ``coordsys='image'`` treats them as original-image
    pixels and divides by ``factor`` for a downsampled display. Segment layout
    (parenthesised ``.ctr`` vs blank-line ``.con``) is auto-detected — see
    :func:`_parse_contour_segments`.
    """
    segments = [s for s in _parse_contour_segments(path) if len(s) >= 2]
    lines: list = []
    if str(coordsys).lower().startswith("sky"):
        # Batch every vertex into ONE WCS transform — a contour file can hold
        # hundreds of segments, and per-segment transforms would be slow.
        flat = [pt for seg in segments for pt in seg]
        if flat:
            a = np.asarray(flat, dtype=float)
            sky = SkyCoord(a[:, 0], a[:, 1], unit="deg", frame="icrs")
            xall, yall = skycoord_to_pixel(sky, wcs, origin=0)
            xall = np.asarray(xall)
            yall = np.asarray(yall)
            i = 0
            for seg in segments:
                j = i + len(seg)
                lines.append((xall[i:j].tolist(), yall[i:j].tolist()))
                i = j
    else:  # image pixels (original resolution) → divide by downsample
        for seg in segments:
            a = np.asarray(seg, dtype=float)
            lines.append(((a[:, 0] / factor).tolist(),
                          (a[:, 1] / factor).tolist()))
    return {"lines": lines, "points": [], "n": len(lines)}
