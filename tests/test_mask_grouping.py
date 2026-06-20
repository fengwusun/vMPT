"""Regression tests for two spec-overlap bugs found with a real GO-6796 mask
under G140H/F070LP:

1. **Non-contiguous opens in one column were merged into one bogus slitlet.**
   Anonymous opens (mask CSV, no target id) at the same (q, d) but in
   physically separate row clusters (e.g. Q4 d=303 s=22-24 AND s=108-110) were
   grouped together; the group's (min,max) row span then masked the whole empty
   gap between the clusters, and closing a single shutter never shrank it.
   Fixed by splitting each column bucket into maximal runs of consecutive rows.

2. **Cross-quadrant conflicts promoted each slitlet's ENTIRE band to purple.**
   `window=None` (unbounded) made every shutter a cross-quad-conflicted slitlet
   contaminated — in both quadrants — a Mask Conflict (tens of thousands for a
   full mask). Fixed by bounding cross-quad purple to the genuine overlap:
   shutters contaminated by BOTH a slitlet and its cross-quad partner.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402

EX_JPG = Path.home() / ".vmpt/examples/example_r0600/JWST_F090W_F200W_F444W.jpg"
EX_WCS = Path.home() / ".vmpt/examples/example_r0600/wcs.fits"

pytestmark = pytest.mark.skipif(
    not (EX_JPG.exists() and EX_WCS.exists()),
    reason="example_r0600 assets missing",
)


def _setup(monkeypatch):
    from vmpt.image_io import load_jpg_with_sidecar
    # No stuck-opens, so only the constructed user opens drive the overlay.
    monkeypatch.setattr(
        m, "_FLAT_REASON", np.where(m._FLAT_REASON == 2, 0, m._FLAT_REASON))
    img = load_jpg_with_sidecar(str(EX_JPG), str(EX_WCS), max_dim=4000)
    H, W = img.shape[:2]
    centre = img.wcs.pixel_to_world(W / 2, H / 2)
    m.state["image"] = img
    m.state["pa_v3"] = 137.0
    m.ra_input.value = repr(float(centre.ra.deg))
    m.dec_input.value = repr(float(centre.dec.deg))
    m.fig.x_range.start, m.fig.x_range.end = -50000.0, 50000.0
    m.fig.y_range.start, m.fig.y_range.end = -50000.0, 50000.0
    m.state["disperser"] = "G140H"
    m.state["filter"] = "F070LP"


def _render(opens):
    m.state["open_shutters"] = {(s.q, s.s, s.d): s for s in opens}
    m.state.pop("_spec_overlap_cache", None)
    m.refresh_overlays()


def _masked_rows(quad):
    rows = set()
    for cds in (m.src_spec_overlap_both, m.src_spec_overlap_user,
                m.src_spec_overlap_stuck):
        for q_, s_ in zip(cds.data["q"], cds.data["s"]):
            if int(q_) == quad:
                rows.add(int(s_))
    return rows


def test_noncontiguous_column_opens_do_not_mask_the_gap(monkeypatch):
    """Two separate slitlets in the same column (Q4 d=303 s=22-24 and
    s=108-110) must mask only their own ±1 rows — NOT the empty gap between
    them. (The old single merged group masked all of rows 22-110.)"""
    _setup(monkeypatch)
    opens = [m.OpenShutter(q=4, s=s, d=303)
             for s in (22, 23, 24, 108, 109, 110)]
    _render(opens)
    rows = _masked_rows(4)
    # The two clusters' ±1 bands are masked …
    assert {22, 23, 24, 25, 21} & rows
    assert {107, 108, 109, 110, 111} & rows
    # … but the whole gap between them is clean.
    gap = sorted(r for r in rows if 26 <= r <= 106)
    assert gap == [], f"gap rows falsely masked: {gap}"


def test_closing_one_cluster_updates_status(monkeypatch):
    """Closing the upper cluster must clear its masking — the core symptom
    (with the merged group, the span never shrank so nothing changed)."""
    _setup(monkeypatch)
    full = [m.OpenShutter(q=4, s=s, d=303)
            for s in (22, 23, 24, 108, 109, 110)]
    _render(full)
    assert {107, 108, 109, 110, 111} & _masked_rows(4)
    # Close s108-110 → its band (rows 107-111) must go clean; lower stays.
    _render([m.OpenShutter(q=4, s=s, d=303) for s in (22, 23, 24)])
    rows = _masked_rows(4)
    assert not ({107, 108, 109, 110, 111} & rows), "upper cluster band not cleared"
    assert {21, 22, 23, 24, 25} & rows, "lower cluster should remain"


def test_cross_quadrant_conflict_is_bounded_not_whole_band(monkeypatch):
    """A Q2↔Q4 cross-quadrant conflict must paint purple only on the genuine
    doubly-hit overlap, leaving singly-hit shutters orange — NOT promote both
    slitlets' entire bands to purple (the old window=None behaviour, which
    left orange=0)."""
    _setup(monkeypatch)
    opens = ([m.OpenShutter(q=2, s=s, d=80) for s in (79, 80, 81)]
             + [m.OpenShutter(q=4, s=s, d=300) for s in (79, 80, 81)])
    _render(opens)
    purple = len(m.src_spec_overlap_both.data["q"])
    orange = len(m.src_spec_overlap_user.data["q"])
    assert purple > 0, "the cross-quad pair should still register a conflict"
    assert orange > 0, "singly-hit shutters must stay orange (purple bounded)"
    # Closing the Q4 slitlet removes the conflict → no purple left.
    _render([m.OpenShutter(q=2, s=s, d=80) for s in (79, 80, 81)])
    assert len(m.src_spec_overlap_both.data["q"]) == 0


def test_cross_quad_conflicted_pick_itself_is_purple(monkeypatch):
    """A picked shutter belonging to a slitlet in a cross-quadrant conflict
    must itself read as Mask Conflict (purple) — its own spectrum is one of
    the two colliding spectra, even though the partner's light lands on
    neighbouring shutters, not this exact one. Regression for the reported
    Q3 s141 d117 pick (conflicts with Q1 d107) that rendered orange/clean."""
    _setup(monkeypatch)
    opens = ([m.OpenShutter(q=2, s=s, d=80) for s in (79, 80, 81)]
             + [m.OpenShutter(q=4, s=s, d=300) for s in (79, 80, 81)])
    _render(opens)
    purple = {(int(q), int(s), int(d)) for q, s, d in zip(
        m.src_spec_overlap_both.data["q"], m.src_spec_overlap_both.data["s"],
        m.src_spec_overlap_both.data["d"])}
    # The picks themselves are purple (they ARE the conflicting slitlets).
    assert (2, 80, 80) in purple
    assert (4, 80, 300) in purple
    # Close the Q4 slitlet → the Q2 picks are no longer in a conflict.
    _render([m.OpenShutter(q=2, s=s, d=80) for s in (79, 80, 81)])
    purple2 = {(int(q), int(s), int(d)) for q, s, d in zip(
        m.src_spec_overlap_both.data["q"], m.src_spec_overlap_both.data["s"],
        m.src_spec_overlap_both.data["d"])}
    assert (2, 80, 80) not in purple2
