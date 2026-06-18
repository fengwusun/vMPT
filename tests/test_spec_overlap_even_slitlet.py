"""Regression for the even-shutter contamination-band skew.

The spec-overlap band used to be modelled as ``s_center ± half_extent``,
symmetric around a *rounded* centre. For an even-shutter slitlet that skews
the band: a 2-shutter slitlet at s=[112,113] rounds its centre to 112 and
reaches one row LOWER than the 3-shutter slitlet [112,113,114] that
contains it — so the 2-shutter pick flagged a "Mask Conflict" (purple) that
its own 3-shutter superset did not. Anchoring the band on the slitlet's
actual rows [s_lo, s_hi] makes it monotonic in the open set.

Reproduces the user's exact case (PRISM/CLEAR): slitlet A on Q1 d225
(s108-110) plus a slitlet B on Q1 d217.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

EX_JPG = Path.home() / ".vmpt/examples/example_r0600/JWST_F090W_F200W_F444W.jpg"
EX_WCS = Path.home() / ".vmpt/examples/example_r0600/wcs.fits"

pytestmark = pytest.mark.skipif(
    not (EX_JPG.exists() and EX_WCS.exists()),
    reason="example_r0600 assets missing",
)


@pytest.fixture()
def app():
    import vmpt.main as m
    from vmpt.image_io import load_jpg_with_sidecar
    img = load_jpg_with_sidecar(str(EX_JPG), str(EX_WCS), max_dim=4000)
    H, W = img.shape[:2]
    c = img.wcs.pixel_to_world(W / 2, H / 2)
    m.state["image"] = img
    m.ra_input.value = repr(float(c.ra.deg))
    m.dec_input.value = repr(float(c.dec.deg))
    m.state["ra_deg"] = float(c.ra.deg)
    m.state["dec_deg"] = float(c.dec.deg)
    m.state["pa_v3"] = 137.0
    m.state["disperser"] = "PRISM"
    m.state["filter"] = "CLEAR"
    # Whole MSA in view so every purple shutter renders into the CDS.
    m.fig.x_range.start, m.fig.x_range.end = -50000.0, 50000.0
    m.fig.y_range.start, m.fig.y_range.end = -50000.0, 50000.0
    return m


def _purple_after(m, opens):
    m.state.pop("_spec_overlap_cache", None)
    m._set_open_shutters(opens)
    m.refresh_overlays()
    d = m.src_spec_overlap_both.data
    return {(int(q), int(s), int(dd))
            for q, s, dd in zip(d.get("q", []), d.get("s", []), d.get("d", []))}


def _os(m, q, s, d, role):
    return m.OpenShutter(q=q, s=s, d=d, role=role, target_id=None)


def test_even_slitlet_no_phantom_conflict(app):
    m = app
    A = {(1, 108, 225): _os(m, 1, 108, 225, "sky"),
         (1, 109, 225): _os(m, 1, 109, 225, "target"),
         (1, 110, 225): _os(m, 1, 110, 225, "sky")}
    B2 = {(1, 112, 217): _os(m, 1, 112, 217, "sky"),
          (1, 113, 217): _os(m, 1, 113, 217, "target")}
    B3 = {**B2, (1, 114, 217): _os(m, 1, 114, 217, "sky")}

    p2 = _purple_after(m, {**A, **B2})
    p3 = _purple_after(m, {**A, **B3})

    # Monotonicity: the 2-shutter slitlet's conflicts must be a subset of the
    # 3-shutter superset's — never the other way round.
    assert p2 <= p3, (
        f"2-shutter pick flagged a conflict its 3-shutter superset doesn't: "
        f"extra = {sorted(p2 - p3)[:10]}"
    )
    # In this geometry the slitlets are >1 row apart → no real collision.
    assert not p3, f"unexpected purple for the 3-shutter case: {sorted(p3)[:10]}"
    assert not p2, f"unexpected purple for the 2-shutter case: {sorted(p2)[:10]}"


def test_even_slitlet_orange_band_stops_at_one_buffer_row(app):
    """A 5-shutter slitlet centred on Q1 s109 d225 with s107 inoperable opens
    s108-111 (4 shutters = even). The orange "Masked" band must reach exactly
    ONE buffer row past the top open shutter (s111 → s112), not two (s113).
    Was broken: the rounded centre + half_extent gave a [107,113] band."""
    m = app
    assert not bool(m.OPERABLE[0, 107 - 1, 225 - 1]), \
        "test premise: Q1 s107 d225 should be inoperable"
    opens = {(1, s, 225): _os(m, 1, s, 225, "target" if s == 109 else "sky")
             for s in (108, 109, 110, 111)}
    m.state.pop("_spec_overlap_cache", None)
    m._set_open_shutters(opens)
    m.refresh_overlays()
    rows = {int(s) for s in m.src_spec_overlap_user.data.get("s", [])}
    assert 112 in rows, "expected orange on the +1 buffer row s112"
    assert 113 not in rows, \
        "orange must not reach s113 — the even-slitlet band overshot by a row"
