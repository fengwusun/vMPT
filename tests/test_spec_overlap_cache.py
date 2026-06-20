"""The spec-overlap memoization in ``vmpt.main.refresh_overlays`` must:

  (1) produce IDENTICAL overlay CDS on a repeat call (cache HIT) as on the
      first call (cache MISS) — proving the cache is a pure optimization;
  (2) invalidate when any input that affects the colours changes
      (open-shutter mask, disperser, filter, or an alpha slider).

This guards the *view-independent* colour cache that fixed the
zoom-dependent purple/orange flicker: the conflict identification and
colour assignment depend only on the open mask + disperser/filter +
alpha sliders, never on the pointing or the zoom, so refresh_overlays
re-runs on every pan/zoom but only re-culls the cached result.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

EX_R0600_JPG = _ROOT / "example_r0600" / "JWST_F090W_F200W_F444W.jpg"
EX_R0600_WCS = _ROOT / "example_r0600" / "wcs.fits"

pytestmark = pytest.mark.skipif(
    not (EX_R0600_JPG.exists() and EX_R0600_WCS.exists()),
    reason="example_r0600 missing",
)


def _snapshot(m):
    """Frozen copy of the three spec-overlap CDS payloads."""
    return {
        name: {k: list(v) for k, v in src.data.items()}
        for name, src in (
            ("stuck", m.src_spec_overlap_stuck),
            ("user", m.src_spec_overlap_user),
            ("both", m.src_spec_overlap_both),
        )
    }


@pytest.fixture()
def app_with_collision():
    """Load vmpt.main, drop in the r0600 image + a known colliding pair
    of Q1 slitlets, and make the whole MSA 'in view' so the coloured
    shutters actually populate the CDS."""
    import vmpt.main as m
    from vmpt.image_io import load_jpg_with_sidecar

    img = load_jpg_with_sidecar(str(EX_R0600_JPG), str(EX_R0600_WCS), max_dim=4000)
    H, W = img.shape[:2]
    center = img.wcs.pixel_to_world(W / 2, H / 2)

    m.state["image"] = img
    m.ra_input.value = repr(float(center.ra.deg))
    m.dec_input.value = repr(float(center.dec.deg))
    m.state["pa_v3"] = 137.0
    m.state["disperser"] = "PRISM"
    m.state["filter"] = "CLEAR"
    # Two touching Q1 slitlets — the canonical purple case (see
    # test_spec_overlap_purple.py): A at (s=100-102, d=322) borders B at
    # (s=97-99, d=323) with no operable row between, so both inner
    # shutters go purple and their bands paint orange.
    m.state["open_shutters"] = {
        (1, 100, 322): m.OpenShutter(q=1, s=100, d=322, target_id="A"),
        (1, 101, 322): m.OpenShutter(q=1, s=101, d=322, target_id="A"),
        (1, 102, 322): m.OpenShutter(q=1, s=102, d=322, target_id="A"),
        (1, 97, 323): m.OpenShutter(q=1, s=97, d=323, target_id="B"),
        (1, 98, 323): m.OpenShutter(q=1, s=98, d=323, target_id="B"),
        (1, 99, 323): m.OpenShutter(q=1, s=99, d=323, target_id="B"),
    }
    # Force the whole MSA into the view so every coloured shutter renders
    # (the cache stores the global colours; the CDS is the in-view cull).
    m.fig.x_range.start, m.fig.x_range.end = -50000.0, 50000.0
    m.fig.y_range.start, m.fig.y_range.end = -50000.0, 50000.0
    # Start from a clean cache each test.
    m.state.pop("_spec_overlap_cache", None)
    return m


def test_cache_hit_matches_miss(app_with_collision):
    m = app_with_collision

    m.refresh_overlays()                       # cache MISS
    assert m.state.get("_spec_overlap_cache") is not None
    sig_miss = m.state["_spec_overlap_cache"][0]
    snap_miss = _snapshot(m)

    # Sanity: the collision actually produced purple ("both") polygons,
    # otherwise the equality check below would pass trivially.
    assert len(snap_miss["both"]["xs"]) > 0, "expected purple polygons in view"

    m.refresh_overlays()                       # cache HIT (identical state)
    assert m.state["_spec_overlap_cache"][0] == sig_miss
    snap_hit = _snapshot(m)

    assert snap_hit == snap_miss, "cache HIT produced different CDS than MISS"


def test_cache_invalidates_on_open_change(app_with_collision):
    m = app_with_collision
    m.refresh_overlays()
    sig_a = m.state["_spec_overlap_cache"][0]

    # Add an unrelated open shutter → mask changes → signature changes.
    m.state["open_shutters"][(2, 50, 100)] = m.OpenShutter(
        q=2, s=50, d=100, target_id="C")
    m.refresh_overlays()
    sig_b = m.state["_spec_overlap_cache"][0]
    assert sig_b != sig_a, "adding an open shutter must invalidate the cache"


def test_cache_invalidates_on_disperser_and_alpha(app_with_collision):
    m = app_with_collision
    m.refresh_overlays()
    sig0 = m.state["_spec_overlap_cache"][0]

    m.state["disperser"] = "G140M"
    m.state["filter"] = "F070LP"
    m.refresh_overlays()
    sig1 = m.state["_spec_overlap_cache"][0]
    assert sig1 != sig0, "changing disperser/filter must invalidate the cache"

    # An alpha-slider change must also bust the cache (colours scale by it).
    m.state["overlap_base_alpha_both"] = 0.5
    m.refresh_overlays()
    sig2 = m.state["_spec_overlap_cache"][0]
    assert sig2 != sig1, "changing a base-alpha must invalidate the cache"


def test_purple_band_bounded_end_to_end(app_with_collision):
    """End-to-end through the REAL refresh_overlays (not a replica): the
    canonical colliding Q1 pair (A 100-102/d322 over B 97-99/d323) must
    paint purple ONLY in the ±2-row conflict window [98..101]. The band
    tails revert to orange, and purple never leaks past the window in any
    quadrant. Guards the bounded-purple (MPT-faithful) rule end to end."""
    m = app_with_collision
    m.refresh_overlays()
    both, user = m.src_spec_overlap_both.data, m.src_spec_overlap_user.data
    picks = set(m.state["open_shutters"])

    def rows(data, q=None):
        return {int(s) for qq, s, d in zip(data["q"], data["s"], data["d"])
                if (int(qq), int(s), int(d)) not in picks
                and (q is None or int(qq) == q)}

    WINDOW = {98, 99, 100, 101}     # [max(100,97)-2 .. min(102,99)+2]
    # Purple is bounded to the window in EVERY quadrant (no full-band leak,
    # no same-detector other-quadrant escalation).
    assert rows(both) <= WINDOW, (
        f"purple leaked past the window: {sorted(rows(both))}")
    assert rows(both, q=1) == WINDOW, "the whole window should be purple in Q1"
    # The Q1 band tails are orange (Masked), disjoint from the window.
    assert rows(user, q=1).isdisjoint(WINDOW)
    assert {96, 97, 102, 103} <= rows(user, q=1)
    # The two touching picks themselves render purple.
    both_keys = set(zip(
        (int(q) for q in both["q"]),
        (int(s) for s in both["s"]),
        (int(d) for d in both["d"])))
    assert (1, 100, 322) in both_keys and (1, 99, 323) in both_keys
