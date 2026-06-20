"""Grating diagonal-step relaxation of the Mask-Conflict (purple) rule.

A no-buffer ADJACENCY conflict (two slitlets exactly 1 row apart, no real row
overlap) is purple under PRISM (matches APT/MPT) but, under a grating, is
demoted to orange (Masked) once the slitlets are far enough apart in COLUMNS —
long grating spectra run parallel, so a column-offset diagonal step isn't a
real collision. Threshold scales with spectral length: H = 20 cols, M = 8.

Real row overlap and PRISM are never relaxed.
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
from vmpt.wavelengths import grating_adjacency_min_colsep  # noqa: E402


# ── pure-logic unit tests (no image / network) ────────────────────────────

def test_min_colsep_classifier():
    # Threshold is 1 for any M/H grating: any nonzero column offset is a
    # deliberate diagonal step (magnitude is physically irrelevant — adjacent
    # rows stay a fixed 1 row apart in cross-dispersion regardless of column).
    assert grating_adjacency_min_colsep("G140H") == 1
    assert grating_adjacency_min_colsep("G235H") == 1
    assert grating_adjacency_min_colsep("G395H") == 1
    assert grating_adjacency_min_colsep("G140M") == 1
    assert grating_adjacency_min_colsep("G235M") == 1
    assert grating_adjacency_min_colsep("PRISM") is None      # ends in M, but not a grating
    assert grating_adjacency_min_colsep("") is None
    assert grating_adjacency_min_colsep(None) is None


def test_relaxed_for_any_grating_diagonal_step():
    H = grating_adjacency_min_colsep("G140H")    # 1
    Mm = grating_adjacency_min_colsep("G235M")   # 1
    # lower slitlet rows 40-42, upper 43-45 → adjacency (gap 0).
    # Any nonzero column offset → relaxed (deliberate diagonal step).
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 102, H) is True   # Δd=2 (user case)
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 101, H) is True   # Δd=1
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 100, H) is False  # Δd=0 same column
    # M grating uses the same threshold.
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 102, Mm) is True   # Δd=2
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 100, Mm) is False  # Δd=0
    # symmetric in slitlet order
    assert m._grating_adjacency_relaxed(43, 45, 102, 40, 42, 100, H) is True


def test_never_relaxed_for_prism_or_real_overlap():
    H = grating_adjacency_min_colsep("G140H")
    # PRISM (min_colsep None) → never relaxed, even far apart.
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 200, None) is False
    # Real row overlap (rows intersect) → always purple, even far apart.
    assert m._grating_adjacency_relaxed(40, 42, 100, 42, 44, 200, H) is False
    # A buffered gap (≥1 row between, gap≠0) is not an adjacency → not relaxed
    # (and wouldn't be a conflict anyway).
    assert m._grating_adjacency_relaxed(40, 42, 100, 44, 46, 200, H) is False
    # Same column adjacency stays purple (Δd=0 < threshold).
    assert m._grating_adjacency_relaxed(40, 42, 100, 43, 45, 100, H) is False


# ── end-to-end through the real refresh_overlays ───────────────────────────

EX_JPG = Path.home() / ".vmpt/examples/example_r0600/JWST_F090W_F200W_F444W.jpg"
EX_WCS = Path.home() / ".vmpt/examples/example_r0600/wcs.fits"

pytestmark_e2e = pytest.mark.skipif(
    not (EX_JPG.exists() and EX_WCS.exists()),
    reason="example_r0600 assets missing",
)


@pytestmark_e2e
def test_relaxation_flips_rendered_purple(monkeypatch):
    """Open two adjacent N=3 slitlets a given column offset apart and read the
    real purple/orange CDS. Stuck-opens are disabled so the A-B user pair is
    the only conflict, isolating the relaxation."""
    from vmpt.image_io import load_jpg_with_sidecar

    # Disable stuck-opens (REASON==2 → 0) so only the user A-B pair conflicts.
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
    q, R, dA = 1, 40, 40

    def purple_count(disp, filt, sep):
        opens = {}
        for i in range(3):
            opens[(q, R + i, dA)] = m.OpenShutter(q=q, s=R + i, d=dA, target_id="A")
        for i in range(3):
            opens[(q, R + 3 + i, dA + sep)] = m.OpenShutter(
                q=q, s=R + 3 + i, d=dA + sep, target_id="B")
        m.state["disperser"] = disp
        m.state["filter"] = filt
        m.state["open_shutters"] = opens
        m.state.pop("_spec_overlap_cache", None)
        m.refresh_overlays()
        picks = set(opens)
        return sum(
            1 for q_, s_, d_ in zip(m.src_spec_overlap_both.data["q"],
                                    m.src_spec_overlap_both.data["s"],
                                    m.src_spec_overlap_both.data["d"])
            if (int(q_), int(s_), int(d_)) not in picks)

    # Grating: any nonzero column offset is a deliberate diagonal step → no
    # purple. (Δd=0 same-column isn't reachable here — vMPT's `different_col`
    # filter never flags two slitlets in the exact same column as colliding;
    # that case is covered by the pure-logic test above.)
    assert purple_count("G140H", "F070LP", 2) == 0    # Δd=2 (user's diagonal step)
    assert purple_count("G140H", "F070LP", 1) == 0    # Δd=1
    assert purple_count("G235M", "F170LP", 2) == 0    # M uses the same threshold
    # PRISM is never relaxed — the same adjacency stays purple.
    assert purple_count("PRISM", "CLEAR", 2) > 0
    assert purple_count("PRISM", "CLEAR", 1) > 0


@pytestmark_e2e
def test_user_diagonal_step_q2_g140h_not_purple(monkeypatch):
    """User's exact report: Q2 s53-55 d315 + Q2 s56-58 d313 (adjacent slitlets,
    Δd=2 columns) under G140H/F100LP must NOT be a Mask Conflict — it's a
    deliberate 2-column diagonal step. Under PRISM the same picks stay purple."""
    from vmpt.image_io import load_jpg_with_sidecar

    # Disable stuck-opens so only the A-B pair can drive purple.
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

    def purple_rows(disp, filt):
        opens = {}
        for s in (53, 54, 55):
            opens[(2, s, 315)] = m.OpenShutter(q=2, s=s, d=315, target_id="A")
        for s in (56, 57, 58):
            opens[(2, s, 313)] = m.OpenShutter(q=2, s=s, d=313, target_id="B")
        m.state["disperser"] = disp
        m.state["filter"] = filt
        m.state["open_shutters"] = opens
        m.state.pop("_spec_overlap_cache", None)
        m.refresh_overlays()
        return {(int(q_), int(s_))
                for q_, s_ in zip(m.src_spec_overlap_both.data["q"],
                                  m.src_spec_overlap_both.data["s"])}

    # G140H: the 2-column diagonal step relaxes → no purple anywhere.
    assert purple_rows("G140H", "F100LP") == set()
    # PRISM: same picks are a genuine no-buffer conflict → purple appears
    # (bounded to the central ±2-row window of the two slitlets).
    prism = purple_rows("PRISM", "CLEAR")
    assert prism, "PRISM must still flag the adjacency as a Mask Conflict"
    assert {s for (q_, s) in prism} <= {54, 55, 56, 57}


# ── optimizer collision-protection side ────────────────────────────────────

def _evaluator(disperser, filt):
    from vmpt.optimizer import PointingEvaluator
    # ra/dec are unused by _apply_collision_drops (it takes quad/s/d directly);
    # source 0 protected, source 1 not. reason=all-zeros disables stuck-opens
    # so rule 1 can't interfere — we test the protected↔unprotected rule.
    ra = np.array([53.16, 53.16])
    dec = np.array([-27.78, -27.78])
    return PointingEvaluator(
        ra, dec, centration="UNCONSTRAINED", slit_length=3,
        protect_mask=np.array([True, False]),
        priorities=np.array([1.0, 2.0]), weights=np.array([1.0, 1.0]),
        disperser=disperser, filt=filt,
        reason=np.zeros((4, 171, 365), dtype=np.int8))


def test_optimizer_row_collide_relaxation():
    ev = _evaluator("G140H", "F070LP")
    assert ev._adj_colsep == 1
    T, F = np.array(True), np.array(False)
    # adjacency (|Δrow|=1), same quad, any nonzero col offset → relaxed
    assert bool(ev._row_collide(np.array(1), np.array(2), T)) is False   # Δd=2 diagonal step
    assert bool(ev._row_collide(np.array(1), np.array(1), T)) is False   # Δd=1
    assert bool(ev._row_collide(np.array(1), np.array(0), T)) is True    # Δd=0 same column
    assert bool(ev._row_collide(np.array(0), np.array(50), T)) is True   # real row overlap
    assert bool(ev._row_collide(np.array(1), np.array(50), F)) is True   # cross-quadrant
    evp = _evaluator("PRISM", "CLEAR")
    assert evp._adj_colsep is None
    assert bool(evp._row_collide(np.array(1), np.array(50), T)) is True  # PRISM never relaxed


def test_optimizer_drops_relax_for_grating_diagonal_step():
    """Rule-3 integration: a protected source at (q1, s100, d100) and an
    unprotected one at (q1, s103, d100+Δd) are row-adjacent (slitlets [99,101]
    vs [102,104], 1 row apart). Under a grating, any nonzero column offset is a
    deliberate diagonal step so the unprotected source is NOT dropped; in the
    SAME column (Δd=0), or under PRISM, it IS dropped."""
    detected = np.array([True, True])
    quad = np.array([1, 1])
    s_frac = np.array([100.0, 103.0])

    def kept_unprotected(disperser, filt, dcol):
        ev = _evaluator(disperser, filt)
        d_frac = np.array([100.0, 100.0 + dcol])
        return bool(ev._apply_collision_drops(detected, quad, s_frac, d_frac)[1])

    assert kept_unprotected("G140H", "F070LP", 2) is True     # Δd=2 diagonal step → kept
    assert kept_unprotected("G140H", "F070LP", 1) is True     # Δd=1 → kept
    assert kept_unprotected("G140H", "F070LP", 0) is False    # Δd=0 same column → dropped
    assert kept_unprotected("G235M", "F170LP", 2) is True     # M same threshold → kept
    assert kept_unprotected("G235M", "F170LP", 0) is False    # Δd=0 same column → dropped
    assert kept_unprotected("PRISM", "CLEAR", 2) is False     # never relaxed → dropped
