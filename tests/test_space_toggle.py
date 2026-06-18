"""Hover + Space toggles exactly ONE shutter, independent of the N-shutter
slitlet setting (`_toggle_single_shutter_at`).

  - open a closed operable shutter → exactly one shutter opens (even with
    slitlet height 3/5);
  - Space on an open shutter closes only that one, leaving slitlet siblings;
  - a non-operable closed shutter is refused.
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


def _operable_with_neighbors():
    """A (q,s,d) (1-based) whose s-1, s, s+1 rows are all operable, with s in
    a safe interior range."""
    op = m.OPERABLE
    for q in range(op.shape[0]):
        for d in range(op.shape[2]):
            col = op[q, :, d]
            for s in range(2, op.shape[1] - 2):
                if col[s - 1] and col[s] and col[s + 1]:
                    return q + 1, s + 1, d + 1
    return None


def test_space_opens_exactly_one_shutter():
    qsd = _operable_with_neighbors()
    assert qsd is not None
    q, s, d = qsd
    m._set_open_shutters({})
    m.state["slitlet_height"] = 5  # would open 5 via a click; Space opens 1
    m._toggle_single_shutter_at(q, s, d)
    assert (q, s, d) in m.state["open_shutters"]
    assert len(m.state["open_shutters"]) == 1, \
        "Space must open exactly one shutter, not an N-shutter slitlet"
    # toggle again → closes it
    m._toggle_single_shutter_at(q, s, d)
    assert (q, s, d) not in m.state["open_shutters"]
    assert len(m.state["open_shutters"]) == 0


def test_space_closes_only_the_hovered_shutter():
    qsd = _operable_with_neighbors()
    q, s, d = qsd
    # Hand-build a 3-shutter slitlet, then Space-close only the centre.
    m._set_open_shutters({
        (q, s - 1, d): m.OpenShutter(q=q, s=s - 1, d=d, role="sky"),
        (q, s, d): m.OpenShutter(q=q, s=s, d=d, role="target"),
        (q, s + 1, d): m.OpenShutter(q=q, s=s + 1, d=d, role="sky"),
    })
    m._toggle_single_shutter_at(q, s, d)
    assert (q, s, d) not in m.state["open_shutters"]
    assert (q, s - 1, d) in m.state["open_shutters"]
    assert (q, s + 1, d) in m.state["open_shutters"]
    assert len(m.state["open_shutters"]) == 2, "siblings must stay open"


def test_space_refuses_non_operable_shutter():
    bad = np.argwhere(~m.OPERABLE)
    assert len(bad), "expected some non-operable shutters"
    q, s, d = (int(bad[0][0]) + 1, int(bad[0][1]) + 1, int(bad[0][2]) + 1)
    m._set_open_shutters({})
    m._toggle_single_shutter_at(q, s, d)
    assert (q, s, d) not in m.state["open_shutters"], \
        "a closed non-operable shutter must not open"
