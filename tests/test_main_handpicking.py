"""Unit tests for the hand-picking helpers in app/main.py.

Importing vmpt.main attaches widgets to curdoc() — harmless during tests.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from vmpt.coords import v2v3_to_radec  # noqa: E402
from vmpt.msa import load_msa_grid  # noqa: E402


def _import_main():
    """Lazy import so any import-time side effects don't pollute other tests."""
    import vmpt.main as m
    return m


def test_sky_to_v2v3_round_trip():
    m = _import_main()
    v2_msa, v3_msa = load_msa_grid()
    fiducial = SkyCoord(53.14 * u.deg, -27.79 * u.deg)
    pa_v3 = 321.0
    # Pick a shutter, forward-transform to sky, then inverse-transform back to V2/V3
    v2c = float(v2_msa[1, 86, 200])
    v3c = float(v3_msa[1, 86, 200])
    sky = v2v3_to_radec(fiducial, pa_v3, np.array([[v2c, v3c]]))
    sky_one = sky[0]
    v2_back, v3_back = m._sky_to_v2v3(sky_one, fiducial, pa_v3)
    assert abs(v2_back - v2c) < 1e-3, f"v2 round-trip drift {abs(v2_back - v2c)}"
    assert abs(v3_back - v3c) < 1e-3, f"v3 round-trip drift {abs(v3_back - v3c)}"


def test_nearest_shutter_finds_known_position():
    m = _import_main()
    v2_msa, v3_msa = load_msa_grid()
    # Target the center of a known operable shutter
    q, s, d = 2, 86, 200
    v2t = float(v2_msa[q - 1, s - 1, d - 1])
    v3t = float(v3_msa[q - 1, s - 1, d - 1])
    result = m._nearest_shutter(v2t, v3t, require_operable=False)
    assert result == (q, s, d)


def test_add_slitlet_builds_three_shutters():
    m = _import_main()
    m.state["open_shutters"] = {}
    m.state["slitlet_height"] = 3
    # Use a Q2 column where many shutters are operable
    q, s_center, d = 2, 86, 200
    added = m._add_slitlet(q, s_center, d, target_id="42")
    assert added >= 1
    # The center shutter (if operable) should be present
    if m.OPERABLE[q - 1, s_center - 1, d - 1]:
        assert (q, s_center, d) in m.state["open_shutters"]
        assert m.state["open_shutters"][(q, s_center, d)].role == "target"
    # Cleanup
    m.state["open_shutters"] = {}


def test_history_undo():
    m = _import_main()
    m.state["open_shutters"] = {}
    m.state["history"] = []
    m._push_history()
    m.state["open_shutters"][(2, 86, 200)] = m.OpenShutter(q=2, s=86, d=200)
    assert len(m.state["history"]) == 1
    m.state["open_shutters"] = m.state["history"].pop()
    assert m.state["open_shutters"] == {}
