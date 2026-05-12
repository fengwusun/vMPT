"""Tests for app/session_io.py."""

from __future__ import annotations

import json

import pytest

from app.empt_io import OpenShutter
from app.session_io import Session, export_session_json, import_session_json


def _sample_session() -> Session:
    return Session(
        pointing_ra_deg=189.12,
        pointing_dec_deg=62.21,
        pa_v3_deg=273.0,
        disperser="G395M",
        filter_name="F290LP",
        slitlet_height=3,
        open_shutters=[
            OpenShutter(q=2, d=200, s=86, target_id="123456", role="target"),
            OpenShutter(q=2, d=200, s=85, target_id="123456", role="sky"),
            OpenShutter(q=2, d=200, s=87, target_id="123456", role="sky"),
            OpenShutter(q=3, d=100, s=42, target_id=None, role="manual"),
        ],
        highlighted=[(2, 80, 200), (1, 10, 50)],
        image_path="/tmp/loaded.fits",
        catalog_path="/tmp/catalog.csv",
    )


def _compare(a: Session, b: Session) -> None:
    assert a.pointing_ra_deg == pytest.approx(b.pointing_ra_deg)
    assert a.pointing_dec_deg == pytest.approx(b.pointing_dec_deg)
    assert a.pa_v3_deg == pytest.approx(b.pa_v3_deg)
    assert a.disperser == b.disperser
    assert a.filter_name == b.filter_name
    assert a.slitlet_height == b.slitlet_height
    assert len(a.open_shutters) == len(b.open_shutters)
    for sa, sb in zip(a.open_shutters, b.open_shutters):
        assert (sa.q, sa.s, sa.d) == (sb.q, sb.s, sb.d)
        assert sa.target_id == sb.target_id
        assert sa.role == sb.role
    assert a.highlighted == b.highlighted
    assert a.image_path == b.image_path
    assert a.catalog_path == b.catalog_path


def test_roundtrip_full(tmp_path):
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    loaded = import_session_json(str(p))
    _compare(s, loaded)
    # manual shutter target_id is preserved as None
    manual = [sh for sh in loaded.open_shutters if sh.role == "manual"][0]
    assert manual.target_id is None


def test_roundtrip_empty(tmp_path):
    s = Session(
        pointing_ra_deg=0.0,
        pointing_dec_deg=0.0,
        pa_v3_deg=0.0,
        disperser="PRISM",
        filter_name="CLEAR",
        slitlet_height=3,
        open_shutters=[],
        highlighted=[],
    )
    p = tmp_path / "empty.json"
    export_session_json(s, str(p))
    loaded = import_session_json(str(p))
    _compare(s, loaded)
    assert loaded.open_shutters == []
    assert loaded.highlighted == []
    assert loaded.image_path is None
    assert loaded.catalog_path is None


def test_missing_pointing_raises(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({
        "instrument": {"disperser": "G395M", "filter": "F290LP"},
        "open_shutters": [],
    }))
    with pytest.raises(ValueError):
        import_session_json(str(p))


def test_unknown_keys_tolerated(tmp_path):
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    # Inject unknown keys at the top level and inside nested objects
    data = json.loads(p.read_text())
    data["future_field"] = {"anything": 42}
    data["pointing"]["wishlist"] = "more precision"
    data["instrument"]["mode"] = "experimental"
    data["open_shutters"][0]["future_attr"] = "ok"
    p.write_text(json.dumps(data, indent=2))
    loaded = import_session_json(str(p))
    _compare(s, loaded)
