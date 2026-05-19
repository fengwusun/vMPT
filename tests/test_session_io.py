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
    # JSON with no `configs` and no `open_shutters` → unrecognized session.
    p = tmp_path / "bad.json"
    p.write_text(json.dumps({"instrument": "JWST/NIRSpec"}))
    with pytest.raises(ValueError):
        import_session_json(str(p))


def test_unknown_keys_tolerated(tmp_path):
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    # Inject unknown keys at the top level and inside the MPT structure.
    data = json.loads(p.read_text())
    data["future_field"] = {"anything": 42}
    data["catalog"]["wishlist"] = "more metadata"
    data["configs"][0]["future_attr"] = "ok"
    p.write_text(json.dumps(data, indent=2))
    # Same for the workspace sidecar.
    side = tmp_path / "vmpt_workspace.json"
    sdata = json.loads(side.read_text())
    sdata["future_field_in_workspace"] = "ok"
    side.write_text(json.dumps(sdata, indent=2))
    loaded = import_session_json(str(p))
    _compare(s, loaded)


def test_session_json_has_no_filesystem_paths(tmp_path):
    """The MPT session.json must NOT contain any image / catalog / sidecar
    paths — APT can't read those, and surfacing them would leak local
    filesystem details into the shared plan file."""
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    text = p.read_text()
    for token in ("image_path", "wcs_sidecar_path", "catalog_path",
                  "/tmp/", s.image_path or "__no_image__",
                  s.catalog_path or "__no_catalog__"):
        if token in ("__no_image__", "__no_catalog__"):
            continue
        assert token not in text, f"session.json leaked {token!r}"


def test_workspace_sidecar_written_and_consumed(tmp_path):
    """Save writes a sibling `vmpt_workspace.json`; load merges it back in."""
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    sidecar = tmp_path / "vmpt_workspace.json"
    assert sidecar.exists(), "workspace sidecar not written"
    sdata = json.loads(sidecar.read_text())
    assert sdata["image_path"] == s.image_path
    assert sdata["catalog_path"] == s.catalog_path
    # Roles + target_ids must be preserved end-to-end via the sidecar.
    loaded = import_session_json(str(p))
    _compare(s, loaded)
    roles = sorted(sh.role for sh in loaded.open_shutters)
    assert roles == ["manual", "sky", "sky", "target"]


def test_mpt_only_load_without_sidecar(tmp_path):
    """If only the MPT session.json is present (e.g. shared without the
    workspace sidecar), load still succeeds: pointing/PA/disperser/
    slitlet structure come back; per-shutter target_id is recovered via
    primaryIds positional alignment, role is reconstructed from slitlet
    geometry."""
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    # Drop the sidecar to simulate the APT-shared scenario.
    (tmp_path / "vmpt_workspace.json").unlink()
    loaded = import_session_json(str(p))
    assert loaded.disperser == s.disperser
    assert loaded.filter_name == s.filter_name
    # Same number of open shutters back out
    assert len(loaded.open_shutters) == len(s.open_shutters)
    # Slitlet geometry preserved
    coords_in = sorted((sh.q, sh.s, sh.d) for sh in s.open_shutters)
    coords_out = sorted((sh.q, sh.s, sh.d) for sh in loaded.open_shutters)
    assert coords_in == coords_out
    # Targeted shutters' target_id recovered (manual ones lose it — APT
    # primaryIds carries int IDs only).
    targeted = [sh for sh in loaded.open_shutters if sh.role == "target"]
    assert any(sh.target_id == "123456" for sh in targeted), (
        "target_id should round-trip via primaryIds positional alignment"
    )


def test_legacy_schema_still_loads(tmp_path):
    """Sessions written by the old (pre-MPT-compatible) format must still load."""
    p = tmp_path / "legacy.json"
    p.write_text(json.dumps({
        "version": "1.0",
        "pointing": {"ra_deg": 90.05, "dec_deg": -20.13, "apa_v3_deg": 63.0},
        "instrument": {"disperser": "PRISM", "filter": "CLEAR", "slitlet_height": 3},
        "open_shutters": [
            {"q": 2, "d": 75, "s": 128, "target_id": None, "role": "manual"},
            {"q": 2, "d": 75, "s": 129, "target_id": None, "role": "manual"},
        ],
        "image_path": "/tmp/some.jpg",
    }))
    s = import_session_json(str(p))
    assert s.disperser == "PRISM" and s.filter_name == "CLEAR"
    assert s.pa_v3_deg == 63.0
    assert len(s.open_shutters) == 2
    assert s.image_path == "/tmp/some.jpg"


def test_new_session_loads_via_parse_mpt_json(tmp_path):
    """A session.json written by the new exporter must be loadable through
    the MPT-side path (Load plan from JSON), since it is now valid MPT
    plan JSON with vMPT extras."""
    from app.mpt_io import parse_mpt_json
    s = _sample_session()
    p = tmp_path / "session.json"
    export_session_json(s, str(p))
    plans = parse_mpt_json(str(p))
    assert len(plans) == 1
    plan = plans[0]
    assert plan.grating == "G395M" and plan.filter_name == "F290LP"
    # The flat 4-shutter list should compress to 2 slitlets: a 3-shutter
    # run at (q=2, d=200, s=85..87) and a 1-shutter at (q=3, d=100, s=42).
    assert len(plan.slitlets) == 2
    assert plan.n_open_shutters == 4
