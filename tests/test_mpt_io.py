"""Tests for the MPT JSON + shutter CSV importer."""
from __future__ import annotations

import json
import os
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from app.mpt_io import (
    list_mpt_plans_in_aptx,
    parse_mpt_json,
    parse_mpt_json_in_aptx,
    parse_shutter_csv,
)

EX = _ROOT / "example_a370"
MPT_JSON = EX / "a370 plan 1, a370 plan 2, a370 plan 3 (Obs 46).json"
SH_CSV = EX / "1208.p1c1 : a370 plan 1e1n1.csv"


@pytest.mark.skipif(not MPT_JSON.exists(), reason="MPT example missing")
def test_parse_mpt_json_3_plans():
    plans = parse_mpt_json(str(MPT_JSON))
    assert len(plans) == 3
    for p in plans:
        # APA is fixed across plans (single observation)
        assert 40 < p.aperture_pa_deg < 50
        # V3 PA derived: 46.1 - 138.575 = -92.47 mod 360 = 267.53
        assert 260 < p.v3_pa_deg < 275
        assert p.grating == "PRISM"
        assert p.filter_name == "CLEAR"
        assert p.slitlets and len(p.slitlets) > 100
        # Most slitlets are 3-shutter for the a370 PRISM plan; allow some 1-shutter "spare sky" entries
        h_values = {sl.h for sl in p.slitlets}
        assert h_values.issubset({1, 3})


@pytest.mark.skipif(not MPT_JSON.exists(), reason="MPT example missing")
def test_to_open_shutters_unfolds_height():
    plans = parse_mpt_json(str(MPT_JSON))
    p0 = plans[0]
    opens = p0.to_open_shutters()
    # 121 slitlets × 3 shutters per slitlet
    assert len(opens) == 121 * 3 == 363
    # The middle shutter of each slitlet has role 'target', others 'sky'
    target_count = sum(1 for o in opens if o.role == "target")
    sky_count = sum(1 for o in opens if o.role == "sky")
    assert target_count == 121
    assert sky_count == 2 * 121


@pytest.mark.skipif(not SH_CSV.exists(), reason="Shutter CSV example missing")
def test_parse_shutter_csv_matches_json_count():
    """The CSV's open-shutter count must match the plan 1 JSON's slitlet
    unwrapping (121 slitlets × 3 = 363 open shutters)."""
    opens = parse_shutter_csv(str(SH_CSV))
    assert len(opens) == 363
    # All have role 'manual' (CSV doesn't preserve slitlet/target info)
    assert all(o.role == "manual" for o in opens)
    # (q, s, d) values inside valid ranges
    for o in opens:
        assert 1 <= o.q <= 4
        assert 1 <= o.s <= 171
        assert 1 <= o.d <= 365


@pytest.mark.skipif(not (MPT_JSON.exists() and SH_CSV.exists()), reason="examples missing")
def test_json_and_csv_open_shutters_agree():
    """The JSON plan 1 unfolded should match the CSV's open-shutter set."""
    plans = parse_mpt_json(str(MPT_JSON))
    json_opens = {(o.q, o.s, o.d) for o in plans[0].to_open_shutters()}
    csv_opens = {(o.q, o.s, o.d) for o in parse_shutter_csv(str(SH_CSV))}
    assert json_opens == csv_opens, (
        f"JSON-derived ({len(json_opens)}) vs CSV ({len(csv_opens)}) disagree; "
        f"only in JSON: {sorted(json_opens - csv_opens)[:5]}, "
        f"only in CSV: {sorted(csv_opens - json_opens)[:5]}"
    )


def test_aptx_list_and_parse(tmp_path):
    """Build a tiny .aptx-like zip with one MPT JSON inside; verify our
    list + extract path works without needing the network."""
    plan = {
        "instrument": "NIRSpec",
        "name": "test plan",
        "aperturePA": 46.0,
        "theta": 43.5,
        "catalog": {"name": "test_catalog"},
        "configs": [{
            "name": "c1",
            "slitlets": [{"q": 2, "d": 200, "s": 86, "h": 3}],
            "exposures": [{"ra": 10.0, "dec": 20.0,
                           "gratingFilter": "PRISM_CLEAR",
                           "sourceIds": [123]}],
            "primaryIds": [123],
        }],
    }
    aptx = tmp_path / "fake.aptx"
    with zipfile.ZipFile(aptx, "w") as zf:
        zf.writestr("manifest", "tiny")
        zf.writestr("plan_one.json", json.dumps(plan))
        # A non-MPT JSON that should be filtered out
        zf.writestr("BESTPointing1.json", json.dumps({"foo": 1}))
    plans = list_mpt_plans_in_aptx(str(aptx))
    assert plans == ["plan_one.json"]  # BESTPointing filtered (no configs)
    extracted = parse_mpt_json_in_aptx(str(aptx), "plan_one.json")
    assert len(extracted) == 1
    assert extracted[0].n_open_shutters == 3


@pytest.mark.skipif(not os.environ.get("VMPT_NETWORK_TESTS"),
                    reason="set VMPT_NETWORK_TESTS=1 to run live STScI download test")
def test_download_apt_program_live():
    """Live network test. Downloads APT 1208 from STScI and confirms the
    archive contains MPT-format JSON files. Skipped by default."""
    from app.mpt_io import download_apt_program
    aptx = download_apt_program("1208")
    plans = list_mpt_plans_in_aptx(aptx)
    assert plans  # APT 1208 has 40+ embedded plans
    # Pick one and parse it
    first = parse_mpt_json_in_aptx(aptx, plans[0])
    assert first and first[0].slitlets is not None
