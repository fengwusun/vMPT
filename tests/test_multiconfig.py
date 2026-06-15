"""v1.4.0 multi-config + multi-source-shutter regression tests.

Covers the pieces that are testable without a live Bokeh server:
  • Catalog.max_configs CSV round-trip + alias matching.
  • OpenShutter.target_ids dataclass field.
  • PointingEvaluator budget_remaining mask + the DROP_BUDGET reason.
  • Session multi-config + multi-source round-trip (export → import).
  • _build_mpt_payload emitting one config block per config, with
    co-shutter sources surfaced in sourceIds.

The interactive state-switching / two-pass driver / viewer live in
vmpt.main (Bokeh) and are exercised by the manual Chrome pass.
"""
import os
import tempfile

import numpy as np

from vmpt.catalog import load_catalog, save_catalog
from vmpt.empt_io import OpenShutter
from vmpt.optimizer import PointingEvaluator, DROP_BUDGET
from vmpt.session_io import (
    Session, export_session_json, import_session_json, _build_mpt_payload,
)


# ---------------------------------------------------------------------------
# Catalog.max_configs
# ---------------------------------------------------------------------------

def _write_csv(text: str) -> str:
    d = tempfile.mkdtemp()
    p = os.path.join(d, "c.csv")
    with open(p, "w") as f:
        f.write(text)
    return p


def test_max_configs_loads_and_aliases():
    # `max_obs` is one of the accepted aliases for max_configs.
    p = _write_csv("ID,RA,DEC,priority,max_obs\n"
                   "1,53.1,-27.8,1,1\n2,53.2,-27.7,2,\n3,53.3,-27.6,1,2\n")
    cat = load_catalog(p)
    assert cat.max_configs.shape == (3,)
    assert cat.max_configs[0] == 1
    assert not np.isfinite(cat.max_configs[1])   # blank → NaN (unset)
    assert cat.max_configs[2] == 2


def test_max_configs_absent_is_all_nan():
    p = _write_csv("ID,RA,DEC\n1,53.1,-27.8\n2,53.2,-27.7\n")
    cat = load_catalog(p)
    assert cat.max_configs.shape == (2,)
    assert not np.isfinite(cat.max_configs).any()


def test_max_configs_csv_round_trip():
    p = _write_csv("ID,RA,DEC,max_configs\n1,53.1,-27.8,1\n2,53.2,-27.7,\n")
    cat = load_catalog(p)
    out = os.path.join(os.path.dirname(p), "out.csv")
    save_catalog(cat, out)
    assert "max_configs" in open(out).readline()
    cat2 = load_catalog(out)
    np.testing.assert_array_equal(
        np.nan_to_num(cat.max_configs, nan=-1.0),
        np.nan_to_num(cat2.max_configs, nan=-1.0),
    )


# ---------------------------------------------------------------------------
# OpenShutter.target_ids
# ---------------------------------------------------------------------------

def test_open_shutter_target_ids_default_empty():
    sh = OpenShutter(q=1, s=50, d=12, target_id="A")
    assert sh.target_ids == []
    sh2 = OpenShutter(q=1, s=50, d=12, target_id="A", target_ids=["A", "B"])
    assert sh2.target_ids == ["A", "B"]


# ---------------------------------------------------------------------------
# PointingEvaluator budget mask
# ---------------------------------------------------------------------------

def _toy_field():
    np.random.seed(7)
    n = 40
    ra = 53 + np.random.rand(n) * 0.08
    dec = -27 + np.random.rand(n) * 0.08
    return ra, dec


def test_budget_none_is_identical_to_no_budget():
    ra, dec = _toy_field()
    ev0 = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det0, _, _ = ev0.evaluate(53.05, -26.97, 30.0)
    ev1 = PointingEvaluator(ra, dec, centration="UNCONSTRAINED",
                            budget_remaining=np.ones(len(ra), bool))
    det1, _, _ = ev1.evaluate(53.05, -26.97, 30.0)
    np.testing.assert_array_equal(det0, det1)


def test_budget_drops_exhausted_sources():
    ra, dec = _toy_field()
    ev0 = PointingEvaluator(ra, dec, centration="UNCONSTRAINED")
    det0, _, _ = ev0.evaluate(53.05, -26.97, 30.0)
    placed = np.where(det0)[0]
    assert placed.size >= 3, "need a few detections for the test"
    budget = np.ones(len(ra), bool)
    budget[placed[:3]] = False
    ev1 = PointingEvaluator(ra, dec, centration="UNCONSTRAINED",
                            budget_remaining=budget)
    det1, _, _, reasons = ev1.evaluate_with_reasons(53.05, -26.97, 30.0)
    assert int(det1.sum()) == int(det0.sum()) - 3
    assert reasons[DROP_BUDGET] == 3
    for i in placed[:3]:
        assert not det1[i]


def test_plain_catalog_required_lam_is_1d():
    """Regression: a catalog with no lam_req column must yield a 1D
    required_lam (object array of empty lists), NOT a 2D (n, 0) array
    whose .size is 0 — which used to make the optimizer raise
    'required_lam size 0 != ra_sources size N' for any plain catalog."""
    p = _write_csv("ID,RA,DEC,priority\n1,53.1,-27.8,1\n2,53.2,-27.7,2\n")
    cat = load_catalog(p)
    arr = np.asarray(cat.required_lam, dtype=object)
    assert arr.ndim == 1 and arr.size == len(cat.ra_deg)
    # The evaluator must accept it without raising.
    ev = PointingEvaluator(
        cat.ra_deg, cat.dec_deg, required_lam=cat.required_lam,
        no_gap=cat.no_gap, extend_blue=cat.extend_blue,
        extend_red=cat.extend_red, disperser="PRISM", filt="CLEAR",
    )
    det, _, _ = ev.evaluate(53.13, -27.78, 30.0)
    assert det.shape == (len(cat.ra_deg),)


def test_budget_size_mismatch_raises():
    ra, dec = _toy_field()
    try:
        PointingEvaluator(ra, dec, budget_remaining=np.ones(3, bool))
    except ValueError:
        return
    raise AssertionError("expected ValueError on budget size mismatch")


# ---------------------------------------------------------------------------
# Session multi-config round-trip
# ---------------------------------------------------------------------------

def _multi_config_session():
    c1 = [OpenShutter(q=1, s=50, d=12, target_id="10274", role="target",
                      target_ids=["10274", "10275"])]
    c2 = [OpenShutter(q=2, s=30, d=40, target_id="10280", role="target",
                      target_ids=["10280"])]
    return Session(
        pointing_ra_deg=53.1, pointing_dec_deg=-27.8, pa_v3_deg=120.0,
        disperser="PRISM", filter_name="CLEAR", slitlet_height=3,
        open_shutters=c1, highlighted=[],
        configs=[
            {"name": "Config 1", "ra_deg": 53.1, "dec_deg": -27.8,
             "pa_v3": 120.0, "open_shutters": c1, "highlighted": []},
            {"name": "Config 2", "ra_deg": 53.2, "dec_deg": -27.7,
             "pa_v3": 121.0, "open_shutters": c2, "highlighted": []},
        ],
        active_config=0,
    )


def test_session_multi_config_round_trip():
    sess = _multi_config_session()
    d = tempfile.mkdtemp()
    p = os.path.join(d, "MPT_plan.json")
    export_session_json(sess, p)
    s2 = import_session_json(p)
    assert len(s2.configs) == 2
    # multi-source shutter preserved
    assert s2.configs[0]["open_shutters"][0].target_ids == ["10274", "10275"]
    # per-config pointing preserved
    assert abs(s2.configs[1]["ra_deg"] - 53.2) < 1e-9
    assert abs(s2.configs[1]["dec_deg"] - (-27.7)) < 1e-9


def test_single_config_session_has_no_configs_key():
    """A single-config bundle must NOT emit the multi-config block, so
    pre-1.4 round-trips are unchanged."""
    import json
    c1 = [OpenShutter(q=1, s=50, d=12, target_id="A", role="target",
                      target_ids=["A"])]
    sess = Session(
        pointing_ra_deg=53.1, pointing_dec_deg=-27.8, pa_v3_deg=120.0,
        disperser="PRISM", filter_name="CLEAR", slitlet_height=3,
        open_shutters=c1,
    )
    d = tempfile.mkdtemp()
    p = os.path.join(d, "MPT_plan.json")
    export_session_json(sess, p)
    ws = json.load(open(os.path.join(d, "vMPT_workspace.json")))
    assert "configs" not in ws


def test_mpt_payload_emits_one_block_per_config_with_all_sources():
    sess = _multi_config_session()
    payload = _build_mpt_payload(sess)
    assert len(payload["configs"]) == 2
    assert payload["stats"][0]["numberOfConfigurations"] == 2
    # co-shutter sources both land in sourceIds
    src1 = payload["configs"][0]["exposures"][0]["sourceIds"]
    assert 10274 in src1 and 10275 in src1
    # per-config exposure pointing
    assert abs(payload["configs"][1]["exposures"][0]["ra"] - 53.2) < 1e-9
    assert payload["plannerSpecification"]["searchParameters"][
        "allowMultiSourceShutters"] is True
