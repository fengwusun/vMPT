"""A loaded shutter mask must auto-populate the MPT catalog — no optimizer
run required.

`parse_shutter_csv` produces anonymous opens (``role='manual'``,
``target_id=None``), so before this feature the MPT-catalog viewer was
empty until the optimizer assigned ids. `_retag_manual_opens_from_catalog`
re-derives each raw-mask shutter's catalog source(s) from the live
footprint match, so the viewer fills straight from the mask. Deliberate
picks (manual clicks / optimizer results) must NOT be touched.
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
EX_CAT = Path.home() / ".vmpt/examples/example_r0600/v01_fsun.cat"

pytestmark = pytest.mark.skipif(
    not (EX_JPG.exists() and EX_WCS.exists() and EX_CAT.exists()),
    reason="example_r0600 assets missing",
)


@pytest.fixture()
def app_with_catalog():
    import vmpt.main as m
    from vmpt.image_io import load_jpg_with_sidecar
    from vmpt.catalog import load_catalog

    img = load_jpg_with_sidecar(str(EX_JPG), str(EX_WCS), max_dim=4000)
    H, W = img.shape[:2]
    center = img.wcs.pixel_to_world(W / 2, H / 2)
    m.state["image"] = img
    m.ra_input.value = repr(float(center.ra.deg))
    m.dec_input.value = repr(float(center.dec.deg))
    m.state["ra_deg"] = float(center.ra.deg)
    m.state["dec_deg"] = float(center.dec.deg)
    m.state["pa_v3"] = 137.0
    m.state["catalog"] = load_catalog(str(EX_CAT))
    # Footprint index for the current pointing (core only — no re-tag yet).
    m._rebuild_shutter_catalog_index_core()
    return m


def test_loaded_mask_populates_mpt_catalog(app_with_catalog):
    m = app_with_catalog
    bucket = m.state["shutter_to_catids"]
    assert bucket, "expected some catalog sources to fall in the MSA"

    # Simulate a loaded shutter mask: anonymous opens on the shutters that
    # contain catalog sources (role='manual', no target id) — exactly what
    # parse_shutter_csv yields.
    keys = list(bucket.keys())[:5]
    m._set_open_shutters({
        (q, s, d): m.OpenShutter(q=q, s=s, d=d, target_id=None, role="manual")
        for (q, s, d) in keys
    })

    # Before re-tag: the viewer sees no source ids → no rows.
    assert m._mpt_view_collect_rows() == [], \
        "raw mask should carry no source ids before the auto-match"

    # Auto-match (this is what loading the mask now triggers).
    matched = m._retag_manual_opens_from_catalog()
    assert matched == len(keys), "every seeded shutter should match a source"

    rows = m._mpt_view_collect_rows()
    assert rows, "MPT catalog must populate from the mask, no optimizer run"
    cat_ids = {str(i) for i in m.state["catalog"].ids}
    assert all(r["id"] in cat_ids for r in rows)
    # The opens themselves now carry the catalog ids.
    for (q, s, d), sh in m.state["open_shutters"].items():
        assert m._open_shutter_ids(sh), f"shutter {(q,s,d)} should be tagged"


def test_retag_does_not_clobber_deliberate_picks(app_with_catalog):
    m = app_with_catalog
    bucket = m.state["shutter_to_catids"]
    occupied = set(bucket.keys())
    # A deliberate pick (role='target') with its own id, on a shutter that
    # holds NO catalog source — re-tag must leave it exactly as-is.
    empty_key = next((q, s, d)
                     for q in (1, 2) for s in range(80, 120) for d in range(150, 220)
                     if (q, s, d) not in occupied
                     and m.OPERABLE[q - 1, s - 1, d - 1])
    q, s, d = empty_key
    m._set_open_shutters({
        empty_key: m.OpenShutter(q=q, s=s, d=d, target_id="DELIBERATE",
                                 role="target", target_ids=["DELIBERATE"]),
    })
    m._retag_manual_opens_from_catalog()
    sh = m.state["open_shutters"][empty_key]
    assert sh.role == "target"
    assert sh.target_id == "DELIBERATE"
    assert m._open_shutter_ids(sh) == ["DELIBERATE"], \
        "deliberate pick must not be re-tagged from the footprint"
