"""`_find_msaoper_json` must pick the MOS-current operability reference
(the NRS_MSASPEC branch, latest USEAFTER <= today) using the cached rmap —
not just the highest-numbered file, which can be a higher-numbered
imaging-only delivery. Falls back to the highest-numbered file when no
rmap is cached.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.msa as msa  # noqa: E402

# An rmap where the MOS branch tops out at _0017, but a *higher-numbered*
# imaging-only reference (_0099) also exists — the naive "newest file"
# heuristic would wrongly pick _0099 for MOS.
_RMAP = """\
header = {'instrument': 'NIRSPEC'}
selector = Match({
    'N/A' : UseAfter({
        '2014-01-01 00:00:00' : 'jwst_nirspec_msaoper_0099.json',
    }),
    'NRS_IFU|NRS_MSASPEC' : UseAfter({
        '2024-08-05 00:00:00' : 'jwst_nirspec_msaoper_0014.json',
        '2026-04-14 00:00:00' : 'jwst_nirspec_msaoper_0017.json',
        '2099-01-01 00:00:00' : 'jwst_nirspec_msaoper_0098.json',
    }),
})
"""


def _make_cache(tmp_path, refs, rmap_text):
    ref_dir = tmp_path / "references" / "jwst" / "nirspec"
    ref_dir.mkdir(parents=True)
    for name in refs:
        (ref_dir / name).write_text('{"msaoper": []}')
    if rmap_text is not None:
        map_dir = tmp_path / "mappings" / "jwst"
        map_dir.mkdir(parents=True)
        (map_dir / "jwst_nirspec_msaoper_0019.rmap").write_text(rmap_text)
    return ref_dir


def test_selects_mos_current_not_highest_numbered(tmp_path):
    _make_cache(tmp_path, [
        "jwst_nirspec_msaoper_0014.json",
        "jwst_nirspec_msaoper_0017.json",
        "jwst_nirspec_msaoper_0099.json",  # imaging-only, higher number
    ], _RMAP)
    picked = os.path.basename(msa._find_msaoper_json(str(tmp_path)))
    assert picked == "jwst_nirspec_msaoper_0017.json", (
        f"should pick the MOS-current _0017, not the higher-numbered "
        f"imaging-only _0099 — got {picked}")


def test_ignores_future_useafter(tmp_path):
    # _0098 has a future USEAFTER (2099) — must not be selected yet even
    # though it's in the MOS branch and present in the cache.
    _make_cache(tmp_path, [
        "jwst_nirspec_msaoper_0017.json",
        "jwst_nirspec_msaoper_0098.json",  # future USEAFTER
    ], _RMAP)
    picked = os.path.basename(msa._find_msaoper_json(str(tmp_path)))
    assert picked == "jwst_nirspec_msaoper_0017.json"


def test_falls_back_to_highest_numbered_without_rmap(tmp_path):
    _make_cache(tmp_path, [
        "jwst_nirspec_msaoper_0014.json",
        "jwst_nirspec_msaoper_0017.json",
    ], rmap_text=None)  # no rmap → glob fallback
    picked = os.path.basename(msa._find_msaoper_json(str(tmp_path)))
    assert picked == "jwst_nirspec_msaoper_0017.json"


def test_ensure_current_operability_disabled_is_noop(monkeypatch):
    # With the auto-update disabled, the startup check is a no-op: no
    # network, no raise, returns None (vMPT just uses the cached reference).
    monkeypatch.setenv("VMPT_OPERABILITY_AUTOUPDATE", "0")
    msa._OPER_CHECK_DONE = False
    assert msa.ensure_current_operability() is None
    # And it never runs twice per process (guard) even when enabled.
    monkeypatch.setenv("VMPT_OPERABILITY_AUTOUPDATE", "1")
    msa._OPER_CHECK_DONE = True
    assert msa.ensure_current_operability() is None


def test_bundled_fallback_is_shipped():
    """The compact operability snapshot must exist in vmpt/data so it ships
    via the package-data `data/*.npz` glob."""
    assert (msa._DATA_DIR / "msaoper_fallback.npz").is_file()


def test_bundled_fallback_loads_when_no_crds(monkeypatch):
    """A fresh deploy with NO CRDS reference (no `crds` package / no
    CRDS_PATH / no ~/crds_cache / no network) must still get the failed-/
    stuck-shutter map from the bundled snapshot — not silently fall through
    to "every shutter operable". Regression for the reported deployment bug.
    """
    import numpy as np
    monkeypatch.setattr(msa, "_find_msaoper_json", lambda crds_path=None: None)
    operable, reason = msa.load_operability()
    assert int((reason == 2).sum()) > 0, "snapshot must carry stuck-opens"
    assert int((reason == 1).sum()) > 0, "snapshot must carry failed-closed"
    assert np.array_equal(operable, reason == 0)
    assert not operable.all(), "must NOT be the degenerate all-operable fallback"
    # The source flag lets the UI warn that this isn't the live CRDS reference.
    assert msa.OPERABILITY_SOURCE.startswith("bundled:")
    assert msa.operability_is_current() is False


def test_corrupt_crds_file_falls_back_to_bundle(monkeypatch, tmp_path):
    """If the resolved CRDS msaoper JSON is unreadable, load the bundled
    snapshot instead of crashing or going all-operable."""
    import numpy as np
    bad = tmp_path / "jwst_nirspec_msaoper_9999.json"
    bad.write_text("{ this is not valid json")
    monkeypatch.setattr(msa, "_find_msaoper_json", lambda crds_path=None: str(bad))
    operable, reason = msa.load_operability()
    assert int((reason == 2).sum()) > 0
    assert not operable.all()
