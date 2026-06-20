"""`_load_catalog_from_path` must not blow up (or dump a traceback) when
handed a non-file path — an empty string, a directory like '.', or a
missing file. A stray catalog-path value of '.' used to reach the astropy
reader and raise IsADirectoryError; the loader now skips non-files
cleanly.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402


def test_load_catalog_path_skips_non_files():
    m.state["catalogs"] = []
    for bad in (".", "", "/no/such/catalog_file_xyz.csv", str(_ROOT)):
        # Must neither raise nor add a catalog entry.
        m._load_catalog_from_path(bad)
        assert m.state["catalogs"] == [], (
            f"non-file path {bad!r} should not add a catalog")


def test_load_catalog_path_accepts_a_real_file():
    cat_file = Path.home() / ".vmpt/examples/example_r0600/v01_fsun.cat"
    if not cat_file.exists():
        return  # example assets not present in this environment
    m.state["catalogs"] = []
    m._load_catalog_from_path(str(cat_file))
    assert len(m.state["catalogs"]) == 1, "a real catalog file should load"
    assert len(m.state["catalogs"][0]["catalog"].ra_deg) > 0
