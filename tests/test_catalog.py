"""Tests for app.catalog."""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.catalog import Catalog, catalog_in_view, load_catalog


def test_load_csv(tmp_path):
    p = tmp_path / "cat.csv"
    p.write_text("ID,RA,DEC,priority\n1,53.0,-27.7,1\n2,53.1,-27.8,2\n3,53.2,-27.9,3\n")
    cat = load_catalog(str(p))
    assert len(cat.ids) == 3
    assert np.allclose(cat.ra_deg, [53.0, 53.1, 53.2])
    assert np.allclose(cat.dec_deg, [-27.7, -27.8, -27.9])
    assert np.allclose(cat.priority, [1, 2, 3])
    assert np.all(np.isnan(cat.mag))
    assert np.all(np.isnan(cat.z))


def test_load_csv_priority_class_strings(tmp_path):
    """Catalogs in the wild use letter-prefixed priority *classes*
    (`P0` = highest, `P1` = next, …). The loader must extract the
    numeric portion rather than throwing ValueError trying to coerce
    `"P0"` to float."""
    p = tmp_path / "classes.csv"
    p.write_text(
        "ID,RA,DEC,priority\n"
        "obj1,53.0,-27.7,P0\n"
        "obj2,53.1,-27.8,P1\n"
        "obj3,53.2,-27.9,P2\n"
    )
    cat = load_catalog(str(p))
    assert len(cat.ra_deg) == 3
    assert np.allclose(cat.priority, [0, 1, 2])


def test_load_csv_masked_mag_and_z_become_nan(tmp_path):
    """When numeric columns have empty cells (astropy masks them on read),
    the loader should yield NaN rather than 0 in the masked positions."""
    p = tmp_path / "masked.csv"
    p.write_text(
        "ID,RA,DEC,mag,z\n"
        "obj1,53.0,-27.7,23.5,6.07\n"
        "obj2,53.1,-27.8,,\n"
        "obj3,53.2,-27.9,21.0,5.5\n"
    )
    cat = load_catalog(str(p))
    assert cat.mag[0] == 23.5
    assert np.isnan(cat.mag[1])
    assert cat.mag[2] == 21.0
    assert cat.z[0] == 6.07
    assert np.isnan(cat.z[1])
    assert cat.z[2] == 5.5


def test_load_ascii_empt_style(tmp_path):
    p = tmp_path / "observed_targets.cat"
    p.write_text(
        "# No   No_sub      No_cat    Pr    RA[deg]     Dec[deg]\n"
        "   1        1        14170   1   53.1633910  -27.7756740\n"
        "   2        1         8821   2   53.1641205  -27.7748813\n"
    )
    cat = load_catalog(str(p))
    assert len(cat.ra_deg) == 2
    assert abs(cat.ra_deg[0] - 53.1633910) < 1e-6
    assert abs(cat.dec_deg[1] - (-27.7748813)) < 1e-6
    # Pr present
    assert np.allclose(cat.priority, [1, 2])


def test_view_bbox():
    cat = Catalog(
        ids=np.array([1, 2, 3, 4]),
        ra_deg=np.array([53.0, 53.5, 54.0, 55.0]),
        dec_deg=np.array([-27.7, -27.8, -27.5, -28.0]),
        priority=np.full(4, np.nan),
        mag=np.full(4, np.nan),
        z=np.full(4, np.nan),
        label=np.array(["", "", "", ""], dtype=object),
        source_path="",
    )
    mask = catalog_in_view(cat, 53.0, 54.0, -27.9, -27.6)
    assert mask.tolist() == [True, True, False, False]


def test_view_bbox_ra_wrap():
    cat = Catalog(
        ids=np.array([1, 2, 3]),
        ra_deg=np.array([359.5, 0.5, 180.0]),
        dec_deg=np.array([0.0, 0.0, 0.0]),
        priority=np.full(3, np.nan),
        mag=np.full(3, np.nan),
        z=np.full(3, np.nan),
        label=np.array(["", "", ""], dtype=object),
        source_path="",
    )
    # RA window 359 → 1 (wraps)
    mask = catalog_in_view(cat, 359.0, 1.0, -1.0, 1.0)
    assert mask.tolist() == [True, True, False]
