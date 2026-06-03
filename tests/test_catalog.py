"""Tests for vmpt.catalog."""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from vmpt.catalog import Catalog, catalog_in_view, load_catalog


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


def test_load_csv_ra_dec_aliases(tmp_path):
    """RA/Dec column names with mixed case + bracketed/underscored unit
    suffixes must all match — APT, eMPT, JADES, and SExtractor each
    use a different spelling."""
    for header in (
        "ID,RA[deg],DEC[deg]",
        "ID,ra_deg,dec_deg",
        "ID,RA (deg),Dec (deg)",
        "ID,RAJ2000,DEJ2000",
        "ID,Right Ascension,Declination",
    ):
        p = tmp_path / f"cat_{abs(hash(header))}.csv"
        p.write_text(header + "\n1,53.0,-27.7\n2,53.1,-27.8\n")
        cat = load_catalog(str(p))
        assert np.allclose(cat.ra_deg, [53.0, 53.1]), header
        assert np.allclose(cat.dec_deg, [-27.7, -27.8]), header


def test_load_csv_id_aliases(tmp_path):
    """ID column accepts several common aliases."""
    for header, want in (
        ("source_id,RA,DEC", [1, 2]),
        ("NO,RA,DEC", [1, 2]),
        ("ObjID,RA,DEC", [1, 2]),
        ("SrcID,RA,DEC", [1, 2]),
    ):
        p = tmp_path / f"cat_{abs(hash(header))}.csv"
        p.write_text(header + "\n1,53.0,-27.7\n2,53.1,-27.8\n")
        cat = load_catalog(str(p))
        assert list(cat.ids) == want, header


def test_load_csv_no_id_synthesises_sequential(tmp_path):
    """When the catalog has no recognizable ID column, vMPT fakes
    sequential IDs 1..N so downstream code still has something to
    auto-tag slitlets with."""
    p = tmp_path / "no_id.csv"
    p.write_text("RA,DEC\n53.0,-27.7\n53.1,-27.8\n53.2,-27.9\n")
    cat = load_catalog(str(p))
    assert list(cat.ids) == [1, 2, 3]


def test_load_csv_name_used_as_id_when_numeric(tmp_path):
    """If the only candidate ID column is `name` (or `label` / `tag`),
    accept it — but only when the values are integers. String names
    must NOT be silently treated as numeric IDs."""
    # Numeric `name` — accepted as ID.
    p1 = tmp_path / "numeric_name.csv"
    p1.write_text("name,RA,DEC\n42,53.0,-27.7\n7,53.1,-27.8\n")
    cat1 = load_catalog(str(p1))
    assert list(cat1.ids) == [42, 7]
    # String `name` — falls through to sequential IDs.
    p2 = tmp_path / "string_name.csv"
    p2.write_text("name,RA,DEC\nNGC-123,53.0,-27.7\nM31,53.1,-27.8\n")
    cat2 = load_catalog(str(p2))
    assert list(cat2.ids) == [1, 2]


def test_load_csv_ids_above_1e7_are_mod_clamped(tmp_path):
    """JADES-style 8–9 digit IDs are mod'd to fit APT's compact integer
    space (anything ≥ 10⁷ → id % 10⁷). Smaller IDs pass through."""
    p = tmp_path / "big_ids.csv"
    p.write_text(
        "ID,RA,DEC\n"
        "12345678,53.0,-27.7\n"      # → 2345678
        "987654321,53.1,-27.8\n"     # → 7654321
        "100,53.2,-27.9\n"           # → 100 (unchanged)
        "10000000,53.3,-28.0\n"      # → 0 (boundary case)
    )
    cat = load_catalog(str(p))
    assert list(cat.ids) == [2345678, 7654321, 100, 0]


def test_load_csv_weight_column_detected(tmp_path):
    """A `weight` column is now a first-class field on `Catalog`, the
    sibling of `priority`. Loader picks it up via several aliases."""
    for header_w in ("weight", "w", "Wt"):
        p = tmp_path / f"w_{header_w}.csv"
        p.write_text(
            f"ID,RA,DEC,priority,{header_w}\n"
            f"1,53.0,-27.7,1,5\n"
            f"2,53.1,-27.8,2,3\n"
            f"3,53.2,-27.9,3,\n"
        )
        cat = load_catalog(str(p))
        assert np.allclose(cat.weight[:2], [5, 3]), header_w
        # Empty cell → NaN
        assert np.isnan(cat.weight[2]), header_w


def test_load_csv_weight_not_claimed_by_priority(tmp_path):
    """A catalog with only a `weight` column (no `priority`) must NOT
    silently treat the weights as priorities. Both fields end up with
    the right arrays — priorities NaN, weights populated."""
    p = tmp_path / "weight_only.csv"
    p.write_text("ID,RA,DEC,weight\n1,53.0,-27.7,5\n2,53.1,-27.8,3\n")
    cat = load_catalog(str(p))
    assert np.all(np.isnan(cat.priority))
    assert np.allclose(cat.weight, [5, 3])


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
