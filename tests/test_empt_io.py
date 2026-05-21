"""Tests for app.empt_io (eMPT-compatible exporters)."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from app.empt_io import (
    OpenShutter,
    Pointing,
    parse_pointing_summary_txt,
    parse_shutter_mask_csv,
    write_observed_targets_cat,
    write_pointing_summary_txt,
    write_shutter_mask_csv,
)

REF_DIR = Path(__file__).resolve().parents[1] / "refs/eMPT_v1/trial_00_ref/m_pick_output/pointing_100"
REF_CSV = REF_DIR / "shutter_mask.csv"
REF_CAT = REF_DIR / "observed_targets.cat"


# ---------------------------------------------------------------------------
# shutter_mask.csv
# ---------------------------------------------------------------------------


def test_csv_header_byte_exact(tmp_path):
    operable = np.ones((4, 171, 365), dtype=bool)
    reason = np.zeros((4, 171, 365), dtype=np.int8)
    out = tmp_path / "shutter_mask.csv"
    write_shutter_mask_csv(str(out), [], operable, reason)
    ours = out.read_text().splitlines()[0]
    ref = REF_CSV.read_text().splitlines()[0]
    assert ours == ref, "header line must match reference byte-for-byte"


def test_csv_grid_shape(tmp_path):
    operable = np.ones((4, 171, 365), dtype=bool)
    reason = np.zeros((4, 171, 365), dtype=np.int8)
    out = tmp_path / "shutter_mask.csv"
    write_shutter_mask_csv(str(out), [], operable, reason)
    lines = out.read_text().splitlines()
    assert len(lines) == 731
    for i, line in enumerate(lines[1:], start=1):
        cells = line.split(",")
        assert len(cells) == 342, f"row {i} has {len(cells)} cells, expected 342"


def test_operability_roundtrip_against_reference(tmp_path):
    """Parse the reference CSV, re-write with the same operability, parse again,
    and verify the operability/reason arrays survive the round-trip."""
    operable_ref, reason_ref, _ = parse_shutter_mask_csv(str(REF_CSV))
    # Sanity check vs the per-quadrant counts established during recon:
    #   Q1 x=15520 s=6 ; Q2 x=14517 s=3 ; Q3 x=15637 s=12 ; Q4 x=18115 s=1
    counts = {}
    for q in range(4):
        nx = int((reason_ref[q] == 1).sum())
        ns = int((reason_ref[q] == 2).sum())
        counts[q + 1] = (nx, ns)
    assert counts[1] == (15520, 6)
    assert counts[2] == (14517, 3)
    assert counts[3] == (15637, 12)
    assert counts[4] == (18115, 1)

    out = tmp_path / "roundtrip.csv"
    write_shutter_mask_csv(str(out), [], operable_ref, reason_ref)
    operable_rt, reason_rt, _ = parse_shutter_mask_csv(str(out))
    assert np.array_equal(operable_rt, operable_ref)
    assert np.array_equal(reason_rt, reason_ref)


def test_csv_writes_open_shutter():
    """Open shutters should appear as '0' in the right (q,s,d) cell."""
    operable = np.ones((4, 171, 365), dtype=bool)
    reason = np.zeros((4, 171, 365), dtype=np.int8)
    open_set = [
        OpenShutter(q=1, s=93, d=148),
        OpenShutter(q=3, s=56, d=46),
        OpenShutter(q=4, s=22, d=333),
    ]
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "m.csv"
        write_shutter_mask_csv(str(out), open_set, operable, reason)
        _, _, recovered = parse_shutter_mask_csv(str(out))
        keys = {(sh.q, sh.s, sh.d) for sh in recovered}
        for sh in open_set:
            assert (sh.q, sh.s, sh.d) in keys


# ---------------------------------------------------------------------------
# observed_targets.cat
# ---------------------------------------------------------------------------


def test_observed_targets_first_row(tmp_path):
    out = tmp_path / "observed_targets.cat"
    write_observed_targets_cat(
        str(out),
        [{"No_sub": 1, "No_cat": 14170, "Pr": 1, "ra_deg": 53.1633910, "dec_deg": -27.7756740}],
    )
    ours = out.read_text().splitlines()
    ref = REF_CAT.read_text().splitlines()
    # Header line should match (ASCII text — no tabs in the reference).
    assert ours[0] == ref[0]
    # First data row must match byte-for-byte.
    assert ours[1] == ref[1], f"\nours={ours[1]!r}\nref={ref[1]!r}"


# ---------------------------------------------------------------------------
# pointing_summary.txt
# ---------------------------------------------------------------------------


def test_pointing_summary_roundtrip(tmp_path):
    p = Pointing(ra_deg=53.1409714, dec_deg=-27.7919712,
                 apa_v3_deg=321.004456, pa_ap_deg=99.579041)
    out = tmp_path / "pointing_summary.txt"
    write_pointing_summary_txt(str(out), p, disperser="PRISM", filter_name="CLEAR",
                               n_targets_total=13119, n_targets_accepted=155)
    parsed = parse_pointing_summary_txt(str(out))
    assert parsed["ra_deg"] == pytest.approx(53.1409714, abs=1e-7)
    assert parsed["dec_deg"] == pytest.approx(-27.7919712, abs=1e-7)
    assert parsed["pa_ap_deg"] == pytest.approx(99.579041, abs=1e-6)
    assert parsed["pa_v3_deg"] == pytest.approx(321.004456, abs=1e-6)


# ---------------------------------------------------------------------------
# MPT-importable catalog (APT format)
# ---------------------------------------------------------------------------


def test_write_mpt_catalog_header_and_data(tmp_path):
    from app.empt_io import write_mpt_catalog
    out = tmp_path / "MPT_catalog.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "Pr": 1, "ra_deg": 39.9826125,   "dec_deg": -1.5916444},
        {"No_cat": 2, "Pr": 5, "ra_deg": 39.98701666,  "dec_deg": -1.5891333,
         "label": "vMPT_synth"},
    ])
    text = out.read_text()
    lines = text.strip().splitlines()
    # Header must start with "#" and use the EXACT column names APT
    # recognizes — bare "RA"/"DEC", not "RA[deg]" / "DEC[deg]".
    assert lines[0].startswith("#")
    header_tokens = lines[0].lstrip("#").strip().split("\t")
    assert header_tokens == ["ID", "RA", "DEC", "Weight", "Primary", "Label"], header_tokens
    # 1 header + 2 data rows; tab-separated so APT's column parser binds
    # each data field to one labelled column.
    assert len(lines) == 3
    for data_line in lines[1:]:
        tokens = data_line.split("\t")
        assert len(tokens) == 6, f"expected 6 tab-separated fields, got {tokens}"
    # First data row carries the ID, RA, Dec, weight, Primary=1, Label=real
    id1, ra1, dec1, w1, prim1, lab1 = lines[1].split("\t")
    assert int(id1) == 1
    assert float(ra1) == 39.9826125
    assert float(dec1) == -1.5916444
    assert int(w1) == 1
    assert int(prim1) == 1
    assert lab1 == "real"
    # Second data row was tagged as synthesized
    _, _, _, w2, _, lab2 = lines[2].split("\t")
    assert int(w2) == 5
    assert lab2 == "vMPT_synth"


def test_write_mpt_catalog_preserves_string_ids(tmp_path):
    """If the input catalog used non-integer IDs (e.g. "RJ0600-10274-P0"),
    the output MPT catalog should preserve those IDs verbatim — the
    user's downstream tooling depends on the same identifiers carrying
    through."""
    from app.empt_io import write_mpt_catalog
    out = tmp_path / "RJ0600-targets.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": "RJ0600-10274-P0", "ra_deg": 90.039845, "dec_deg": -20.136384,
         "label": "real"},
        {"No_cat": "RJ0600-8846-P0",  "ra_deg": 90.038126, "dec_deg": -20.140741,
         "label": "real"},
        # A synthesized row with an int ID
        {"No_cat": 9001, "ra_deg": 90.04, "dec_deg": -20.13, "label": "vMPT_synth"},
    ])
    lines = out.read_text().strip().splitlines()
    assert lines[1].split("\t")[0] == "RJ0600-10274-P0"
    assert lines[2].split("\t")[0] == "RJ0600-8846-P0"
    assert lines[3].split("\t")[0] == "9001"


def test_write_mpt_catalog_distinguishes_real_vs_synth(tmp_path):
    """The output catalog must let downstream tools tell which rows came
    from the user's input catalog and which were synthesized by vMPT
    for unmatched slitlets."""
    from app.empt_io import write_mpt_catalog
    out = tmp_path / "MPT_catalog.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 100, "ra_deg": 53.16, "dec_deg": -27.77, "label": "real"},
        {"No_cat": 101, "ra_deg": 53.16, "dec_deg": -27.77, "label": "vMPT_synth"},
        {"No_cat": 102, "ra_deg": 53.16, "dec_deg": -27.77},  # default to "real"
    ])
    data = out.read_text().strip().splitlines()[1:]
    labels = [line.split("\t")[5] for line in data]
    assert labels == ["real", "vMPT_synth", "real"]


def test_write_mpt_catalog_uses_only_jdox_column_names(tmp_path):
    """Guard against accidentally reintroducing unit-suffixed column
    names like `RA[deg]` — APT's column-detection matches the bare
    labels listed on the JDox MPT Catalogs page."""
    from app.empt_io import write_mpt_catalog
    out = tmp_path / "MPT_catalog.cat"
    write_mpt_catalog(str(out), [{"No_cat": 1, "ra_deg": 0.0, "dec_deg": 0.0}])
    text = out.read_text()
    assert "[deg]" not in text, "do not write `[deg]` suffixes — not a recognized APT label"
    assert "Priority" not in text, "`Weight` is the JDox-recognized name, not `Priority`"
