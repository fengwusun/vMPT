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
