"""Tests for the APT/MPT export-bundle refinements:
  • target-prefixed filenames + matching catalog.name
  • plan name `vmpt-<target>-<stamp>` and 5-dp aperturePA
  • .cat Weight (weight, fallback priority) / Primary 1-0 / Magnitude / Redshift
  • generated bundle README content
"""

from __future__ import annotations

from vmpt.coords import V3_IDL_Y_ANGLE
from vmpt.empt_io import write_mpt_catalog
from vmpt.session_io import (
    Session,
    apt_catalog_basename,
    apt_plan_basename,
    bundle_readme_text,
    _build_mpt_payload,
)


# ── filename helpers ──────────────────────────────────────────────────────

def test_apt_basenames_target_prefixed():
    assert apt_catalog_basename("/data/rxcj2211_targets.cat") == \
        "rxcj2211_targets_APT_catalog"
    assert apt_plan_basename("/data/rxcj2211_targets.cat") == \
        "rxcj2211_targets_MPT_plan"


def test_apt_basenames_fallback_when_no_catalog():
    assert apt_catalog_basename(None).endswith("_APT_catalog")
    assert apt_plan_basename(None).endswith("_MPT_plan")


# ── plan payload: name + aperturePA + catalog.name ────────────────────────

def _session(**kw):
    base = dict(
        pointing_ra_deg=53.1, pointing_dec_deg=-27.8, pa_v3_deg=70.43,
        disperser="PRISM", filter_name="CLEAR", slitlet_height=3,
        open_shutters=[],
        catalog_path="/data/rxcj2211_targets.cat",
        created="2026-06-15T12:19:06Z",
        name=None,
    )
    base.update(kw)
    return Session(**base)


def test_plan_name_is_informative():
    payload = _build_mpt_payload(_session())
    assert payload["name"] == "vmpt-rxcj2211_targets-20260615T12:19:06Z"


def test_explicit_session_name_is_preserved():
    payload = _build_mpt_payload(_session(name="my-custom-plan"))
    assert payload["name"] == "my-custom-plan"


def test_aperture_pa_rounded_to_5dp():
    sess = _session()
    payload = _build_mpt_payload(sess)
    expected = round((sess.pa_v3_deg + V3_IDL_Y_ANGLE) % 360.0, 5)
    assert payload["aperturePA"] == expected
    # No float-arithmetic noise: rounding is idempotent.
    assert payload["aperturePA"] == round(payload["aperturePA"], 5)


def test_catalog_name_matches_apt_basename():
    payload = _build_mpt_payload(_session())
    assert payload["catalog"]["name"] == "rxcj2211_targets_APT_catalog"
    assert payload["catalog"]["primariesName"] == "rxcj2211_targets_APT_catalog"


# ── README ────────────────────────────────────────────────────────────────

def test_bundle_readme_has_the_three_steps():
    txt = bundle_readme_text(
        catalog_filename="rxcj2211_targets_APT_catalog.cat",
        catalog_name="rxcj2211_targets_APT_catalog",
        plan_filename="rxcj2211_targets_MPT_plan.json",
        n_configs=2,
    )
    assert "Import MSA Source Catalog" in txt
    assert "Whitespace Separated" in txt
    assert "Import Plan(s)" in txt
    assert "Create Observation" in txt
    assert "3-Shutter Slitlet" in txt or "3-shutter" in txt.lower()
    # references the actual files + catalog name
    assert "rxcj2211_targets_APT_catalog.cat" in txt
    assert "rxcj2211_targets_MPT_plan.json" in txt
    # multi-config note present when n_configs > 1
    assert "config_1/" in txt and "config_2/" in txt


def test_bundle_readme_single_config_has_no_config_subdir_row():
    txt = bundle_readme_text(
        catalog_filename="a_APT_catalog.cat", catalog_name="a_APT_catalog",
        plan_filename="a_MPT_plan.json", n_configs=1,
    )
    assert "config_1/" not in txt


# ── .cat columns: Weight / Primary / Magnitude / Redshift ─────────────────

def test_weight_uses_weight_key_not_priority(tmp_path):
    out = tmp_path / "c.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "ra_deg": 1.0, "dec_deg": 2.0,
         "Weight": 42, "Pr": 7, "Primary": 1},
    ])
    cols = out.read_text().strip().splitlines()[1].split("\t")
    # ID RA DEC Weight Primary Label  (no mag/z here)
    assert int(cols[3]) == 42          # Weight column = the Weight value
    assert int(cols[4]) == 1           # Primary


def test_weight_falls_back_to_priority(tmp_path):
    out = tmp_path / "c.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "ra_deg": 1.0, "dec_deg": 2.0, "Pr": 7},  # no Weight
    ])
    cols = out.read_text().strip().splitlines()[1].split("\t")
    assert int(cols[3]) == 7           # falls back to priority


def test_primary_flag_distinguishes_catalog_vs_synth(tmp_path):
    out = tmp_path / "c.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "ra_deg": 1.0, "dec_deg": 2.0, "Pr": 3, "Primary": 1},
        {"No_cat": 2, "ra_deg": 1.0, "dec_deg": 2.0, "Pr": 5, "Primary": 0,
         "label": "vMPT_synth"},
    ])
    rows = [ln.split("\t") for ln in out.read_text().strip().splitlines()[1:]]
    assert int(rows[0][4]) == 1
    assert int(rows[1][4]) == 0


def test_magnitude_redshift_columns_appear_when_present(tmp_path):
    out = tmp_path / "c.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "ra_deg": 1.0, "dec_deg": 2.0, "Weight": 1, "Primary": 1,
         "Magnitude": 27.3, "Redshift": 6.5},
        {"No_cat": 2, "ra_deg": 1.0, "dec_deg": 2.0, "Weight": 1, "Primary": 0,
         "label": "vMPT_synth"},   # missing mag/z -> NaN
    ])
    lines = out.read_text().strip().splitlines()
    header = lines[0].lstrip("#").strip().split("\t")
    assert header == ["ID", "RA", "DEC", "Weight", "Primary",
                      "Magnitude", "Redshift", "Label"], header
    r1 = lines[1].split("\t")
    assert float(r1[5]) == 27.3 and float(r1[6]) == 6.5
    r2 = lines[2].split("\t")
    # Missing cells use FINITE sentinels (APT rejects NaN), not "NaN".
    assert r2[5] == "99.9" and r2[6] == "-1.0"
    assert float(r2[5]) == 99.9 and float(r2[6]) == -1.0
    assert len(r1) == len(r2) == 8


def test_no_mag_z_columns_when_absent(tmp_path):
    """Back-compat: with no Magnitude/Redshift anywhere, the header stays
    the original 6-column form (so existing importers are unaffected)."""
    out = tmp_path / "c.cat"
    write_mpt_catalog(str(out), [
        {"No_cat": 1, "ra_deg": 1.0, "dec_deg": 2.0, "Pr": 1},
    ])
    header = out.read_text().splitlines()[0].lstrip("#").strip().split("\t")
    assert header == ["ID", "RA", "DEC", "Weight", "Primary", "Label"]
