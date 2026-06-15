"""Tests for vmpt.catalog.evaluate_catalog_condition — the sandboxed
row-selection expression evaluator used by the catalog editor's
condition+value max_configs rule."""

from __future__ import annotations

import numpy as np
import pytest

from vmpt.catalog import evaluate_catalog_condition as ev


COLS = {
    "mag": ["26.0", "27.5", "28.2", ""],     # blank → NaN
    "z": ["5.0", "6.5", "7.1", "2.0"],
    "role": ["science", "sky", "science", "filler"],
}


# ── happy paths ──────────────────────────────────────────────────────────

def test_simple_numeric_comparison():
    mask = ev("mag > 27", COLS)
    assert mask.dtype == bool
    assert mask.tolist() == [False, True, True, False]  # blank mag → False


def test_compound_and():
    mask = ev("(mag > 27) & (z > 6)", COLS)
    assert mask.tolist() == [False, True, True, False]


def test_compound_or():
    mask = ev("(z > 7) | (role == 'sky')", COLS)
    assert mask.tolist() == [False, True, True, False]


def test_invert():
    mask = ev("~(z > 6)", COLS)
    assert mask.tolist() == [True, False, False, True]


def test_string_equality():
    mask = ev("role == 'science'", COLS)
    assert mask.tolist() == [True, False, True, False]


def test_membership_isin():
    mask = ev("isin(role, ('sky', 'filler'))", COLS)
    assert mask.tolist() == [False, True, False, True]


def test_whitelisted_function():
    cols = {"flux": ["10", "100", "1000", "1"]}
    mask = ev("log10(flux) >= 2", cols)
    assert mask.tolist() == [False, True, True, False]


def test_blank_is_nan_not_match():
    # Blank numeric entries must never satisfy a numeric comparison.
    mask = ev("mag < 100", COLS)
    assert mask.tolist() == [True, True, True, False]


def test_scalar_result_broadcasts():
    # A constant expression still yields one value per row.
    mask = ev("z > 0", {"z": ["1", "2", "3"]})
    assert mask.tolist() == [True, True, True]


# ── error reporting (raised BEFORE any edit is applied) ───────────────────

def test_empty_condition_errors():
    with pytest.raises(ValueError, match="empty"):
        ev("", COLS)


def test_syntax_error_message():
    with pytest.raises(ValueError, match="[Ss]yntax"):
        ev("mag > ", COLS)


def test_unknown_column_message():
    with pytest.raises(ValueError, match="Unknown column 'mag_f444w'"):
        ev("mag_f444w > 27", COLS)


def test_attribute_access_blocked():
    with pytest.raises(ValueError):
        ev("mag.__class__", COLS)


def test_non_whitelisted_call_blocked():
    with pytest.raises(ValueError):
        ev("foo(mag)", COLS)


def test_dunder_import_blocked():
    # Even if it parsed, builtins are stripped; the name check fires first.
    with pytest.raises(ValueError):
        ev("__import__('os')", COLS)


def test_and_or_keyword_gives_friendly_hint():
    with pytest.raises(ValueError, match=r"&|\|"):
        ev("(mag > 27) and (z > 6)", COLS)


def test_non_boolean_result_rejected():
    # A bare arithmetic expression is not a boolean selection.
    with pytest.raises(ValueError):
        ev("mag + z", {"mag": ["1", "2"], "z": ["3", "4"]})
