"""Tests for vmpt.catalog_ops — the weight ↔ priority helpers used by
the catalog editor."""

from __future__ import annotations

import numpy as np
import pytest

from vmpt.catalog_ops import (
    compute_priorities_from_weights,
    compute_weights_from_priorities,
)


# ---------------------------------------------------------------------
# compute_weights_from_priorities
# ---------------------------------------------------------------------


def test_w_from_p_returns_none_when_no_finite_priorities():
    assert compute_weights_from_priorities(["", "nan", "  "]) is None


def test_w_from_p_lowest_class_gets_one():
    # 4 sources, all in the same priority class — they get weight 1.
    out = compute_weights_from_priorities(["1", "1", "1", "1"])
    assert out == ["1", "1", "1", "1"]


def test_w_from_p_two_classes_strict_dominance():
    """Single P0 must outweigh many P1 — strict-dominance constraint.

    With 3 P1s and 1 P0:  w(P1)=1 → w(P0) > 1 (strict gt) AND
    1*w(P0) > 3*1=3 → w(P0) ≥ 4. We compute 4."""
    # Priority strings: rows 1,2,3 → P1 (lower); row 4 → P0 (higher).
    pris = ["1", "1", "1", "0"]
    out = compute_weights_from_priorities(pris)
    # First three rows are P=1 → weight 1; fourth is P=0 → weight 4.
    assert out[:3] == ["1", "1", "1"]
    assert out[3] == "4"


def test_w_from_p_inequalities_hold_across_all_classes():
    """For a 3-class catalog, both inequalities must hold pairwise."""
    # 5 sources at P=2 (lowest), 3 at P=1, 2 at P=0 (highest).
    pris = ["2", "2", "2", "2", "2", "1", "1", "1", "0", "0"]
    out = compute_weights_from_priorities(pris)
    # Pull the class-level weights back out.
    w_by_class = {}
    for s, p in zip(out, pris):
        w_by_class[int(p)] = int(s)
    # Strictly increasing as p decreases.
    assert w_by_class[0] > w_by_class[1] > w_by_class[2]
    # Class sum strictly increases as p decreases.
    counts = {2: 5, 1: 3, 0: 2}
    sums = {p: counts[p] * w_by_class[p] for p in (2, 1, 0)}
    assert sums[0] > sums[1] > sums[2]


def test_w_from_p_handles_nan_priorities():
    """Rows with NaN priority get empty-string weight, others normal."""
    out = compute_weights_from_priorities(["0", "1", "", "nan", "1"])
    assert out[2] == "" and out[3] == ""
    # The two P=1 rows share a weight; the P=0 row has strictly more.
    p1_w = int(out[1])
    p0_w = int(out[0])
    assert p0_w > p1_w


# ---------------------------------------------------------------------
# compute_priorities_from_weights
# ---------------------------------------------------------------------


def test_p_from_w_returns_none_when_no_finite_weights():
    assert compute_priorities_from_weights(["", "nan"]) is None


def test_p_from_w_largest_weight_becomes_priority_1():
    out = compute_priorities_from_weights(["5", "3", "1", "5"])
    # Three distinct weights → three priority classes.
    # 5 → 1, 3 → 2, 1 → 3
    assert out == ["1", "2", "3", "1"]


def test_p_from_w_dense_priority_index():
    """Unique weights map to contiguous 1..K priorities."""
    out = compute_priorities_from_weights(["10", "10", "5", "2", "2", "2"])
    # 10 → 1, 5 → 2, 2 → 3
    assert out == ["1", "1", "2", "3", "3", "3"]


def test_p_from_w_blanks_pass_through():
    out = compute_priorities_from_weights(["5", "", "3", "nan"])
    assert out[1] == "" and out[3] == ""
    assert out[0] == "1" and out[2] == "2"


# ---------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------


def test_round_trip_p_to_w_to_p_preserves_class_ordering():
    """Compute w from p, then p from w. The resulting priorities have
    the SAME class ordering (larger weight = smaller p) even if the
    absolute priority numbers change.
    """
    pris = ["2", "2", "2", "1", "1", "0"]
    w = compute_weights_from_priorities(pris)
    p2 = compute_priorities_from_weights(w)
    # Each class collapses to a single weight, so:
    #   p=0 → highest weight → p2=1
    #   p=1 → middle         → p2=2
    #   p=2 → lowest         → p2=3
    expected = ["3", "3", "3", "2", "2", "1"]
    assert p2 == expected
