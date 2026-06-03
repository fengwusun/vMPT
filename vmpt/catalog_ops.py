"""Pure helpers for the catalog editor — computing weight from priority
and vice versa. Lives outside `main.py` so unit tests don't need to
import Bokeh.
"""

from __future__ import annotations

import numpy as np


def compute_weights_from_priorities(priorities) -> list[str] | None:
    """Iterative weight formula.

    For sources at each priority class (smaller p = higher priority),
    find the smallest INTEGER w(p) such that simultaneously

        w(p)         >  w(p+1)
        N(p) * w(p)  >  N(p+1) * w(p+1)

    iterating from the LOWEST priority class (largest p) upward.
    Returns a per-row list of stringified ints (NaN priority → "")
    or None if no finite priorities exist.
    """
    pris: list[float] = []
    for v in priorities:
        try:
            f = float(str(v).strip())
        except (ValueError, TypeError):
            f = float("nan")
        pris.append(f)
    finite_pris = [p for p in pris if np.isfinite(p)]
    if not finite_pris:
        return None
    # Unique priority values, descending (lowest priority class first).
    classes = sorted(set(finite_pris), reverse=True)
    counts = {c: sum(1 for p in pris if p == c) for c in classes}
    weights: dict = {classes[0]: 1}
    w_prev = 1
    n_prev = counts[classes[0]]
    for q in range(1, len(classes)):
        c = classes[q]
        n_q = counts[c]
        # smallest int w satisfying both inequalities:
        # w > w_prev AND n_q * w > n_prev * w_prev
        candidate = max(w_prev + 1,
                        n_prev * w_prev // n_q + 1)
        weights[c] = candidate
        w_prev = candidate
        n_prev = n_q
    out: list[str] = []
    for p in pris:
        if not np.isfinite(p):
            out.append("")
        else:
            out.append(str(int(weights[p])))
    return out


def compute_priorities_from_weights(weights) -> list[str] | None:
    """Group rows by unique weight value, descending. Largest weight
    → priority 1; next → 2; …  Rows with non-finite weight get "" (NaN
    priority)."""
    ws: list[float] = []
    for v in weights:
        try:
            f = float(str(v).strip())
        except (ValueError, TypeError):
            f = float("nan")
        ws.append(f)
    finite_ws = [w for w in ws if np.isfinite(w)]
    if not finite_ws:
        return None
    classes = sorted(set(finite_ws), reverse=True)  # largest weight first
    p_of: dict = {c: i + 1 for i, c in enumerate(classes)}
    out: list[str] = []
    for w in ws:
        if not np.isfinite(w):
            out.append("")
        else:
            out.append(str(int(p_of[w])))
    return out
