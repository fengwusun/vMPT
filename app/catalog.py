"""Target catalog loader (CSV / ASCII / FITS table)."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass

import numpy as np
from astropy.io import ascii as ioascii
from astropy.table import Table

# Pattern that picks up the *numeric portion* of a value like "P0",
# "P1", "class-3", "1.5e-2". Used by _as_float so common JWST priority-
# class encodings (P0 = highest, P1 = …) flow through as numeric 0, 1,
# etc. without forcing the user to hand-edit their catalog.
_NUM_RE = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


@dataclass
class Catalog:
    ids: np.ndarray
    ra_deg: np.ndarray
    dec_deg: np.ndarray
    priority: np.ndarray
    mag: np.ndarray
    z: np.ndarray
    label: np.ndarray
    source_path: str


_ID_KEYS = ("id", "no", "no_cat", "objid", "objectid", "source_id")
_RA_KEYS = ("ra", "ra_deg", "ra[deg]", "raj2000", "alpha_j2000")
_DEC_KEYS = ("dec", "dec_deg", "dec[deg]", "decj2000", "delta_j2000")
_PRI_KEYS = ("priority", "pr", "pri", "prio")
_MAG_KEYS = ("mag", "magnitude", "f444w_mag", "mag_f444w", "f356w_mag", "mag_f356w", "f200w_mag", "mag_f200w")
_Z_KEYS = ("z", "zspec", "zphot", "z_spec", "z_phot", "redshift", "z_best", "z_use")
_LABEL_KEYS = ("label", "name", "tag")


def _find_col(table: Table, candidates) -> str | None:
    lc_map = {c.lower(): c for c in table.colnames}
    for cand in candidates:
        if cand in lc_map:
            return lc_map[cand]
    return None


def _as_float(table: Table, name: str | None) -> np.ndarray:
    """Coerce a column to float, tolerantly.

    Catalogs in the wild use a few non-numeric conventions for fields
    that vMPT wants as numbers — the most common is the **priority
    class** (`P0`, `P1`, …). Rather than throwing, we:

      • try the fast path (`np.asarray(..., dtype=float)`);
      • on failure fall back to row-by-row parsing — empty strings and
        masked values become NaN, and the *numeric portion* of any
        string is extracted (so `"P0"` → 0.0, `"class-3"` → 3.0,
        `"high-mu"` → NaN).
    """
    n = len(table)
    if name is None:
        return np.full(n, np.nan, dtype=float)
    col = table[name]
    # Numeric column with astropy masks → fill masked entries with NaN.
    # (np.asarray on a MaskedArray drops the mask and exposes the
    # underlying buffer, which usually has 0 in the masked slots — not
    # what we want for empty `mag` / `z` cells.)
    if np.issubdtype(getattr(col, "dtype", np.dtype("O")), np.floating):
        try:
            return np.ma.filled(np.ma.asarray(col), np.nan).astype(float)
        except (ValueError, TypeError):
            pass
    # Numeric (int) column → straight cast is fine.
    try:
        return np.asarray(col, dtype=float)
    except (ValueError, TypeError):
        pass
    # Non-numeric column → row-by-row parse, extracting trailing digits.
    out = np.full(n, np.nan, dtype=float)
    mask = getattr(col, "mask", None)
    for i, v in enumerate(col):
        if mask is not None and mask is not False:
            try:
                if mask[i]:
                    continue
            except (TypeError, IndexError):
                pass
        if v is None:
            continue
        s = str(v).strip()
        if not s or s.lower() in ("--", "nan", "none", "null"):
            continue
        try:
            out[i] = float(s)
            continue
        except ValueError:
            pass
        m = _NUM_RE.search(s)
        if m is not None:
            try:
                out[i] = float(m.group(0))
            except ValueError:
                pass
    return out


def _as_str(table: Table, name: str | None) -> np.ndarray:
    n = len(table)
    if name is None:
        return np.array([""] * n, dtype=object)
    return np.asarray([str(v) for v in table[name]], dtype=object)


def load_catalog(path: str) -> Catalog:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".fits", ".fit", ".fz"):
        table = Table.read(path)
    elif ext == ".csv":
        table = ioascii.read(path, format="csv")
    else:
        table = ioascii.read(path)

    id_col = _find_col(table, _ID_KEYS)
    ra_col = _find_col(table, _RA_KEYS)
    dec_col = _find_col(table, _DEC_KEYS)
    if ra_col is None or dec_col is None:
        raise ValueError(f"Catalog at {path} missing RA/Dec columns. Have: {table.colnames}")
    if id_col is None:
        ids = np.arange(1, len(table) + 1)
    else:
        raw = table[id_col]
        try:
            ids = np.asarray(raw, dtype=np.int64)
        except (ValueError, TypeError):
            ids = np.asarray([str(v) for v in raw], dtype=object)

    pri_col = _find_col(table, _PRI_KEYS)
    mag_col = _find_col(table, _MAG_KEYS)
    z_col = _find_col(table, _Z_KEYS)
    label_col = _find_col(table, _LABEL_KEYS)

    return Catalog(
        ids=ids,
        ra_deg=np.asarray(table[ra_col], dtype=float),
        dec_deg=np.asarray(table[dec_col], dtype=float),
        priority=_as_float(table, pri_col),
        mag=_as_float(table, mag_col),
        z=_as_float(table, z_col),
        label=_as_str(table, label_col),
        source_path=path,
    )


def catalog_in_view(cat: Catalog, ra_min, ra_max, dec_min, dec_max) -> np.ndarray:
    ra = cat.ra_deg
    dec = cat.dec_deg
    in_dec = (dec >= dec_min) & (dec <= dec_max)
    if ra_min <= ra_max:
        in_ra = (ra >= ra_min) & (ra <= ra_max)
    else:
        # RA range wraps across 0/360
        in_ra = (ra >= ra_min) | (ra <= ra_max)
    return in_ra & in_dec
