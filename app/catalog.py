"""Target catalog loader (CSV / ASCII / FITS table)."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np
from astropy.io import ascii as ioascii
from astropy.table import Table

# Pattern that picks up the *numeric portion* of a value like "P0",
# "P1", "class-3", "1.5e-2". Used by _as_float so common JWST priority-
# class encodings (P0 = highest, P1 = …) flow through as numeric 0, 1,
# etc. without forcing the user to hand-edit their catalog.
_NUM_RE = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")

# Catalog IDs above this threshold are taken mod ID_MOD before being
# stored. JADES-style IDs can run to 8–9 digits, but APT MPT and the
# eMPT pipeline both expect compact integer source numbers — anything
# beyond ~10⁷ tends to be silently truncated or rejected downstream.
# Collisions after the mod are vanishingly rare in real catalogs;
# we accept that trade-off in exchange for a clean integer space.
ID_MOD = 10_000_000


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
    # `weight` is a sibling of `priority`: float64 with NaN for missing.
    # Used by the optimizer's Meritocracy mode (sum of placed weights),
    # and by the Hierarchy mode internally to break ties within a
    # priority tier.
    weight: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=float)
    )
    # Original-column → values for every column the loader did NOT
    # claim as one of the canonical fields above. Stored as object
    # arrays so we don't lose mixed-type information (e.g., priority-
    # class strings, free-text notes). The catalog editor exposes
    # these via the column picker; we never mutate them algorithmically
    # so they round-trip back into Save-as-CSV unchanged unless the
    # user edits them.
    extras: dict = field(default_factory=dict)


# Lookup tables for the loose column-matcher (`_find_col`). Each
# candidate is normalised with `_norm` (lowercase + strip bracketed
# units + collapse to alphanumeric + strip trailing unit tokens). The
# normalisation makes `RA`, `ra`, `RA[deg]`, `RA(deg)`, `RA_deg`,
# `ra_J2000`, `ALPHA_J2000`, `R.A.` all map to the same key.
_ID_KEYS = (
    "id", "no", "nocat", "objid", "objectid", "sourceid", "source",
    "src", "srcid", "targetid", "targid", "ident",
)
# Permissive ID fallbacks: accepted only when the column's values are
# numeric (else we'd silently sort sources by their human-readable name).
_ID_FALLBACK_KEYS = ("name", "label", "tag", "target", "targetname", "#")
_RA_KEYS = (
    "ra", "rightascension", "raj2000", "alpha", "alphaj2000",
    "rad", "radeg",
)
_DEC_KEYS = (
    "dec", "declination", "decj2000", "delta", "deltaj2000",
    "decd", "decdeg",
    # Vizier-style "DEJ2000" normalises to "de" once "J2000" is
    # stripped as a unit/epoch token. Adding "de" keeps that catalog
    # convention working.
    "de",
)
_PRI_KEYS = ("priority", "pr", "pri", "prio", "priorityclass")
_WEIGHT_KEYS = ("weight", "w", "wt", "weights")
_MAG_KEYS = (
    "mag", "magnitude", "f444wmag", "magf444w", "f356wmag", "magf356w",
    "f200wmag", "magf200w",
)
_Z_KEYS = ("z", "zspec", "zphot", "redshift", "zbest", "zuse")
_LABEL_KEYS = ("label", "name", "tag")

# Trailing tokens that look like *units* on an otherwise-clean column
# name — stripped after lowercasing + alphanumeric collapse so
# `RA[deg]`, `RA_deg`, `ra (deg)`, `RAJ2000` all collapse to `ra`.
_UNIT_SUFFIX_TOKENS = (
    "degrees", "degree", "deg",
    "radians", "radian", "rad",
    "arcseconds", "arcsec", "asec",
    "j2000", "icrs", "fk5",
)


def _norm(name: str) -> str:
    """Normalise a column name for loose matching."""
    if name is None:
        return ""
    s = str(name).lower()
    # Strip bracketed / parenthesised unit suffixes ("RA[deg]" → "RA").
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^)]*\)", "", s)
    # Collapse remaining non-alphanumerics ("ra_deg" → "radeg", "R.A." → "ra").
    s = re.sub(r"[^a-z0-9]+", "", s)
    # Strip trailing unit tokens ("radeg" → "ra"). Loop so chained
    # suffixes (e.g. "decjsiomdeg") peel off one by one.
    changed = True
    while changed:
        changed = False
        for tok in _UNIT_SUFFIX_TOKENS:
            if len(s) > len(tok) and s.endswith(tok):
                s = s[: -len(tok)]
                changed = True
                break
    return s


def _find_col(table: Table, candidates) -> str | None:
    """Return the original column name matching any normalised candidate."""
    norm_map: dict[str, str] = {}
    for c in table.colnames:
        norm_map.setdefault(_norm(c), c)
    for cand in candidates:
        n = _norm(cand)
        if n and n in norm_map:
            return norm_map[n]
    return None


def _find_id_col(table: Table) -> tuple[str | None, bool]:
    """Locate the catalog's ID column.

    Returns `(name, is_numeric_fallback)`. The fallback flag is True
    when we accepted a permissive candidate (`name`, `label`, …)
    *because* its values coerced to integers — used downstream to
    decide whether to preserve the original token alongside the int ID.
    """
    name = _find_col(table, _ID_KEYS)
    if name is not None:
        return name, False
    # Permissive: accept name/label/tag only if values look like integers.
    for cand in _ID_FALLBACK_KEYS:
        col_name = _find_col(table, (cand,))
        if col_name is None:
            continue
        col = table[col_name]
        try:
            arr = np.asarray(col, dtype=np.int64)
        except (ValueError, TypeError):
            continue
        # Sanity: empty / all-zero columns are unlikely to be IDs.
        if arr.size > 0:
            return col_name, True
    return None, False


def _coerce_int_ids(raw, nrows: int) -> np.ndarray:
    """Return an int64 ID array of length `nrows`, with mod ID_MOD
    applied to any source ID at or above 10⁷.

    If `raw` can't be coerced to int (string IDs like "RJ0600-x-P0"),
    we return the raw values as an object array — the integer
    extraction happens later in the exporter's `_to_int_id`."""
    try:
        ids = np.asarray(raw, dtype=np.int64)
    except (ValueError, TypeError):
        return np.asarray([str(v) for v in raw], dtype=object)
    big = np.abs(ids) >= ID_MOD
    if big.any():
        ids = ids.copy()
        ids[big] = np.mod(ids[big], ID_MOD)
    return ids


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
    # what we want for empty `mag` / `z` cells.) Apply this for any
    # numeric dtype, not just floats — empty cells in an int column
    # also need NaN handling.
    if np.issubdtype(getattr(col, "dtype", np.dtype("O")), np.number):
        try:
            arr = np.ma.asarray(col)
            # Cast to float FIRST so the NaN fill is representable,
            # then fill. Otherwise filling an int-typed masked array
            # with `np.nan` silently coerces to 0.
            return np.ma.filled(arr.astype(float), np.nan)
        except (ValueError, TypeError):
            pass
    # Non-numeric column — fall through to row-by-row parse below.
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

    ra_col = _find_col(table, _RA_KEYS)
    dec_col = _find_col(table, _DEC_KEYS)
    if ra_col is None or dec_col is None:
        raise ValueError(
            f"Catalog at {path} missing RA/Dec columns. Have: {table.colnames}"
        )

    id_col, id_from_fallback = _find_id_col(table)
    if id_col is None:
        # Catalog has no ID-like column — fake sequential IDs 1..N so
        # downstream code (slitlet auto-tag, MPT export) still works.
        ids = np.arange(1, len(table) + 1, dtype=np.int64)
    else:
        ids = _coerce_int_ids(table[id_col], len(table))

    pri_col = _find_col(table, _PRI_KEYS)
    weight_col = _find_col(table, _WEIGHT_KEYS)
    mag_col = _find_col(table, _MAG_KEYS)
    z_col = _find_col(table, _Z_KEYS)
    # If `name`/`label` was used as the ID fallback, don't ALSO claim it
    # as the label column — that would just duplicate the ID.
    label_candidates = _LABEL_KEYS
    if id_from_fallback and id_col is not None:
        label_candidates = tuple(
            k for k in _LABEL_KEYS if _norm(k) != _norm(id_col)
        )
    label_col = _find_col(table, label_candidates)

    # Every column the loader didn't claim above survives in `extras`
    # so the catalog editor can show it. Stored as object arrays so we
    # don't lose mixed-type information.
    claimed = {c for c in (id_col, ra_col, dec_col, pri_col, weight_col,
                           mag_col, z_col, label_col) if c}
    extras: dict = {}
    for col_name in table.colnames:
        if col_name in claimed:
            continue
        col = table[col_name]
        try:
            extras[col_name] = np.asarray([
                ("" if v is None else str(v))
                for v in (np.ma.filled(col, "") if hasattr(col, "mask") else col)
            ], dtype=object)
        except (TypeError, ValueError):
            # Fall back to a plain object copy; if even that fails we
            # silently drop the column rather than crashing the loader.
            try:
                extras[col_name] = np.asarray(list(col), dtype=object)
            except Exception:  # noqa: BLE001
                pass

    return Catalog(
        ids=ids,
        ra_deg=np.asarray(table[ra_col], dtype=float),
        dec_deg=np.asarray(table[dec_col], dtype=float),
        priority=_as_float(table, pri_col),
        weight=_as_float(table, weight_col),
        mag=_as_float(table, mag_col),
        z=_as_float(table, z_col),
        label=_as_str(table, label_col),
        source_path=path,
        extras=extras,
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
