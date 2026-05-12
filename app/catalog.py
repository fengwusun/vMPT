"""Target catalog loader (CSV / ASCII / FITS table)."""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
from astropy.io import ascii as ioascii
from astropy.table import Table


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
    n = len(table)
    if name is None:
        return np.full(n, np.nan, dtype=float)
    return np.asarray(table[name], dtype=float)


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
