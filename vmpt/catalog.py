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
    # ---- Per-target spectral constraints (v1.3.0+) ------------------
    # Each row optionally constrains how its spectrum must fall on the
    # detector given the current (disperser, filter). At every
    # candidate pointing the optimizer fetches the source's centre-
    # shutter wavelength endpoints via `vmpt.wavelengths.cutoffs` and
    # drops the target if any constraint fails.
    #
    # required_lam[i] is a list of (lam_lo, lam_hi) tuples in μm — the
    # **spectral coverage** the user requires for source i. Empty list
    # = no requirement. Stored as a dtype=object array so the ragged
    # per-row list lengths work in NumPy.
    required_lam: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=object)
    )
    # no_gap[i] = True → the NRS1/NRS2 detector gap must NOT fall
    # inside [lam_blue, lam_red] (i.e. cutoffs() must return NaN for
    # both gap_lo and gap_hi). Strict interpretation per v1.3.0.
    no_gap: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=bool)
    )
    # extend_blue[i] = True → the centre shutter's lam_blue must
    # reach the disperser/filter's nominal blue limit (no left-edge
    # truncation due to where in V2 the target sits).
    extend_blue: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=bool)
    )
    # extend_red[i] = True → same, for the red end.
    extend_red: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=bool)
    )
    # protect[i] = True → this target's spectrum is collision-
    # protected (same semantics as the v1.2 catalog-wide protect_mask
    # — the optimizer takes the logical OR of this flag and the v1.2
    # cutoff-derived mask).
    protect: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=bool)
    )
    # centration[i] is a per-target source-centering override. Empty
    # string ("") means "use the optimizer's global Source centering
    # setting" — the v1.3.0 default. Otherwise must be one of the
    # five canonical labels (UNCONSTRAINED, ENTIRE_OPEN, MIDPOINT,
    # CONSTRAINED, TIGHTLY_CONSTRAINED) — anything else gets coerced
    # to "" at load time. The optimizer reads this verbatim, looks up
    # CENTRATION_BUFFERS, and the per-target buffer wins **uncondition-
    # ally** — even when it's laxer than the global. v1.3.1+.
    centration: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=object)
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
# ---- v1.3.0 per-target spectral-constraint columns -------------------
# Loose-matched the same way as the canonical fields above. `lam_req`
# stores a string per row ("1.0-1.3; 1.5-1.8") parsed at load via
# `_parse_lam_req_str` into a list[(float, float)]. The four boolean
# columns accept any of the truthy-ish text values _BOOL_TRUE_TOKENS.
_LAM_REQ_KEYS = (
    "lamreq", "lambdareq", "lambdarequired", "wavelengthrequired",
    "requiredlam", "reqlam", "requiredwavelength",
)
_NO_GAP_KEYS = ("nogap", "gapless", "nodetectorgap")
_EXT_BLUE_KEYS = ("extendblue", "extendsblue", "blueextends", "bluest")
_EXT_RED_KEYS = ("extendred", "extendsred", "redextends", "reddest")
_PROTECT_KEYS = (
    "protect", "protected", "protectcollision", "collisionprotect",
)
# Per-target source-centering override (v1.3.1+). Cell value is one of
# `_VALID_CENTRATION_LEVELS` (case-insensitive); anything else becomes
# "" (no override). The match is loose — see `_as_centration_str` for
# the alias rules (e.g. "tight" → "TIGHTLY_CONSTRAINED").
_CENTRATION_KEYS = (
    "centration", "centering", "sourcecentration", "sourcecentering",
    "centerclass", "centeringclass",
)
_VALID_CENTRATION_LEVELS = (
    "UNCONSTRAINED",
    "ENTIRE_OPEN",
    "MIDPOINT",
    "CONSTRAINED",
    "TIGHTLY_CONSTRAINED",
)
# Short aliases users are likely to type, mapping back to the canonical
# label. The normalised key is lowercase + non-alnum stripped, matching
# `_norm()`'s shape (so "tightly-constrained" works too).
_CENTRATION_ALIASES = {
    "unconstrained": "UNCONSTRAINED", "unc": "UNCONSTRAINED",
    "none": "UNCONSTRAINED", "off": "UNCONSTRAINED",
    "entireopen": "ENTIRE_OPEN", "entire": "ENTIRE_OPEN",
    "open": "ENTIRE_OPEN",
    "midpoint": "MIDPOINT", "mid": "MIDPOINT", "middle": "MIDPOINT",
    "constrained": "CONSTRAINED", "con": "CONSTRAINED",
    "tight": "TIGHTLY_CONSTRAINED",
    "tightly": "TIGHTLY_CONSTRAINED",
    "tightlyconstrained": "TIGHTLY_CONSTRAINED",
    "tightconstrained": "TIGHTLY_CONSTRAINED",
}
_BOOL_TRUE_TOKENS = frozenset(
    ("1", "true", "yes", "y", "t", "✓", "✔", "on")
)

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


def _as_bool(table: Table, name: str | None) -> np.ndarray:
    """Coerce a column to bool, recognising any value in
    :data:`_BOOL_TRUE_TOKENS` as True. Empty / NaN / 0 / "false" are
    False. Used to read the per-target boolean constraint columns."""
    n = len(table)
    if name is None:
        return np.zeros(n, dtype=bool)
    col = table[name]
    out = np.zeros(n, dtype=bool)
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
        s = str(v).strip().lower()
        if not s or s in ("nan", "none", "null", "--"):
            continue
        out[i] = s in _BOOL_TRUE_TOKENS
    return out


def _parse_lam_req_str(s: str) -> list[tuple[float, float]]:
    """Parse the user-facing wavelength-range string format.

    Format: zero or more ``"lo-hi"`` ranges in μm, semicolon- or
    comma-separated. Examples:

      ""                  → []
      "1.0-1.3"           → [(1.0, 1.3)]
      "1.0-1.3; 1.5-1.8"  → [(1.0, 1.3), (1.5, 1.8)]
      "0.9 - 1.0, 2 - 3"  → [(0.9, 1.0), (2.0, 3.0)]

    Invalid fragments are silently dropped — the popover UI shows a
    yellow warning on save when it spots them; from the loader's
    perspective they just become missing constraints.
    """
    if s is None:
        return []
    s = str(s).strip()
    if not s or s.lower() in ("nan", "none", "null", "--"):
        return []
    out: list[tuple[float, float]] = []
    for chunk in re.split(r"[;,]", s):
        c = chunk.strip()
        if not c:
            continue
        # "1.0-1.3" or "1.0 — 1.3" or "1.0 to 1.3"
        m = re.match(
            r"^\s*([\d.eE+-]+)\s*(?:-|–|—|to)\s*([\d.eE+-]+)\s*$",
            c,
        )
        if not m:
            continue
        try:
            lo, hi = float(m.group(1)), float(m.group(2))
        except ValueError:
            continue
        if lo > hi:
            lo, hi = hi, lo
        if np.isfinite(lo) and np.isfinite(hi):
            out.append((lo, hi))
    return out


def _format_lam_req(ranges) -> str:
    """Inverse of :func:`_parse_lam_req_str` — serialise back to the
    string format used in the CSV and in the popover input."""
    if ranges is None:
        return ""
    parts: list[str] = []
    for r in ranges:
        try:
            lo, hi = float(r[0]), float(r[1])
        except (TypeError, ValueError, IndexError):
            continue
        if np.isfinite(lo) and np.isfinite(hi):
            parts.append(f"{lo:g}-{hi:g}")
    return "; ".join(parts)


def _as_lam_req(table: Table, name: str | None) -> np.ndarray:
    """Read a wavelength-required column from the catalog. Each cell is
    parsed by :func:`_parse_lam_req_str` into a `list[tuple]`; the
    result is wrapped in a `dtype=object` array (ragged-friendly)."""
    n = len(table)
    if name is None:
        return np.array([[] for _ in range(n)], dtype=object)
    col = table[name]
    out = np.empty(n, dtype=object)
    mask = getattr(col, "mask", None)
    for i, v in enumerate(col):
        masked = False
        if mask is not None and mask is not False:
            try:
                masked = bool(mask[i])
            except (TypeError, IndexError):
                masked = False
        out[i] = [] if masked else _parse_lam_req_str(v)
    return out


def _normalise_centration(value) -> str:
    """Coerce a single centration cell to one of the canonical labels.

    Returns ``""`` (no override) for empty / unrecognised values, or
    one of the strings in :data:`_VALID_CENTRATION_LEVELS`. Matching
    is case-insensitive and tolerant of underscores / hyphens (e.g.
    ``"tightly-constrained"`` and ``"Tight"`` both resolve to
    ``"TIGHTLY_CONSTRAINED"``).
    """
    if value is None:
        return ""
    s = str(value).strip()
    if not s or s.lower() in ("nan", "none", "null", "--"):
        return ""
    # Direct canonical match (handles the common case fast).
    upper = s.upper().replace("-", "_").replace(" ", "_")
    if upper in _VALID_CENTRATION_LEVELS:
        return upper
    # Loose-match against the alias table (lowercase + alnum-only).
    key = re.sub(r"[^a-z0-9]+", "", s.lower())
    if key in _CENTRATION_ALIASES:
        return _CENTRATION_ALIASES[key]
    return ""


def _as_centration_str(table: Table, name: str | None) -> np.ndarray:
    """Read a centration-override column. Cells coerce via
    :func:`_normalise_centration`; unrecognised values silently become
    ``""`` (i.e. "use the optimizer's global setting")."""
    n = len(table)
    if name is None:
        return np.array([""] * n, dtype=object)
    col = table[name]
    out = np.empty(n, dtype=object)
    mask = getattr(col, "mask", None)
    for i, v in enumerate(col):
        masked = False
        if mask is not None and mask is not False:
            try:
                masked = bool(mask[i])
            except (TypeError, IndexError):
                masked = False
        out[i] = "" if masked else _normalise_centration(v)
    return out


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
    # Per-target constraint columns (v1.3.0+). All optional; defaults
    # leave constraints unset and v1.2.x behaviour is preserved.
    lam_req_col = _find_col(table, _LAM_REQ_KEYS)
    no_gap_col = _find_col(table, _NO_GAP_KEYS)
    extend_blue_col = _find_col(table, _EXT_BLUE_KEYS)
    extend_red_col = _find_col(table, _EXT_RED_KEYS)
    protect_col = _find_col(table, _PROTECT_KEYS)
    # Per-target source-centering override (v1.3.1+).
    centration_col = _find_col(table, _CENTRATION_KEYS)
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
                           mag_col, z_col, label_col,
                           lam_req_col, no_gap_col, extend_blue_col,
                           extend_red_col, protect_col,
                           centration_col) if c}
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
        required_lam=_as_lam_req(table, lam_req_col),
        no_gap=_as_bool(table, no_gap_col),
        extend_blue=_as_bool(table, extend_blue_col),
        extend_red=_as_bool(table, extend_red_col),
        protect=_as_bool(table, protect_col),
        centration=_as_centration_str(table, centration_col),
        source_path=path,
        extras=extras,
    )


def save_catalog(cat: Catalog, path: str, *,
                 include_constraints: str = "auto") -> None:
    """Write a :class:`Catalog` back to CSV.

    Emits the standard eight columns (``ID, RA, DEC, priority,
    weight, mag, z, label``) followed optionally by the six v1.3.x
    per-target constraint columns (``lam_req, no_gap, extend_blue,
    extend_red, protect, centration``), then any ``extras`` columns
    the catalog is carrying. The output is round-trip-compatible
    with :func:`load_catalog` — write, reload, and the resulting
    :class:`Catalog` matches the input modulo dtype.

    Parameters
    ----------
    cat : Catalog
        The catalog to save.
    path : str
        Destination CSV path. Parent directories are NOT created
        automatically — callers should ensure the parent exists.
    include_constraints : {"auto", "always", "never"}
        Controls whether the six constraint columns appear in the
        output:

        - ``"auto"`` (default): emit the columns iff at least one
          row has a non-default value (the same rule the catalog
          editor's Save-as-CSV button uses, so v1.2.x catalogs that
          never picked up constraints get the same CSV format they
          had before).
        - ``"always"``: always emit the columns, even when every
          row is at defaults. Useful when you want a "template"
          CSV the user can hand-edit.
        - ``"never"``: omit them. Use when you specifically want
          to drop the constraint metadata.

    Notes
    -----
    Wavelength-range cells round-trip via the same string format
    the catalog editor uses (``"1.0-1.3; 1.5-1.8"``), parsed back
    by :func:`_parse_lam_req_str` on the next load.

    NaN / missing values render as empty cells in the CSV; on
    reload they come back as NaN (for float columns) or empty
    string (for label / extras).
    """
    import csv

    if include_constraints not in ("auto", "always", "never"):
        raise ValueError(
            f"include_constraints must be one of "
            f"'auto', 'always', 'never'; got {include_constraints!r}"
        )

    n = len(cat.ra_deg)

    def _fmt_int_or_blank(v) -> str:
        try:
            f = float(v)
            if not np.isfinite(f):
                return ""
            return str(int(round(f)))
        except (TypeError, ValueError):
            return "" if v is None else str(v)

    def _fmt_float_or_blank(v) -> str:
        try:
            f = float(v)
            if not np.isfinite(f):
                return ""
            # Drop trailing zeros / trailing decimal so the CSV is
            # tidy for hand-editing.
            s = f"{f:.6f}".rstrip("0").rstrip(".")
            return s or "0"
        except (TypeError, ValueError):
            return "" if v is None else str(v)

    def _fmt_id(v) -> str:
        # Catalog.ids can be int64 or object dtype (string IDs).
        try:
            return str(int(v))
        except (TypeError, ValueError):
            return "" if v is None else str(v)

    # Decide constraint-column emission policy. `"auto"` mode emits
    # only when at least one row has a non-default value.
    required_lam = getattr(cat, "required_lam", None)
    no_gap = np.asarray(getattr(cat, "no_gap", []), dtype=bool)
    extend_blue = np.asarray(getattr(cat, "extend_blue", []), dtype=bool)
    extend_red = np.asarray(getattr(cat, "extend_red", []), dtype=bool)
    protect = np.asarray(getattr(cat, "protect", []), dtype=bool)
    centration = np.asarray(getattr(cat, "centration", []), dtype=object)
    has_lam = (required_lam is not None
               and len(required_lam) == n
               and any(bool(r) and len(r) > 0 for r in required_lam))
    has_centration = (
        centration.size == n
        and any(bool(str(v).strip()) for v in centration)
    )
    has_constraints = (
        has_lam
        or (no_gap.size == n and no_gap.any())
        or (extend_blue.size == n and extend_blue.any())
        or (extend_red.size == n and extend_red.any())
        or (protect.size == n and protect.any())
        or has_centration
    )
    emit_constraints = (
        include_constraints == "always"
        or (include_constraints == "auto" and has_constraints)
    )

    extras = getattr(cat, "extras", {}) or {}

    constraint_cols = (["lam_req", "no_gap", "extend_blue",
                        "extend_red", "protect", "centration"]
                       if emit_constraints else [])
    header = (["ID", "RA", "DEC", "priority", "weight", "mag", "z",
               "label", *constraint_cols, *extras.keys()])

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i in range(n):
            row = [
                _fmt_id(cat.ids[i]),
                _fmt_float_or_blank(cat.ra_deg[i]),
                _fmt_float_or_blank(cat.dec_deg[i]),
                _fmt_int_or_blank(cat.priority[i])
                    if i < len(cat.priority) else "",
                _fmt_int_or_blank(cat.weight[i])
                    if i < len(cat.weight) else "",
                _fmt_float_or_blank(cat.mag[i])
                    if i < len(cat.mag) else "",
                _fmt_float_or_blank(cat.z[i])
                    if i < len(cat.z) else "",
                str(cat.label[i]) if i < len(cat.label) else "",
            ]
            if emit_constraints:
                rl = (required_lam[i] if required_lam is not None
                      and i < len(required_lam) else [])
                row.append(_format_lam_req(rl) if rl else "")
                row.append("1" if (no_gap.size > i
                                   and bool(no_gap[i])) else "")
                row.append("1" if (extend_blue.size > i
                                   and bool(extend_blue[i])) else "")
                row.append("1" if (extend_red.size > i
                                   and bool(extend_red[i])) else "")
                row.append("1" if (protect.size > i
                                   and bool(protect[i])) else "")
                # Empty cell when no override; otherwise the canonical
                # label (already normalised by the loader or `_normalise_
                # centration`). Writing an unrecognised label here is
                # not a hard error — the next load just resets it to "".
                if centration.size > i:
                    c = str(centration[i]).strip()
                    row.append(c if c else "")
                else:
                    row.append("")
            for k, vals in extras.items():
                row.append(str(vals[i]) if i < len(vals) else "")
            w.writerow(row)


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
