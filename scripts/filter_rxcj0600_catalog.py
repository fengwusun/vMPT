"""Slim down the RXCJ0600 example catalog for public distribution.

The original ``example_r0600/v01_fsun.cat`` ships with ~28 600 rows.
Loading that into the Bokeh canvas is slow on a fresh launch — most
users just want a representative scatter to try the picker against.

This script applies two filters and rewrites the file in place:

  1. F200W_mag or F444W_mag inside the JWST/NIRSpec sweet spot
     (24 ≤ mag ≤ 27 in either band — the OR, so dropout sources
     in one filter still make it through).
  2. Random subsample to **2000** rows (`numpy.random.default_rng(42)`
     for reproducibility) — keeps a representative magnitude /
     spatial distribution without biasing toward the brightest or
     faintest tail.

Usage::

    python scripts/filter_rxcj0600_catalog.py
        [--source PATH] [--dest PATH] [--cap N] [--mag-lo M] [--mag-hi M]

Defaults reproduce the public-release filter.  Re-run with
``--source path/to/full_catalog.cat`` to regenerate if you have the
original 28 k catalog stashed elsewhere; the in-repo example file is
the FILTERED one, so re-running on it just leaves it unchanged.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_PATH = _REPO_ROOT / "example_r0600" / "v01_fsun.cat"


def _parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Slim down the RXCJ0600 example catalog "
                    "for distribution."
    )
    p.add_argument(
        "--source", type=Path, default=_DEFAULT_PATH,
        help="Path to the full source catalog (default: in-repo file).",
    )
    p.add_argument(
        "--dest", type=Path, default=_DEFAULT_PATH,
        help="Where to write the filtered catalog. Default = source "
             "(in-place rewrite).",
    )
    p.add_argument("--mag-lo", type=float, default=24.0,
                   help="Lower magnitude bound (default: 24).")
    p.add_argument("--mag-hi", type=float, default=27.0,
                   help="Upper magnitude bound (default: 27).")
    p.add_argument("--cap", type=int, default=2000,
                   help="Maximum row count (default: 2000).")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for the subsample (default: 42).")
    return p.parse_args(argv)


def filter_catalog(
    df: pd.DataFrame,
    *,
    mag_lo: float, mag_hi: float, cap: int, seed: int,
) -> pd.DataFrame:
    """Return a slimmed-down dataframe.

    Rows kept iff ``F200W_mag`` or ``F444W_mag`` is in
    ``[mag_lo, mag_hi]`` (NaN treated as "not in range"). The
    survivors get subsampled to at most ``cap`` rows via a
    deterministic ``numpy.random.default_rng(seed)``.
    """
    f200 = pd.to_numeric(df["F200W_mag"], errors="coerce")
    f444 = pd.to_numeric(df["F444W_mag"], errors="coerce")
    in_200 = (f200 >= mag_lo) & (f200 <= mag_hi)
    in_444 = (f444 >= mag_lo) & (f444 <= mag_hi)
    keep = (in_200 | in_444).fillna(False)
    filtered = df[keep].copy()
    n_kept = len(filtered)

    if n_kept > cap:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_kept, size=cap, replace=False)
        # Sort the sampled indices so the output preserves the
        # original catalog's RA / Dec / ID ordering (easier diffing
        # when the script is re-run with different params).
        idx.sort()
        filtered = filtered.iloc[idx].reset_index(drop=True)
    else:
        filtered = filtered.reset_index(drop=True)
    return filtered


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    if not args.source.exists():
        sys.stderr.write(f"error: source catalog not found: {args.source}\n")
        return 2

    df = pd.read_csv(args.source)
    n0 = len(df)
    filtered = filter_catalog(
        df,
        mag_lo=args.mag_lo, mag_hi=args.mag_hi,
        cap=args.cap, seed=args.seed,
    )
    n1 = len(filtered)

    args.dest.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.dest, index=False)
    print(
        f"vmpt: filtered {args.source.name} "
        f"{n0:,} rows → {n1:,} rows "
        f"(F200W or F444W in [{args.mag_lo:g}, {args.mag_hi:g}], "
        f"cap={args.cap}, seed={args.seed})"
    )
    print(f"  wrote {args.dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
