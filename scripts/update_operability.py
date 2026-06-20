#!/usr/bin/env python
"""Fetch the *current* NIRSpec MSA operability reference into the CRDS cache.

vMPT reads the MSA failed-/stuck-shutter map from the MOS ``msaoper``
reference in your local CRDS cache (``~/crds_cache/references/jwst/nirspec/``).
As of v1.6.x vMPT already refreshes this automatically on startup
(best-effort; see ``vmpt.msa.ensure_current_operability``), but this script
forces the check on demand — handy after a long offline stretch, or to see
which reference is current.

It downloads the reference the *operational* CRDS context selects for MOS
(``EXP_TYPE = NRS_MSASPEC``, USEAFTER = today) — the same operability
APT/MPT uses — so vMPT matches it after a restart.

Usage:
    python scripts/update_operability.py

Needs network access to https://jwst-crds.stsci.edu. Safe to re-run.
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def main() -> int:
    os.environ.setdefault("CRDS_PATH", os.path.expanduser("~/crds_cache"))
    cache_dir = os.path.join(
        os.environ["CRDS_PATH"], "references", "jwst", "nirspec")
    before = {os.path.basename(p) for p in
              glob.glob(os.path.join(cache_dir, "jwst_nirspec_msaoper_*.json"))}
    print(f"CRDS cache: {cache_dir}")
    print(f"  msaoper files before: {sorted(before) or '(none)'}")

    try:
        from vmpt.msa import ensure_current_operability
    except ImportError as e:
        print(f"ERROR: cannot import vmpt ({e}). Run from the repo with the "
              f"stenv environment active.", file=sys.stderr)
        return 2

    # force=True (ignore the once-per-process guard / disable env var),
    # no timeout (this is a deliberate, blocking CLI), verbose output.
    path = ensure_current_operability(timeout=None, force=True, verbose=True)
    if not path:
        print("Could not resolve the current operability reference "
              "(offline, or crds unavailable).", file=sys.stderr)
        return 1

    name = os.path.basename(path)
    if name in before:
        print("  → already in your cache; you are up to date.")
    else:
        print("  → downloaded. Restart vMPT to load it "
              "(operability is read at startup).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
