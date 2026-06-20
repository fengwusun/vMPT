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
    python scripts/update_operability.py            # refresh your CRDS cache
    python scripts/update_operability.py --bundle   # ALSO refresh the shipped
                                                     # fallback snapshot (maintainers)

`--bundle` regenerates ``vmpt/data/msaoper_fallback.npz`` — the compact
operability snapshot vMPT ships and loads when no CRDS reference is available
(fresh installs with no ``crds``/``CRDS_PATH``/network). Run it before cutting
a release so the bundled floor stays current.

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

    if "--bundle" in sys.argv:
        _write_bundle(path)
    return 0


def _write_bundle(json_path: str) -> None:
    """Regenerate the shipped fallback snapshot from a msaoper JSON."""
    import json

    import numpy as np

    reason = np.zeros((4, 171, 365), dtype=np.int8)
    with open(json_path) as f:
        data = json.load(f)
    for e in data.get("msaoper", []):
        q = int(e["Q"]); d = int(e["x"]); s = int(e["y"])
        if not (1 <= q <= 4 and 1 <= s <= 171 and 1 <= d <= 365):
            continue
        st = str(e.get("state", "")).lower()
        if "open" in st:
            reason[q - 1, s - 1, d - 1] = 2
        elif "closed" in st:
            reason[q - 1, s - 1, d - 1] = 1
    out = _ROOT / "vmpt" / "data" / "msaoper_fallback.npz"
    np.savez_compressed(out, reason=reason, source=os.path.basename(json_path))
    print(f"  → regenerated bundled fallback {out} "
          f"(stuck-open={int((reason == 2).sum())}, "
          f"failed-closed={int((reason == 1).sum())}).")


if __name__ == "__main__":
    raise SystemExit(main())
