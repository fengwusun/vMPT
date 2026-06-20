"""NIRSpec MSA shutter grid loader and operability parser."""

import datetime as _dt
import json
import os
import re
from glob import glob
from pathlib import Path

import numpy as np

_DATA_DIR = Path(__file__).resolve().parent / "data"
_GRID_PATH = _DATA_DIR / "nirspec_msa_v2v3.npz"

_grid_cache: dict = {}


def load_msa_grid() -> tuple[np.ndarray, np.ndarray]:
    if "v2" not in _grid_cache:
        d = np.load(_GRID_PATH)
        _grid_cache["v2"] = d["v2_msa"]
        _grid_cache["v3"] = d["v3_msa"]
    return _grid_cache["v2"], _grid_cache["v3"]


def shutter_center_v2v3(q: int, s: int, d: int) -> tuple[float, float]:
    v2, v3 = load_msa_grid()
    return float(v2[q - 1, s - 1, d - 1]), float(v3[q - 1, s - 1, d - 1])


def shutters_in_bbox(v2_min: float, v2_max: float, v3_min: float, v3_max: float) -> np.ndarray:
    v2, v3 = load_msa_grid()
    mask = (v2 >= v2_min) & (v2 <= v2_max) & (v3 >= v3_min) & (v3 <= v3_max)
    qs, ss, ds = np.where(mask)
    out = np.empty((qs.size, 3), dtype=np.int32)
    out[:, 0] = qs + 1
    out[:, 1] = ss + 1
    out[:, 2] = ds + 1
    return out


def _select_current_msaoper(ref_dir: str) -> str | None:
    """Pick the msaoper reference CRDS would use for **MOS** observations
    *now* — the ``NRS_MSASPEC`` branch of the newest cached msaoper rmap,
    latest USEAFTER ≤ today, among the reference files actually present in
    ``ref_dir``.

    CRDS reference *numbers* are delivery order, not operability date, and
    the rmap also carries an imaging-only default (``msaoper_0001``); so
    "newest file by name" can pick the wrong reference. Reading the rmap's
    UseAfter table matches how APT/MPT selects operability. Returns ``None``
    (and the caller falls back to the highest-numbered file) if the rmap
    isn't cached or can't be parsed.
    """
    try:
        present = {
            os.path.basename(p)
            for p in glob(os.path.join(ref_dir, "jwst_nirspec_msaoper_*.json"))
        }
        if not present:
            return None
        # references/jwst/nirspec → <CRDS>/mappings/jwst
        map_dir = os.path.normpath(
            os.path.join(ref_dir, "..", "..", "..", "mappings", "jwst"))
        rmaps = sorted(glob(os.path.join(map_dir, "jwst_nirspec_msaoper_*.rmap")))
        if not rmaps:
            return None
        with open(rmaps[-1]) as f:
            text = f.read()
        # Isolate the MOS UseAfter block: '…NRS_MSASPEC' : UseAfter({ … })
        block = re.search(
            r"MSASPEC['\"]?\s*:\s*UseAfter\(\{(.*?)\}\)", text, re.S)
        if not block:
            return None
        today = _dt.datetime.now()
        best: tuple[_dt.datetime, str] | None = None
        for date_str, fname in re.findall(
            r"'([0-9][0-9\-: ]+)'\s*:\s*'(jwst_nirspec_msaoper_\d+\.json)'",
            block.group(1),
        ):
            try:
                when = _dt.datetime.strptime(date_str.strip(), "%Y-%m-%d %H:%M:%S")
            except ValueError:
                continue
            if when <= today and fname in present:
                if best is None or when > best[0]:
                    best = (when, fname)
        if best is None:
            return None
        return os.path.join(ref_dir, best[1])
    except Exception:  # noqa: BLE001 — selection is best-effort; fall back
        return None


def _find_msaoper_json(crds_path: str | None) -> str | None:
    candidates: list[str] = []
    if crds_path:
        candidates.append(os.path.join(crds_path, "references", "jwst", "nirspec"))
    env_crds = os.environ.get("CRDS_PATH")
    if env_crds:
        candidates.append(os.path.join(env_crds, "references", "jwst", "nirspec"))
    candidates.append(os.path.expanduser("~/crds_cache/references/jwst/nirspec"))
    for cand in candidates:
        if not os.path.isdir(cand):
            continue
        # Prefer the rmap-driven MOS-current reference; fall back to the
        # highest-numbered file if the rmap isn't available.
        sel = _select_current_msaoper(cand)
        if sel:
            return sel
        files = sorted(glob(os.path.join(cand, "jwst_nirspec_msaoper_*.json")))
        if files:
            return files[-1]
    return None


_OPER_CHECK_DONE = False


def ensure_current_operability(timeout: float | None = 8.0,
                               force: bool = False,
                               verbose: bool = False) -> str | None:
    """Best-effort refresh of the MSA operability reference at startup.

    Make sure the *current* MOS ``msaoper`` reference — the one the
    operational CRDS context selects for ``EXP_TYPE = NRS_MSASPEC`` by
    USEAFTER, i.e. the same operability APT/MPT uses — is present in the
    local CRDS cache, fetching it from CRDS if it isn't. vMPT then loads it
    on startup, so the failed-/stuck-shutter map stays current.

    Non-fatal and offline-safe: if ``crds`` isn't installed, there's no
    network, or the check runs longer than ``timeout`` seconds, it quietly
    gives up and vMPT uses whatever is already cached. Runs at most once per
    process unless ``force=True``. Set ``VMPT_OPERABILITY_AUTOUPDATE=0`` to
    disable. Returns the resolved local path, or ``None``.
    """
    global _OPER_CHECK_DONE
    if not force:
        if _OPER_CHECK_DONE:
            return None
        _OPER_CHECK_DONE = True
        flag = os.environ.get("VMPT_OPERABILITY_AUTOUPDATE", "1").strip().lower()
        if flag in ("0", "false", "no", "off"):
            return None

    out: dict = {}

    def _fetch() -> None:
        try:
            os.environ.setdefault("CRDS_SERVER_URL", "https://jwst-crds.stsci.edu")
            os.environ.setdefault(
                "CRDS_PATH", os.path.expanduser("~/crds_cache"))
            os.environ.setdefault("CRDS_VERBOSITY", "-1")
            import logging
            logging.getLogger("CRDS").setLevel(logging.ERROR)
            import crds
            header = {
                "META.INSTRUMENT.NAME": "NIRSPEC",
                "META.EXPOSURE.TYPE": "NRS_MSASPEC",   # MOS branch of the rmap
                "META.OBSERVATION.DATE": _dt.date.today().isoformat(),
                "META.OBSERVATION.TIME": "00:00:00",
            }
            refs = crds.getreferences(
                header, reftypes=["msaoper"], observatory="jwst",
                context="jwst-operational")
            out["path"] = refs.get("msaoper")
        except BaseException as e:  # noqa: BLE001 — must never break startup
            out["error"] = e

    if timeout is None:
        _fetch()
    else:
        import threading
        t = threading.Thread(target=_fetch, daemon=True)
        t.start()
        t.join(timeout)
        if t.is_alive():
            print(f"[msa] operability check is taking >{timeout:.0f}s; using "
                  f"the cached reference for now (it will be ready next start).")
            return None

    if out.get("path"):
        print("[msa] operability checked — current MOS reference is "
              f"{os.path.basename(out['path'])}.")
        return out["path"]
    err = out.get("error")
    print("[msa] could not refresh operability from CRDS "
          f"({type(err).__name__ if err else 'unavailable'}); "
          "using the cached reference.")
    return None


def load_operability(crds_path: str | None = None) -> tuple[np.ndarray, np.ndarray]:
    operable = np.ones((4, 171, 365), dtype=bool)
    reason = np.zeros((4, 171, 365), dtype=np.int8)

    path = _find_msaoper_json(crds_path)
    if path is None:
        print("[msa] No CRDS msaoper JSON found; treating all shutters as operable.")
        return operable, reason

    with open(path) as f:
        data = json.load(f)

    for entry in data.get("msaoper", []):
        q = int(entry["Q"])
        d = int(entry["x"])
        s = int(entry["y"])
        if not (1 <= q <= 4 and 1 <= s <= 171 and 1 <= d <= 365):
            continue
        state = str(entry.get("state", "")).lower()
        # "stuck closed"/"closed" -> failed_closed; "stuck open"/"open" -> failed_open
        if "open" in state:
            operable[q - 1, s - 1, d - 1] = False
            reason[q - 1, s - 1, d - 1] = 2
        elif "closed" in state:
            operable[q - 1, s - 1, d - 1] = False
            reason[q - 1, s - 1, d - 1] = 1
    print(f"[msa] Loaded operability from {path}")
    return operable, reason
