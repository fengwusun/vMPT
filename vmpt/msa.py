"""NIRSpec MSA shutter grid loader and operability parser."""

import json
import os
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
        files = sorted(glob(os.path.join(cand, "jwst_nirspec_msaoper_*.json")))
        if files:
            return files[-1]
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
