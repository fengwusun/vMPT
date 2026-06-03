"""Session JSON save/load.

`session.json` is written as a pure APT MPT plan JSON — no vMPT-only
keys, no file paths — so APT's MPT loader accepts it directly. A sibling
`vmpt_workspace.json` (same parent directory) carries the bits MPT
doesn't preserve: per-shutter target_id / role, highlighted set, image
+ catalog paths, slitlet height. vMPT reads both on import; APT only
sees the MPT file.

Old-style sessions (single file with a flat top-level `open_shutters`
list and a `pointing` block) are still accepted on import.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from vmpt.coords import V3_IDL_Y_ANGLE
from vmpt.empt_io import OpenShutter


SESSION_TOOL_VERSION = "1.4"
# Bundle filenames — chosen so the prefix telegraphs the role of each file:
#   MPT_*   → load these into APT MPT (plan JSON + primaries catalog)
#   vMPT_*  → vMPT-only state (image / catalog paths, target_id+role per shutter)
#   eMPT_*  → use these with the European eMPT pipeline (or any tool that
#             reads the eMPT shutter-mask / observed-targets / pointing-summary)
MPT_PLAN_FILENAME = "MPT_plan.json"
MPT_CATALOG_FILENAME = "MPT_catalog.cat"   # primaries catalog APT can import
WORKSPACE_FILENAME = "vMPT_workspace.json"
EMPT_OBSERVED_FILENAME = "eMPT_observed_targets.cat"
EMPT_POINTING_FILENAME = "eMPT_pointing_summary.txt"
EMPT_SHUTTER_MASK_FILENAME = "eMPT_shutter_mask.csv"
# Legacy / back-compat — older bundles used these names; the importer
# falls back to them if the current names are missing.
_LEGACY_MPT_PLAN_FILENAMES = ("session_MPT_plan.json",)
_LEGACY_WORKSPACE_FILENAMES = ("vmpt_workspace.json",)

# Maps slitlet height → APT's `msaSlitlet` enum value (per the reference
# G140H+G235H+G395H and a370 plans). 3 is the common case; vMPT only
# exposes 1, 3, 5 in the UI dropdown.
_MSA_SLITLET_ENUM = {
    1: "ONE_SHUTTER",
    2: "TWO_SHUTTER",
    3: "THREE_SHUTTER",
    5: "FIVE_SHUTTER",
}


@dataclass
class Session:
    pointing_ra_deg: float
    pointing_dec_deg: float
    pa_v3_deg: float                       # V3 PA, NOT the aperture PA
    disperser: str
    filter_name: str
    slitlet_height: int
    open_shutters: list[OpenShutter]
    highlighted: list[tuple[int, int, int]] = field(default_factory=list)
    image_path: Optional[str] = None
    wcs_sidecar_path: Optional[str] = None
    catalog_path: Optional[str] = None
    # vMPT 1.1+: multi-catalog support. Each entry: {"path": str,
    # "enabled": bool}. `catalog_path` (single) is preserved for
    # backward compatibility with vMPT 1.0 bundles — set to the first
    # entry's path when catalog_paths is non-empty.
    catalog_paths: list = field(default_factory=list)
    tool_version: str = SESSION_TOOL_VERSION
    created: Optional[str] = None
    name: Optional[str] = None


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _group_into_slitlets(
    open_shutters: list[OpenShutter],
) -> list[tuple[dict, Optional[str]]]:
    """Compress flat OpenShutter list back into MPT slitlets.

    Returns a list of (slitlet_dict, primary_id) tuples. Slitlets whose
    member shutters carry a non-None target_id are placed FIRST in the
    output (so positional `primaryIds[j] ↔ slitlets[j]` alignment works
    the way parse_mpt_json expects); manual / target-less slitlets come
    after, with primary_id=None.
    """
    # Bucket shutters by (q, d). For each bucket, sort by s and walk runs.
    buckets: dict[tuple[int, int], list[OpenShutter]] = {}
    for sh in open_shutters:
        buckets.setdefault((int(sh.q), int(sh.d)), []).append(sh)

    targeted: list[tuple[dict, str]] = []
    manual: list[dict] = []
    for (q, d), members in buckets.items():
        members.sort(key=lambda sh: int(sh.s))
        # Group consecutive s into runs.
        run: list[OpenShutter] = []
        runs: list[list[OpenShutter]] = []
        for sh in members:
            if run and int(sh.s) == int(run[-1].s) + 1:
                run.append(sh)
            else:
                if run:
                    runs.append(run)
                run = [sh]
        if run:
            runs.append(run)
        for r in runs:
            sl = {"q": q, "d": d, "s": int(r[0].s), "h": len(r)}
            # If any shutter in the run has a target_id, attach the most
            # common one as the slitlet's primary.
            tids = [sh.target_id for sh in r if sh.target_id is not None]
            if tids:
                primary = Counter(tids).most_common(1)[0][0]
                targeted.append((sl, str(primary)))
            else:
                manual.append(sl)
    # Targeted first → primaryIds positional alignment holds for them.
    out: list[tuple[dict, Optional[str]]] = [(sl, t) for sl, t in targeted]
    out.extend((sl, None) for sl in manual)
    return out


def _unfold_slitlets(slitlets: list[dict]) -> list[OpenShutter]:
    """Expand each {q,d,s,h} into h OpenShutter rows. Middle one →
    'target', others → 'sky'. Matches APT MPT semantics."""
    out: list[OpenShutter] = []
    for sl in slitlets:
        q = int(sl["q"]); d = int(sl["d"])
        s0 = int(sl["s"]); h = int(sl.get("h", 1))
        mid = h // 2
        for off in range(h):
            role = "target" if off == mid else "sky"
            out.append(OpenShutter(q=q, s=s0 + off, d=d, target_id=None, role=role))
    return out


def _spectral_offset_map(disperser: str) -> str:
    """APT's spectralOverlapShutterOffsetMap value, per disperser.

    APT names the map per individual grating (e.g. JWST_NIRSPEC_G395H),
    not per resolution class — verified against the reference G395H,
    G140H, and PRISM plan exports.
    """
    d = (disperser or "").upper()
    if d == "PRISM":
        return "JWST_NIRSPEC_PRISM"
    if d in {"G140M", "G235M", "G395M", "G140H", "G235H", "G395H"}:
        return f"JWST_NIRSPEC_{d}"
    return "JWST_NIRSPEC_PRISM"


def _build_mpt_payload(session: Session) -> dict:
    """Pure APT MPT plan JSON — no vMPT-only keys. Structure mirrors the
    reference G140H+G235H+G395H and a370 plan exports byte-for-byte
    schema-wise (we don't reproduce the full plannerSpecification, but
    the top-level shape is identical)."""
    apa = (float(session.pa_v3_deg) + V3_IDL_Y_ANGLE) % 360.0
    pairs = _group_into_slitlets(session.open_shutters)
    slitlets = [sl for sl, _ in pairs]
    primary_ids: list[int] = []
    for _, tid in pairs:
        if tid is None:
            continue
        try:
            primary_ids.append(int(tid))
        except (TypeError, ValueError):
            # Non-numeric target_ids can't sit in primaryIds (APT uses int);
            # silently skip — the workspace sidecar carries them losslessly.
            pass

    grating_filter = f"{session.disperser}_{session.filter_name}"
    msa_slitlet = _MSA_SLITLET_ENUM.get(int(session.slitlet_height), "THREE_SHUTTER")

    created = session.created if session.created is not None else _utc_now_iso()
    name = session.name or f"vMPT session — {created}"

    # `catalog_basename` is the identifier APT will look for in its Target
    # List database. We always export a primaries catalog file in the
    # bundle whose name MATCHES this basename (so importing the .cat under
    # its default name in APT lines up automatically with the plan).
    # Default is the stem of MPT_CATALOG_FILENAME (e.g. "MPT_catalog");
    # if a user-loaded catalog file is available we use its stem instead.
    from pathlib import Path as _P
    catalog_basename = _P(MPT_CATALOG_FILENAME).stem
    if session.catalog_path:
        catalog_basename = _P(session.catalog_path).stem

    n_targets = len(primary_ids)
    return {
        "instrument": "JWST/NIRSpec",
        "name": name,
        "aperturePA": apa,
        "theta": 0.0,
        "catalog": {
            "name": catalog_basename,
            "primariesName": catalog_basename,
            "fillersName": None,
            "primaries": None,
            "fillers": None,
        },
        "referencePointing": {
            "ra": float(session.pointing_ra_deg),
            "dec": float(session.pointing_dec_deg),
        },
        "configs": [{
            "name": "c1",
            "version": f"vMPT-{SESSION_TOOL_VERSION}",
            "info": {"fixedSlit": None},
            "masterBackground": False,
            "slitlets": slitlets,
            "exposures": [{
                "name": "c1e1",
                "gratingFilter": grating_filter,
                "msaSlitlet": msa_slitlet,
                "ra": float(session.pointing_ra_deg),
                "dec": float(session.pointing_dec_deg),
                "sourceIds": primary_ids,
            }],
            "primaryIds": primary_ids,
            "fillerIds": [],
        }],
        "stats": [{
            "name": "s0",
            "score": float(n_targets),
            "numberOfConfigurations": 1,
            "numberOfTargets": n_targets,
            "duration": 0.0,
            "totalDuration": 0.0,
        }],
        "errors": [],
        "plannerSpecification": {
            "gratingSpecification": {
                "gratings": [grating_filter],
                "allowContamination": False,
                "multiplexLimit": None,
                "multiplexingMinimum": None,
            },
            "planName": name,
            "planAngle": apa,
            "theta": None,
            "candidates": {
                "fillers": None,
                "primaries": catalog_basename,
                "catalog": catalog_basename,
                "slitSources": None,
            },
            "slitSpecification": {
                "sweetSpot": "PERCENT_0",
                "slitlet": msa_slitlet,
            },
            "searchParameters": {
                "useWeights": False,
                "enableMonteCarlo": False,
                "monteCarloShuffles": None,
                "ignoreStuckOpen": False,
                "spectralOverlapThreshold": 1.5,
                "numberOfConfigurations": 1,
                "allowMultiSourceShutters": False,
                "spectralOverlapShutterOffsetMap": _spectral_offset_map(session.disperser),
            },
            "slitSearchSpecification": None,
            "maskingSpecification": {
                "fillerMask": None,
                "primaryMask": None,
                "noGapFiller": False,
                "noGapPrimary": False,
                "noRedCutoffPrimary": False,
                "noBlueCutoffPrimary": False,
                "noRedCutoffFiller": False,
                "noBlueCutoffFiller": False,
            },
            "pointingSpecification": {
                "ditherType": "NONE",
                "shouldNod": False,
                "pointingMode": "GRID_SEARCH",
                "numberOfNods": 5,
                "fixedDitherOffsets": [{"spatial": 0, "dispersion": 0}],
                "fixedPointings": None,
                "partiallyCompletedPrimaries": False,
                "minPrimaryDitherPoints": 0,
                "partiallyCompletedFillers": False,
                "minFillerDitherPoints": None,
            },
            "searchGridSpecification": {
                "searchArea": {
                    "width": 10.0,
                    "height": 10.0,
                    "center": {"x": 0.0, "y": 0.0},
                    "ylength": 10.0,
                    "xlength": 10.0,
                    "corners": [],
                    "offsetToBottomCorner": {"x": -5.0, "y": -5.0},
                },
                "searchAreaCenter": {
                    "ra": float(session.pointing_ra_deg),
                    "dec": float(session.pointing_dec_deg),
                },
                "searchAreaHeight": 10.0,
                "searchAreaWidth": 10.0,
                "searchStepSize": 0.5,
            },
            "wavelengthRangeSpecification": {
                "primaryRange1": False,
                "primaryRange2": False,
                "primaryRange3": False,
                "primaryRange4": False,
                "primaryRange5": False,
                "fillerRange1": False,
                "fillerRange2": False,
                "fillerRange3": False,
                "fillerRange4": False,
                "fillerRange5": False,
            },
        },
    }


def _build_workspace_payload(session: Session) -> dict:
    """vMPT-only extras: per-shutter target_id/role, highlighted set,
    image + catalog paths. Written next to the MPT session.json so APT
    never sees it."""
    return {
        "vmpt_version": session.tool_version,
        "created": session.created if session.created is not None else _utc_now_iso(),
        "pa_v3_deg": float(session.pa_v3_deg),
        "slitlet_height": int(session.slitlet_height),
        "open_shutters": [
            {
                "q": int(sh.q), "d": int(sh.d), "s": int(sh.s),
                "target_id": (str(sh.target_id) if sh.target_id is not None else None),
                "role": sh.role,
            }
            for sh in session.open_shutters
        ],
        "highlighted": [[int(q), int(s), int(d)] for (q, s, d) in session.highlighted],
        "image_path": session.image_path,
        "wcs_sidecar_path": session.wcs_sidecar_path,
        "catalog_path": session.catalog_path,
        "catalog_paths": [
            {"path": str(e.get("path")), "enabled": bool(e.get("enabled", True))}
            for e in (session.catalog_paths or [])
            if e.get("path")
        ],
    }


def export_session_json(session: Session, path: str) -> None:
    """Write the session as an MPT-format plan JSON at `path`, AND a
    sibling `vmpt_workspace.json` carrying the vMPT-only extras."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        f.write(json.dumps(_build_mpt_payload(session), indent=2))
    workspace_path = p.parent / WORKSPACE_FILENAME
    with open(workspace_path, "w") as f:
        f.write(json.dumps(_build_workspace_payload(session), indent=2))


def _parse_grating_filter(gf: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    if not isinstance(gf, str):
        return None, None
    for sep in ("_", "/"):
        if sep in gf:
            a, b = gf.split(sep, 1)
            return a, b
    return None, None


def _parse_catalog_paths(
    raw: object,
    legacy_single: Optional[str] = None,
) -> list:
    """Normalise the workspace's catalog_paths field into a clean list of
    `{"path": str, "enabled": bool}` dicts.

    Older bundles (vMPT ≤ 1.0) only stored a single `catalog_path`; if
    the new field is absent or empty, synthesise a single-entry list
    from that fallback so multi-catalog UIs still have one row to show.
    """
    out: list = []
    if isinstance(raw, list):
        for entry in raw:
            if isinstance(entry, dict) and entry.get("path"):
                out.append({
                    "path": str(entry["path"]),
                    "enabled": bool(entry.get("enabled", True)),
                })
            elif isinstance(entry, str) and entry:
                out.append({"path": entry, "enabled": True})
    if not out and legacy_single:
        out.append({"path": str(legacy_single), "enabled": True})
    return out


def _import_mpt(data: dict, sidecar: dict) -> Session:
    """Combine an MPT-format session.json + vmpt_workspace.json into a Session."""
    configs = data.get("configs") or []
    if not configs or not isinstance(configs[0], dict):
        raise ValueError("session has no usable configs")
    cfg = configs[0]
    # Pick first DISPERSED exposure (same logic as mpt_io.parse_mpt_json).
    exps = cfg.get("exposures") or []
    primary_exp = None
    for e in exps:
        if isinstance(e, dict) and e.get("gratingFilter"):
            primary_exp = e
            break
    if primary_exp is None and exps and isinstance(exps[0], dict):
        primary_exp = exps[0]
    ra = dec = None
    grating = filt = None
    if primary_exp is not None:
        try:
            ra = float(primary_exp["ra"]); dec = float(primary_exp["dec"])
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(f"malformed exposure ra/dec: {e}") from e
        grating, filt = _parse_grating_filter(primary_exp.get("gratingFilter"))

    # PA from the workspace's exact V3 PA, else derived from APA.
    try:
        if "pa_v3_deg" in sidecar:
            pa_v3 = float(sidecar["pa_v3_deg"])
        else:
            apa = float(data["aperturePA"])
            pa_v3 = (apa - V3_IDL_Y_ANGLE) % 360.0
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(f"malformed aperturePA / pa_v3_deg: {e}") from e

    slitlet_height = int(sidecar.get("slitlet_height", 3))

    # Prefer the workspace's lossless open_shutters list (carries target_id
    # and role). Fall back to unfolding the MPT slitlets, which loses
    # those — but at least restores the shutter grid.
    if isinstance(sidecar.get("open_shutters"), list):
        opens: list[OpenShutter] = []
        for i, sh in enumerate(sidecar["open_shutters"]):
            try:
                opens.append(OpenShutter(
                    q=int(sh["q"]), s=int(sh["s"]), d=int(sh["d"]),
                    target_id=sh.get("target_id"),
                    role=sh.get("role", "target"),
                ))
            except (KeyError, TypeError, ValueError) as e:
                raise ValueError(f"malformed workspace open_shutters[{i}]: {e}") from e
    else:
        opens = _unfold_slitlets(cfg.get("slitlets") or [])
        # If primaryIds positionally aligns with slitlets, recover target_id
        # for the middle ('target') shutter of each slitlet.
        primary_ids = cfg.get("primaryIds") or []
        slitlets = cfg.get("slitlets") or []
        if primary_ids and len(primary_ids) <= len(slitlets):
            cursor = 0
            for j, sl in enumerate(slitlets):
                h = int(sl.get("h", 1))
                tid = (str(primary_ids[j]) if j < len(primary_ids) else None)
                if tid is not None:
                    for off in range(h):
                        opens[cursor + off].target_id = tid
                cursor += h

    highlighted: list[tuple[int, int, int]] = []
    for i, hl in enumerate(sidecar.get("highlighted") or []):
        try:
            q, s, d = hl
            highlighted.append((int(q), int(s), int(d)))
        except (TypeError, ValueError) as e:
            raise ValueError(f"malformed workspace highlighted[{i}]: {e}") from e

    if grating is None or filt is None:
        raise ValueError("could not determine disperser/filter from session")

    catalog_paths = _parse_catalog_paths(
        sidecar.get("catalog_paths"), sidecar.get("catalog_path"),
    )

    return Session(
        pointing_ra_deg=ra if ra is not None else 0.0,
        pointing_dec_deg=dec if dec is not None else 0.0,
        pa_v3_deg=pa_v3,
        disperser=str(grating),
        filter_name=str(filt),
        slitlet_height=slitlet_height,
        open_shutters=opens,
        highlighted=highlighted,
        image_path=sidecar.get("image_path"),
        wcs_sidecar_path=sidecar.get("wcs_sidecar_path"),
        catalog_path=sidecar.get("catalog_path"),
        catalog_paths=catalog_paths,
        tool_version=str(sidecar.get("vmpt_version", SESSION_TOOL_VERSION)),
        created=sidecar.get("created") or data.get("created"),
        name=data.get("name"),
    )


def _import_legacy(data: dict) -> Session:
    """Import the original vMPT-only session schema (flat `open_shutters`)."""
    for key in ("pointing", "instrument", "open_shutters"):
        if key not in data:
            raise ValueError(f"missing required key: {key!r}")

    pointing = data["pointing"]
    instrument = data["instrument"]
    try:
        ra = float(pointing["ra_deg"])
        dec = float(pointing["dec_deg"])
        pa = float(pointing["apa_v3_deg"])
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(f"malformed pointing: {e}") from e

    try:
        disperser = str(instrument["disperser"])
        filter_name = str(instrument["filter"])
    except (KeyError, TypeError) as e:
        raise ValueError(f"malformed instrument: {e}") from e
    slitlet_height = int(instrument.get("slitlet_height", 3))

    opens: list[OpenShutter] = []
    for i, sh in enumerate(data["open_shutters"]):
        try:
            opens.append(OpenShutter(
                q=int(sh["q"]), s=int(sh["s"]), d=int(sh["d"]),
                target_id=sh.get("target_id"),
                role=sh.get("role", "target"),
            ))
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(f"malformed open_shutters[{i}]: {e}") from e

    highlighted: list[tuple[int, int, int]] = []
    for i, hl in enumerate(data.get("highlighted", [])):
        try:
            q, s, d = hl
            highlighted.append((int(q), int(s), int(d)))
        except (TypeError, ValueError) as e:
            raise ValueError(f"malformed highlighted[{i}]: {e}") from e

    catalog_paths = _parse_catalog_paths(
        data.get("catalog_paths"), data.get("catalog_path"),
    )

    return Session(
        pointing_ra_deg=ra,
        pointing_dec_deg=dec,
        pa_v3_deg=pa,
        disperser=disperser,
        filter_name=filter_name,
        slitlet_height=slitlet_height,
        open_shutters=opens,
        highlighted=highlighted,
        image_path=data.get("image_path"),
        wcs_sidecar_path=data.get("wcs_sidecar_path"),
        catalog_path=data.get("catalog_path"),
        catalog_paths=catalog_paths,
        tool_version=str(data.get("version", "1.0")),
        created=data.get("created"),
    )


def _load_json_or_empty(p: Path) -> dict:
    """Return {} if the file is missing or malformed."""
    if not p.exists():
        return {}
    try:
        with open(p) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def import_session_json(path: str) -> Session:
    """Parse a session JSON back into a Session. The user can point at
    EITHER file in a bundle:

      • `session_MPT_plan.json` → pure MPT plan; we look for a sibling
        `vmpt_workspace.json` to merge in target_ids, roles, image path.
      • `vmpt_workspace.json` → vMPT extras; we look for a sibling
        `session_MPT_plan.json` (or any `*plan*.json` matching MPT shape)
        to pull pointing / PA / disperser / slitlet geometry.

    Legacy single-file sessions (`open_shutters` at top level) still load.
    """
    p = Path(path)
    try:
        with open(p) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise ValueError(f"could not read session JSON: {e}") from e

    if not isinstance(data, dict):
        raise ValueError("session JSON root must be an object")

    # Case A: user pointed at the workspace sidecar — find the MPT plan sibling.
    if "vmpt_version" in data or "open_shutters" in data and "pointing" not in data:
        sidecar = data
        mpt_data: dict = {}
        for fname in (MPT_PLAN_FILENAME, *_LEGACY_MPT_PLAN_FILENAMES):
            mpt_data = _load_json_or_empty(p.parent / fname)
            if mpt_data:
                break
        if not mpt_data:
            # Last resort: any *.json sibling whose shape says MPT plan.
            for candidate in sorted(p.parent.glob("*.json")):
                if candidate == p:
                    continue
                d = _load_json_or_empty(candidate)
                if "configs" in d and "aperturePA" in d:
                    mpt_data = d
                    break
        if not mpt_data:
            raise ValueError(
                f"workspace at {p.name} needs a sibling MPT plan JSON "
                f"(expected {MPT_PLAN_FILENAME}); none found"
            )
        return _import_mpt(mpt_data, sidecar)

    # Case B: user pointed at the MPT plan — find the workspace sidecar.
    if "configs" in data:
        sidecar: dict = {}
        for fname in (WORKSPACE_FILENAME, *_LEGACY_WORKSPACE_FILENAMES):
            sidecar = _load_json_or_empty(p.parent / fname)
            if sidecar:
                break
        return _import_mpt(data, sidecar)

    # Case C: legacy single-file vMPT session.
    if "open_shutters" in data:
        return _import_legacy(data)

    raise ValueError(
        "not a recognized session JSON (no 'configs', 'open_shutters', "
        "or 'vmpt_version' top-level key)"
    )
