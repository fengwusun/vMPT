"""Internal session JSON save/load (round-trips picks inside our tool, not for APT)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Optional

from app.empt_io import OpenShutter


@dataclass
class Session:
    pointing_ra_deg: float
    pointing_dec_deg: float
    pa_v3_deg: float
    disperser: str
    filter_name: str
    slitlet_height: int
    open_shutters: list[OpenShutter]
    highlighted: list[tuple[int, int, int]] = field(default_factory=list)
    image_path: Optional[str] = None
    catalog_path: Optional[str] = None
    tool_version: str = "1.0"
    created: Optional[str] = None


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def export_session_json(session: Session, path: str) -> None:
    """Serialize a Session to a JSON file. Schema documented in CONTEXT.md."""
    payload = {
        "version": session.tool_version,
        "created": session.created if session.created is not None else _utc_now_iso(),
        "pointing": {
            "ra_deg": float(session.pointing_ra_deg),
            "dec_deg": float(session.pointing_dec_deg),
            "apa_v3_deg": float(session.pa_v3_deg),
        },
        "instrument": {
            "disperser": session.disperser,
            "filter": session.filter_name,
            "slitlet_height": int(session.slitlet_height),
        },
        "open_shutters": [
            {
                "q": int(sh.q),
                "d": int(sh.d),
                "s": int(sh.s),
                "target_id": (str(sh.target_id) if sh.target_id is not None else None),
                "role": sh.role,
            }
            for sh in session.open_shutters
        ],
        "highlighted": [[int(q), int(s), int(d)] for (q, s, d) in session.highlighted],
        "image_path": session.image_path,
        "catalog_path": session.catalog_path,
    }
    with open(path, "w") as f:
        f.write(json.dumps(payload, indent=2))


def import_session_json(path: str) -> Session:
    """Parse a session JSON back into a Session dataclass. Raises ValueError on bad input."""
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise ValueError(f"could not read session JSON: {e}") from e

    if not isinstance(data, dict):
        raise ValueError("session JSON root must be an object")

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

    open_shutters: list[OpenShutter] = []
    for i, sh in enumerate(data["open_shutters"]):
        try:
            open_shutters.append(
                OpenShutter(
                    q=int(sh["q"]),
                    s=int(sh["s"]),
                    d=int(sh["d"]),
                    target_id=sh.get("target_id"),
                    role=sh.get("role", "target"),
                )
            )
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(f"malformed open_shutters[{i}]: {e}") from e

    highlighted: list[tuple[int, int, int]] = []
    for i, hl in enumerate(data.get("highlighted", [])):
        try:
            q, s, d = hl
            highlighted.append((int(q), int(s), int(d)))
        except (TypeError, ValueError) as e:
            raise ValueError(f"malformed highlighted[{i}]: {e}") from e

    return Session(
        pointing_ra_deg=ra,
        pointing_dec_deg=dec,
        pa_v3_deg=pa,
        disperser=disperser,
        filter_name=filter_name,
        slitlet_height=slitlet_height,
        open_shutters=open_shutters,
        highlighted=highlighted,
        image_path=data.get("image_path"),
        catalog_path=data.get("catalog_path"),
        tool_version=str(data.get("version", "1.0")),
        created=data.get("created"),
    )
