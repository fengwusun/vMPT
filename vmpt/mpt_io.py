"""Import MSA plans from APT/MPT exports (JSON plan files, shutter CSVs, .aptx archives)."""

from __future__ import annotations

import io
import json
import re
import tempfile
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from vmpt.coords import V3_IDL_Y_ANGLE
from vmpt.empt_io import OpenShutter

# STScI public APT downloader URL pattern. e.g. .../apt/1208/ → 1208.aptx.
APT_URL_TEMPLATE = "https://www.stsci.edu/jwst-program-info/download/jwst/apt/{program_id}/"


@dataclass
class MPTSlitlet:
    """A single slitlet entry from MPT JSON: starting shutter (q, s, d)
    and slitlet height h. Unfolds to h shutters at rows s, s+1, …, s+h-1."""
    q: int
    s: int
    d: int
    h: int = 1
    primary_id: Optional[int] = None  # primary target source ID, if known


@dataclass
class MPTPlan:
    """One MSA plan / config from an APT MPT JSON file."""
    name: str
    aperture_pa_deg: float          # APT's "aperturePA" — NIRSpec APA
    v3_pa_deg: float                # derived: APA - V3IdlYAngle (mod 360)
    ra_deg: Optional[float] = None  # from primary exposure if present
    dec_deg: Optional[float] = None
    grating: Optional[str] = None
    filter_name: Optional[str] = None
    slitlets: list[MPTSlitlet] = None
    catalog_name: Optional[str] = None
    primary_ids: list[int] = None
    n_open_shutters: int = 0

    def to_open_shutters(self) -> list[OpenShutter]:
        """Unfold slitlets into a flat list of OpenShutter entries.

        Each slitlet at (q, s, d, h) maps to h shutters at rows
        s, s+1, …, s+h-1 in column d of quadrant q. The middle shutter
        (s + h//2) is the "target" row; the others are "sky".
        """
        out: list[OpenShutter] = []
        for sl in self.slitlets or []:
            mid_offset = sl.h // 2
            for off in range(sl.h):
                role = "target" if off == mid_offset else "sky"
                tid = str(sl.primary_id) if sl.primary_id is not None else None
                out.append(OpenShutter(
                    q=int(sl.q), s=int(sl.s + off), d=int(sl.d),
                    target_id=tid, role=role,
                ))
        return out


def parse_mpt_json(path: str) -> list[MPTPlan]:
    """Parse an APT/MPT JSON plan file into one MPTPlan per `configs` entry.

    APT's MSA-Planner exports JSON with this shape::

        {
          "aperturePA": <degrees>,         # NIRSpec aperture PA on sky
          "theta": <degrees>,              # MSA roll (we don't use it directly)
          "catalog": {"name": ..., "primariesName": ..., ...},
          "configs": [
              {"name": "c1 : foo plan 1",
               "slitlets": [{"q":1,"d":12,"s":50,"h":3}, ...],
               "exposures": [{"ra": ..., "dec": ..., "gratingFilter": "PRISM_CLEAR",
                              "sourceIds": [...], ...}],
               "primaryIds": [...], "fillerIds": [...]
              },
              ...
          ]
        }

    Returns one MPTPlan per config. The plan's RA/Dec/grating/filter come
    from the first exposure of the config. The V3 PA is derived from
    aperturePA - V3IdlYAngle so vMPT can drive the overlay directly.

    Raises ValueError on malformed input.
    """
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise ValueError(f"Could not read MPT JSON {path}: {e}") from e

    if not isinstance(data, dict) or "configs" not in data:
        raise ValueError(f"{path}: not an MPT plan JSON (missing 'configs')")

    try:
        apa = float(data["aperturePA"])
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(f"{path}: missing or invalid 'aperturePA'") from e
    v3pa = (apa - V3_IDL_Y_ANGLE) % 360.0

    cat = data.get("catalog", {}) or {}
    catalog_name = cat.get("name") or cat.get("primariesName")

    plans: list[MPTPlan] = []
    for i, cfg in enumerate(data["configs"]):
        if not isinstance(cfg, dict):
            continue
        name = str(cfg.get("name", f"config_{i}"))
        # Parse slitlets
        raw_sl = cfg.get("slitlets") or []
        slitlets: list[MPTSlitlet] = []
        primary_ids = cfg.get("primaryIds") or []
        # Map each slitlet to a primary target id if there's a positional
        # correspondence (APT typically lists slitlets in primaryIds order).
        for j, sl in enumerate(raw_sl):
            try:
                slitlets.append(MPTSlitlet(
                    q=int(sl["q"]),
                    s=int(sl["s"]),
                    d=int(sl["d"]),
                    h=int(sl.get("h", 1)),
                    primary_id=(int(primary_ids[j]) if j < len(primary_ids) else None),
                ))
            except (KeyError, TypeError, ValueError) as e:
                raise ValueError(
                    f"{path}: malformed slitlet [{i}][{j}]: {e}"
                ) from e

        # Derive RA/Dec and grating from the first DISPERSED exposure.
        # APT plans interleave target-acquisition/imaging steps (which have
        # `gratingFilter: null`) with the science exposure that carries the
        # grating + filter — sometimes the dispersed exposure isn't first.
        # Fall back to exposures[0] for plans with no dispersed step at all
        # (e.g. an MSA shutter-mask preview).
        exps = cfg.get("exposures") or []
        ra = dec = None
        grating = filt = None
        primary = None
        for e in exps:
            if isinstance(e, dict) and e.get("gratingFilter"):
                primary = e
                break
        if primary is None and exps:
            primary = exps[0] if isinstance(exps[0], dict) else None
        if primary is not None:
            try:
                ra = float(primary.get("ra")) if primary.get("ra") is not None else None
                dec = float(primary.get("dec")) if primary.get("dec") is not None else None
            except (TypeError, ValueError):
                pass
            gf = primary.get("gratingFilter") or ""
            if "_" in gf:
                grating, filt = gf.split("_", 1)
            elif "/" in gf:
                grating, filt = gf.split("/", 1)

        # Count physical open shutters (each slitlet contributes h)
        n_open = sum(sl.h for sl in slitlets)
        plans.append(MPTPlan(
            name=name, aperture_pa_deg=apa, v3_pa_deg=v3pa,
            ra_deg=ra, dec_deg=dec,
            grating=grating, filter_name=filt,
            slitlets=slitlets, catalog_name=catalog_name,
            primary_ids=list(primary_ids), n_open_shutters=n_open,
        ))

    return plans


def list_mpt_plans_in_aptx(aptx_path: str) -> list[str]:
    """Return the names of MPT-style JSON files embedded in an .aptx archive.

    .aptx files are zip archives containing one or more <plan>.json files
    in MPT format, alongside the proposal XML, manifest, and pointing
    files. We filter to JSON files that look like MPT plans (top-level
    keys include 'configs' and 'aperturePA').
    """
    out: list[str] = []
    with zipfile.ZipFile(aptx_path) as zf:
        for info in zf.infolist():
            if not info.filename.lower().endswith(".json"):
                continue
            try:
                with zf.open(info) as fh:
                    data = json.load(fh)
            except (json.JSONDecodeError, OSError):
                continue
            if isinstance(data, dict) and "configs" in data and "aperturePA" in data:
                out.append(info.filename)
    return sorted(out)


def parse_mpt_json_in_aptx(aptx_path: str, member_name: str) -> list[MPTPlan]:
    """Extract one MPT JSON entry from an .aptx archive into a list of
    MPTPlan (delegating to parse_mpt_json after extraction)."""
    with zipfile.ZipFile(aptx_path) as zf, tempfile.NamedTemporaryFile(
        suffix=".json", delete=False, mode="wb",
    ) as out:
        out.write(zf.read(member_name))
        tmp_path = out.name
    try:
        return parse_mpt_json(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def download_apt_program(program_id: int | str, dest_path: Optional[str] = None) -> str:
    """Fetch <program_id>.aptx from STScI's public proposal-info URL.

    Returns the local filesystem path of the downloaded archive. If
    `dest_path` is None, writes to a temp file. Raises ValueError if the
    fetch fails or the response isn't a recognizable .aptx.
    """
    pid = str(program_id).strip()
    if not re.fullmatch(r"\d+", pid):
        raise ValueError(f"program_id must be an integer; got {program_id!r}")
    url = APT_URL_TEMPLATE.format(program_id=pid)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            payload = resp.read()
    except Exception as e:  # noqa: BLE001
        raise ValueError(f"Could not fetch APT {pid} from {url}: {e}") from e
    # Sanity check — zip magic bytes
    if not payload.startswith(b"PK"):
        raise ValueError(f"APT {pid}: server didn't return a zip file (got "
                         f"{len(payload)} bytes starting {payload[:8]!r})")
    if dest_path is None:
        tmp = tempfile.NamedTemporaryFile(prefix=f"apt_{pid}_", suffix=".aptx", delete=False)
        tmp.write(payload)
        tmp.close()
        return tmp.name
    Path(dest_path).write_bytes(payload)
    return str(dest_path)


def parse_shutter_csv(path: str) -> list[OpenShutter]:
    """Parse a shutter_mask.csv exported by APT/MPT (or eMPT).

    Format: 731 lines (1 header + 730 data rows × 342 cells). Cells
    encode shutter state: '0' = commanded open (user's pick), '1' =
    closed, 'x' = failed-closed (operability), 's' = failed-open.
    We return only the '0' cells.

    Tiling matches our writer: CSV row 1..365 = Q1 + Q2 by d, CSV row
    366..730 = Q3 + Q4 by d (d = row - 365); CSV col 1..171 = s for
    Q1/Q3, CSV col 172..342 = s for Q2/Q4 (s = col or col - 171).
    """
    try:
        with open(path) as f:
            lines = f.read().splitlines()
    except OSError as e:
        raise ValueError(f"Could not read shutter CSV {path}: {e}") from e

    # Strip header (first line that doesn't start with a digit/letter)
    data_lines = [ln for ln in lines if ln and not ln.lstrip().startswith("#")]
    if len(data_lines) != 730:
        raise ValueError(
            f"{path}: expected 730 data rows in shutter CSV, got {len(data_lines)}"
        )

    out: list[OpenShutter] = []
    for ir, line in enumerate(data_lines):
        cells = [c.strip() for c in line.split(",")]
        if len(cells) != 342:
            raise ValueError(
                f"{path}: row {ir+1} has {len(cells)} cells; expected 342"
            )
        # Determine quadrant pair and d-index from row index
        if ir < 365:
            d = ir + 1
            q_left, q_right = 1, 2
        else:
            d = ir - 365 + 1
            q_left, q_right = 3, 4
        for ic, c in enumerate(cells):
            if c != "0":
                continue
            if ic < 171:
                q = q_left
                s = ic + 1
            else:
                q = q_right
                s = ic - 171 + 1
            out.append(OpenShutter(q=q, s=s, d=d, target_id=None, role="manual"))
    return out
