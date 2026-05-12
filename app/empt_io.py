"""eMPT-compatible exporters (observed_targets.cat, pointing_summary.txt, shutter_mask.csv).

The shutter_mask.csv tiling reproduces what ``make_csv_file`` in
``refs/eMPT_v1/reference_files/shutter_routines_new.f90`` writes (lines 535-678).

Fortran code summary (1-indexed throughout):

    msamap(kk, ii, jj)     ! kk=quadrant 1..4, ii=dispersion 1..365, jj=shutter row 1..171

    do ir = 1, 365                 ! csv data row, top half
        kk=1, ii=ir, jj=1..171  -> chars  1..341 step 2  (cells  1..171)  Q1
        kk=2, ii=ir, jj=1..171  -> chars 343..683 step 2 (cells 172..342) Q2
    do ir = 366, 730               ! csv data row, bottom half
        kk=3, ii=ir-365, jj=1..171 -> cells   1..171   Q3
        kk=4, ii=ir-365, jj=1..171 -> cells 172..342   Q4

So with our ``(q, s, d)`` convention (q in 1..4, s in 1..171, d in 1..365):

    csv_row (1..730), csv_col (1..342):
        top half    (csv_row 1..365):   d = csv_row;        q = 1 if csv_col<=171 else 2
        bottom half (csv_row 366..730): d = csv_row - 365;  q = 3 if csv_col<=171 else 4
        s = csv_col                    if csv_col <= 171
        s = csv_col - 171              otherwise

i.e. each CSV data row holds a single dispersion column ``d``; within that row the
171 cells of Q{1,3} sit side-by-side with the 171 cells of Q{2,4}, indexed by
shutter-row ``s``. There is **no** transpose or reverse.

Cell alphabet (precedence top-down):

    'x' failed-closed   (operability)
    's' failed-open     (operability)
    '0' commanded open  (user's pick)
    '1' commanded closed / functional

Each line (header + 730 data rows) is exactly 683 characters wide; the header text
is padded with spaces to that width to match the reference byte-for-byte.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------


@dataclass
class Pointing:
    """Single pointing for the export bundle."""

    ra_deg: float
    dec_deg: float
    apa_v3_deg: float
    pa_ap_deg: float | None = None


@dataclass
class OpenShutter:
    """A user-commanded open shutter."""

    q: int
    s: int
    d: int
    target_id: int | str | None = None
    role: str = "target"


# ---------------------------------------------------------------------------
# observed_targets.cat
# ---------------------------------------------------------------------------


_OBSCAT_HEADER = "# No   No_sub      No_cat    Pr    RA[deg]     Dec[deg]"


def write_observed_targets_cat(path: str, targets: list[dict]) -> None:
    """Write eMPT-style observed_targets.cat.

    Each ``targets`` dict must contain at least ``No_cat``, ``Pr``, ``ra_deg``
    and ``dec_deg``. ``No_sub`` is optional (defaults to the running ``No``).
    """
    lines = [_OBSCAT_HEADER]
    for i, t in enumerate(targets, start=1):
        no = i
        no_sub = int(t.get("No_sub", no))
        no_cat = int(t["No_cat"])
        pr = int(t["Pr"])
        ra = float(t["ra_deg"])
        dec = float(t["dec_deg"])
        # Field widths chosen to match the reference layout:
        #   "   1        1        14170   1   53.1633910  -27.7756740"
        line = (
            f"{no:4d}{no_sub:9d}{no_cat:13d}{pr:4d}"
            f"   {ra:10.7f}  {dec:11.7f}"
        )
        lines.append(line)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# pointing_summary.txt
# ---------------------------------------------------------------------------


def write_pointing_summary_txt(
    path: str,
    pointing: Pointing,
    disperser: str,
    filter_name: str,
    n_targets_total: int = 0,
    n_targets_accepted: int = 0,
) -> None:
    """Write a pointing_summary.txt matching the reference layout."""
    pa_ap = pointing.pa_ap_deg if pointing.pa_ap_deg is not None else pointing.apa_v3_deg
    pa_v3 = pointing.apa_v3_deg
    mask = f"{disperser}/{filter_name}"
    lines = [
        "",
        "",
        " ------- Summary -----------",
        "",
        " Single Pointing",
        "",
        " Pointing information:",
        "",
        " RA, Dec of Central Pointing:",
        f" Nod 0:    {pointing.ra_deg:11.7f}  {pointing.dec_deg:12.7f}",
        "",
        " Official Assigned APT/MPT roll angle:",
        f" PA_AP:    {pa_ap:9.6f}",
        "",
        " Actual MSA roll angle:",
        f" PA_AP:    {pa_ap:9.6f}",
        f" PA_V3:   {pa_v3:10.6f}",
        "",
        " RA, Dec of Nodded Pointings:",
        f" Nod 1:    {pointing.ra_deg:11.7f}  {pointing.dec_deg:12.7f}",
        f" Nod 2:    {pointing.ra_deg:11.7f}  {pointing.dec_deg:12.7f}",
        "",
        f" MSA mask intended for: {mask:<20s}",
        "",
        f" Total number of targets in input catalog:    {n_targets_total:12d}",
        f" Number of accepted targets:                  {n_targets_accepted:12d}",
        "",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def parse_pointing_summary_txt(path: str) -> dict:
    """Tiny round-trip helper used by the test suite."""
    out: dict = {}
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("Nod 0:"):
                parts = line.split()
                # "Nod 0: RA Dec"
                out["ra_deg"] = float(parts[2])
                out["dec_deg"] = float(parts[3])
            elif line.startswith("PA_AP:") and "pa_ap_deg" not in out:
                out["pa_ap_deg"] = float(line.split()[1])
            elif line.startswith("PA_V3:"):
                out["pa_v3_deg"] = float(line.split()[1])
    return out


# ---------------------------------------------------------------------------
# shutter_mask.csv
# ---------------------------------------------------------------------------


_CSV_HEADER_TEXT = (
    "# This CSV indicates which shutters should be open/closed on the MSA"
    " - created by ESA NIRSpec Team"
)
_CSV_LINE_WIDTH = 683  # 342 cells + 341 commas


def _csv_header_line() -> str:
    return _CSV_HEADER_TEXT.ljust(_CSV_LINE_WIDTH)


def write_shutter_mask_csv(
    path: str,
    open_shutters: list[OpenShutter],
    operable: np.ndarray,
    reason: np.ndarray,
) -> None:
    """Write the 730 x 342 shutter-mask grid.

    Parameters
    ----------
    path:
        Output path.
    open_shutters:
        Commanded-open shutters; written as ``'0'`` (overrides ``'1'`` but not
        operability failures).
    operable:
        ``(4, 171, 365)`` bool array. ``True`` = functional.
    reason:
        ``(4, 171, 365)`` int8 array. ``1`` = failed-closed (-> 'x'),
        ``2`` = failed-open (-> 's'); other values fall back to operable/open.
    """
    if operable.shape != (4, 171, 365):
        raise ValueError(f"operable must be (4,171,365), got {operable.shape}")
    if reason.shape != (4, 171, 365):
        raise ValueError(f"reason must be (4,171,365), got {reason.shape}")

    # Start every cell at '1' (closed, functional), then layer in overrides.
    # cells[q-1, s-1, d-1] is a single character.
    cells = np.full((4, 171, 365), "1", dtype="<U1")
    cells[(reason == 2) | (~operable & (reason == 0))] = "1"  # safety default
    cells[reason == 2] = "s"
    cells[reason == 1] = "x"
    for sh in open_shutters:
        if reason[sh.q - 1, sh.s - 1, sh.d - 1] in (1, 2):
            # Don't override a failed shutter — eMPT precedence puts operability first.
            continue
        cells[sh.q - 1, sh.s - 1, sh.d - 1] = "0"

    out_lines = [_csv_header_line()]
    # Top half: d = csv_row (1..365); cells[:,:171] from q=1, cells[:,171:] from q=2
    for d in range(1, 366):
        left = cells[0, :, d - 1]   # Q1, indexed by s (171 cells)
        right = cells[1, :, d - 1]  # Q2
        out_lines.append(",".join(np.concatenate([left, right]).tolist()))
    # Bottom half: d = csv_row - 365; left=Q3, right=Q4
    for d in range(1, 366):
        left = cells[2, :, d - 1]
        right = cells[3, :, d - 1]
        out_lines.append(",".join(np.concatenate([left, right]).tolist()))

    with open(path, "w") as f:
        f.write("\n".join(out_lines) + "\n")


def parse_shutter_mask_csv(path: str) -> tuple[np.ndarray, np.ndarray, list[OpenShutter]]:
    """Inverse of :func:`write_shutter_mask_csv` — used by the test suite.

    Returns ``(operable, reason, open_shutters)`` with the same shapes/conventions
    as the writer expects.
    """
    with open(path) as f:
        lines = f.read().splitlines()
    if len(lines) != 731:
        raise ValueError(f"expected 731 lines, got {len(lines)}")
    data = lines[1:]
    operable = np.ones((4, 171, 365), dtype=bool)
    reason = np.zeros((4, 171, 365), dtype=np.int8)
    open_shutters: list[OpenShutter] = []
    for ri, line in enumerate(data):
        cells = line.split(",")
        if len(cells) != 342:
            raise ValueError(f"row {ri} has {len(cells)} cells, expected 342")
        if ri < 365:
            d = ri + 1
            q_left, q_right = 1, 2
        else:
            d = ri - 365 + 1
            q_left, q_right = 3, 4
        for s in range(1, 172):
            for q, c in ((q_left, cells[s - 1]), (q_right, cells[171 + s - 1])):
                if c == "x":
                    operable[q - 1, s - 1, d - 1] = False
                    reason[q - 1, s - 1, d - 1] = 1
                elif c == "s":
                    operable[q - 1, s - 1, d - 1] = False
                    reason[q - 1, s - 1, d - 1] = 2
                elif c == "0":
                    open_shutters.append(OpenShutter(q=q, s=s, d=d))
                # '1' is the default; nothing to do
    return operable, reason, open_shutters
