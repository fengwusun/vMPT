# Quick start

A two-minute tour. Assumes you've already run `pip install jwst-vmpt`
(see [Installation](installation.md)).

## 1. Launch the app

```bash
vmpt
# browser opens at http://localhost:5099/vmpt by default
# pass --port to use a different port
```

If a tab doesn't open automatically, visit the URL yourself.

## 2. Load an example

In the **Input** tab (left sidebar), click one of the two example
buttons:

- **Load Abell 370 example** — JWST/NIRCam F182M + F200W + F210M
  three-band FITS (~42 MB). Includes a target catalog and an APT
  MPT plan from GTO-1208.
- **Load RXCJ0600 example** — JPG + WCS-sidecar pair (~17 MB).
  Demonstrates the JPG + WCS workflow. Includes a 28k-source
  target catalog.

If you haven't downloaded the examples yet, run

```bash
vmpt examples download
```

(from any working directory — since v1.3.1 the tarball lands in
`~/.vmpt/examples/`, which the buttons search automatically). See
[Troubleshooting](troubleshooting.md) if the buttons stay dead
after the download finishes.

## 3. Aim the MSA

Switch to the **Pointing** tab. The pointing centre auto-fills to
the image centre. Drag the V3 PA slider (or type into the APA box)
to rotate the MSA quadrants on the image.

Pick a Disperser / Filter combo — e.g. PRISM/CLEAR for low-res,
G395H/F290LP for high-res near-IR.

## 4. Pick shutters

Switch to the **Settings** tab → choose `N-shutter slitlet` (default
3). Then on the image:

- **Click** → opens an N-shutter slitlet at the nearest operable
  shutter.
- **Click an open shutter** → closes it (and its siblings).
- **Double-click** → toggles a cyan visual highlight (not exported).
- **Shift-click** → moves the pointing centre to that location.
- **Mouse wheel** → zoom.

Loaded a catalog? vMPT auto-tags every opened slitlet with the
catalog source ID whose footprint falls inside any opened shutter.
The status bar (bottom of the page) names the matched source.

## 5. Or run the optimizer

In the **Pointing** tab, click `Open optimizer…`. Choose a method
(Democracy / Meritocracy / Hierarchy — see
[MSA pointing optimizer](optimizer.md)), set your search radius and
optional collision protection, then `Run optimization`. The top-10
solutions appear in a results table; click `Apply #N` to set the
pointing and auto-open slitlets for the matched targets.

Planning **two exposures**? Bump *Number of configs* to 2 (beneath the
optimizer) — one Run then returns a complete 2-config plan, and *Max
configs per source* = 1 keeps the second config from re-observing the
first's targets. See [Multiple configurations](multiconfig.md).

## 6. Save / share

In the **MPT** tab:

- `Save session` → writes a self-contained bundle (`MPT_plan.json` +
  `vMPT_workspace.json` + the eMPT three CSVs) to a directory you
  choose. Mail or git that directory; the recipient does
  `Load session` and gets your exact picks.
- `Export eMPT bundle` → same writer, into a fresh timestamped
  `exports/` subfolder.

See [Exporting & sharing](exporting.md) for the on-disk format.
