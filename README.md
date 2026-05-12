# vMPT — visual MSA Planning Tool

Interactive Bokeh app for hand-picking JWST/NIRSpec MSA shutter
configurations on an image of the target field. Lets you (and your
collaborators) pick shutters one at a time, see the spectral
conflicts in real time, and export a bundle that loads into APT.

It mirrors the workflow of MPT and
[eMPT](https://github.com/esdc-esac-esa-int/eMPT_v1) but **without**
the automated optimization step — you keep full control. Save the
state as a JSON file, send it to a collaborator, and they pick up
where you left off.

![status](https://img.shields.io/badge/tests-50%20passed-brightgreen)
![python](https://img.shields.io/badge/python-3.11-blue)
![license](https://img.shields.io/badge/license-MIT-blue)

---

## Installation

vMPT is a local-only tool: it runs on your machine, files stay on
your disk, computation uses your local Python. There are two install
paths — pick whichever matches your environment.

### Option A — STScI's `stenv` (recommended for JWST users)

If you already use the [STScI JWST/HST pipeline environment](https://stenv.readthedocs.io/),
most dependencies are already present and you only need to add Bokeh
and `jwst_gtvt`.

```bash
git clone https://github.com/fengwusun/vMPT.git
cd vMPT
conda activate stenv
pip install bokeh jwst_gtvt
./run.sh
```

The browser should open at `http://localhost:5006/app`.

### Option B — fresh conda env

If you don't have `stenv` or want a clean environment:

```bash
git clone https://github.com/fengwusun/vMPT.git
cd vMPT
conda create -n vmpt python=3.11
conda activate vmpt
pip install -r requirements.txt
./run.sh
```

### Option C — plain pip (no conda)

```bash
git clone https://github.com/fengwusun/vMPT.git
cd vMPT
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
./run.sh
```

### Verify the install

```bash
pytest tests/    # 50/50 should pass; ~6 seconds
```

If everything's green, the tool is ready. If `pytest` complains about
a missing module, check that you're in the right environment.

---

## First-time use — two-minute tour

The repo ships with two example fields you can load with one click.

1. **Start the server**:
   ```bash
   ./run.sh
   ```
   Browser opens at `http://localhost:5006/app`. If it doesn't, open
   that URL yourself.

2. **Load an example** from the **Load** tab (left sidebar):
   - **"Load Abell 370 example"** — a 44 MB three-band FITS at
     `example_a370/a370_f182m_f200w_f210m.fits`. Fastest to try.
   - **"Load RXCJ0600 example"** — a 251 MB JPG + WCS-sidecar pair
     (`example_r0600/*`). Demonstrates the JPG+sidecar workflow.

3. **Aim the MSA** in the **Aim** tab:
   - The pointing center auto-fills to the image center. Drag the
     **V3 PA** slider to rotate the MSA pattern.
   - Optional: enter a date and click **"Compute allowed V3 PA"** —
     queries `jwst_gtvt` and reports the valid V3 PA window for that
     date. Snaps the slider to the nominal angle.

4. **Pick shutters** in the **Pick** tab + on the image:
   - **Click on the image** → opens the nearest operable shutter as
     a 3-shutter slitlet (the clicked shutter + the two adjacent rows
     for nod-and-shuffle).
   - **Click an open shutter** → closes it.
   - **Double-click** → toggles a cyan highlight (a visual flag, not
     exported).
   - **Shift-click** → moves the pointing center to that location.
   - **Wheel** → zoom both axes equally.
   - **Drag** → pan.

5. **Save and export** in the **Save** tab:
   - **Save session** → writes a JSON snapshot of the full state.
     Share this file with a collaborator to hand off the work.
   - **Export eMPT bundle** → writes the three files APT needs into
     a timestamped subfolder of `exports/`.

---

## Color legend

What each color means on the figure:

| Color | Meaning |
|---|---|
| <span style="color:dodgerblue">**Dodgerblue rectangles**</span> | The four MSA quadrant outlines |
| <span style="color:gold">**Gold polygons**</span> | The 5 NIRSpec fixed slits (always visible) |
| <span style="color:lime">**Lime cross**</span> | Current pointing center (shift-click to move) |
| Faint white grid | Operable shutters (toggle in Pick → Layers) |
| <span style="color:red">**Red edge**</span> | Stuck-open shutter (always visible) |
| <span style="color:red">**Red filled**</span> | Open (you-selected) shutter |
| <span style="color:cyan">**Cyan edge**</span> | Highlighted shutter (double-click marker) |
| <span style="color:orange">**Orange tint**</span> | Spectral-conflict warning — these shutters' spectra would collide with an open shutter's on the detector. Darker orange = multiple opens contribute. |
| <span style="color:yellow">**Yellow circles**</span> | Catalog targets (loaded in Load tab) |

Failed-closed shutters are not drawn at all — they don't exist for
the user's purposes.

---

## Loading your own data

### Image: FITS

In the **Load** tab, paste the absolute path into "FITS path (local)"
and hit Enter. First HDU with image data is auto-selected; the WCS
comes from that HDU's header.

### Image: JPG + sidecar FITS

For fields where you have a pretty RGB JPG but the WCS lives in a
separate FITS header (this is what tools like
[fitsmap](https://github.com/ryanhausen/fitsmap) produce), put the
WCS-only FITS path in "Sidecar FITS path" and the JPG path in "JPG
path". Order matters — set the sidecar first.

The JPG can be tens of millions of pixels; vMPT downsamples to ≤6000
on the longest edge and rescales the WCS accordingly.

### Catalog (optional)

CSV, whitespace-ASCII, or FITS table. Required columns
(case-insensitive): **ID, RA, DEC**. Optional and used if present:
`priority`/`Pr`, `mag`/`mag_F444W`/`F444W_mag`, `z`/`zspec`/`zphot`,
`label`/`name`. Yellow circles appear on the image; the **Pick** tab
has compact text inputs to filter by priority class or magnitude.

A small example catalog lives at `tests/fixtures/tiny_catalog.csv`.

---

## Collaborating on a target list

The intended team workflow:

```
You                                    Collaborator
───                                    ────────────

1. Open vMPT, load image + catalog
2. Pick shutters
3. Save session  ──── session.json ──> Load session
                                       (vMPT loads the same image +
                                        catalog + picks)
                                       Add / remove / adjust picks
                                       Save session
8.  Load session  <── session.json ────
9.  Continue picking
…

When done:
   Export eMPT bundle  →  3 files for APT
```

The session JSON is small (a few KB) and contains:
- Pointing RA / Dec / V3 PA
- Disperser + filter
- Every open shutter with its `(q, s, d)` triple and target ID
- Highlighted (flagged) shutters
- Paths to the image and catalog (so the receiver loads the same
  files automatically — **paths must be reachable on their disk**)

To make the path part work, your team should either:
- Use a shared filesystem mount (Dropbox / Google Drive / network
  share / `git lfs`-tracked data folder), so the paths inside the JSON
  resolve identically on every machine, **or**
- Open the JSON manually and edit the `image_path` / `catalog_path`
  fields before each collaborator loads it.

Bundled session: when you click **Export eMPT bundle**, vMPT also
writes `session.json` *inside* the same bundle directory. Hand someone
the whole folder and they can both reload the picks and re-import
into APT without juggling separate files.

---

## Exporting to APT

The **Export eMPT bundle** button writes three files into a
timestamped subfolder of your export directory (default
`./exports/empt_bundle_YYYYMMDDTHHMMSS/`):

| File | What APT does with it |
|---|---|
| `observed_targets.cat` | Whitespace-separated. *Form Editor → Targets → Import MSA Source Catalog* (Whitespace Separated format). |
| `pointing_summary.txt` | Free-form text. **Copy** the `PA_AP`, RA, and Dec values into APT's *MSA Planner → Search Grid* (set search box width/height to 0 and N=1). Click **Generate Plan**. |
| `shutter_mask.csv` | 730 × 342 grid. For each nod row in the generated plan: *Edit Configuration → Edit → Import CSV*. Overwrites APT's auto mask with vMPT's. |

The `shutter_mask.csv` is byte-compatible with eMPT's format — the
writer was reverse-engineered from `reference_files/shutter_routines_new.f90`
in the eMPT repo and round-trips identical against the bundled
`trial_00_ref/` example.

---

## Troubleshooting

**`bokeh: command not found`**
You're not in the right Python environment. Activate the env where
you ran `pip install`:
```bash
conda activate stenv      # or vmpt, depending on which option you used
```

**Port 5006 already in use**
Another Bokeh process is running. Find and kill it:
```bash
pkill -f "bokeh serve"
```
Or pass `--port 5007` to `bokeh serve`.

**Image upload fails or stops with "No 2D image HDU found"**
The Bokeh WebSocket has a default 20 MB cap; large uploads get
truncated. `./run.sh` raises the cap to 500 MB, but the right move
is to **use the path input** ("FITS path (local)") instead of the
"Choose File" upload — the file is read directly from disk, no
WebSocket size limit applies.

**Clicking a shutter opens the wrong one**
This was a Dec≠0 WCS-Jacobian bug; fixed in commit `4682571`. Pull
the latest:
```bash
git pull
pkill -f "bokeh serve"
./run.sh
```

**`jwst_gtvt` query takes forever the first time**
First call downloads JWST's ephemeris file (~30 MB). Subsequent calls
in the same session are fast.

**Wavelength tooltip shows λ > 5.3 μm for PRISM**
Fixed in `b96f126` — pull and restart. The cutoff is now clamped to
the grating's intrinsic range.

---

## Tool architecture

```
app/
├── main.py            Bokeh server entry; UI wiring
├── coords.py          V2/V3 ↔ RA/Dec transforms (pysiaf-backed)
├── msa.py             MSA shutter grid + CRDS operability loader
├── wavelengths.py     Analytic per-grating dispersion + cutoffs
├── image_io.py        FITS + JPG-with-sidecar loaders
├── catalog.py         CSV/ASCII/FITS catalog reader
├── empt_io.py         eMPT-format export writers
└── session_io.py      JSON save/load of picking session

data/
└── nirspec_msa_v2v3.npz   Per-shutter V2/V3 coordinates (4×171×365)

tests/                 pytest suite (50 tests, ~6 s)
example_a370/          Abell 370 cluster FITS (44 MB)
example_r0600/         RXCJ0600 JPG + sidecar (240 MB)
```

### Performance

`refresh_overlays` runs in ~10 ms when the operable-shutter layer is
toggled off (the default), and ~70 ms when it's on. The hot path is
pure-numpy: precomputed V2/V3 offsets for all 249,660 shutters, a
single WCS inverse-Jacobian computed at the pointing pixel per
refresh, and two matmuls (rotation by PA, then sky→pixel). PA slider
drag is real-time (light refresh during drag, full refresh on
release).

---

## Known limitations

- **V2 dispersion calibration**: per-disperser spectrum extents are
  approximated; PRISM is calibrated against eMPT's `prism_sep.dat`
  (35″ V2 half-extent), M/H gratings are approximations
  (200″ and 500″ respectively). For research-quality numbers,
  replace with JDox-sourced or CRDS-derived constants.
- **Bokeh single-session state**: opening the same server in two
  browser tabs lets picks bleed across them. Use one tab per user.
- **JPG/sidecar dimension mismatch**: a `>10 %` mismatch warns but
  doesn't refuse. Verify your sidecar.

---

## References

- **eMPT** (export format inspiration): Bonaventura et al. 2023,
  A&A 672 A40 — [arXiv:2302.10957](https://arxiv.org/abs/2302.10957) /
  [GitHub](https://github.com/esdc-esac-esa-int/eMPT_v1)
- **JWST PA conventions**: [JDox PA reference](https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-position-angles-ranges-and-offsets)
- **NIRSpec MOS / MPT**: [JDox MPT page](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/apt-workflow-articles/apt-mosaic-spectroscopy/mos-mode-msa-planning-tool)
- **MSA operability**: STScI CRDS `jwst_nirspec_msaoper_*.json` (auto-loaded if `CRDS_PATH` is set)
- **jwst_gtvt** (visibility): [GitHub](https://github.com/spacetelescope/jwst_gtvt)

## License

MIT. See [LICENSE](LICENSE).

## Citation

If vMPT helps you plan an observation that ends up in a paper, a
mention is appreciated. The export-bundle format is calibrated
against eMPT and please cite Bonaventura et al. 2023 if you use the
eMPT export path.
