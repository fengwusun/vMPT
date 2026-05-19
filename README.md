# vMPT — visual MSA Planning Tool

Interactive Bokeh app for hand-picking JWST/NIRSpec MSA shutter
configurations on an image of the target field. Lets you (and your
collaborators) pick shutters one at a time, see the spectral
conflicts in real time, and export a bundle that loads into APT
and/or the [eMPT pipeline](https://github.com/esdc-esac-esa-int/eMPT_v1).

vMPT mirrors the workflow of MPT and eMPT but **without** the
automated optimization step — you keep full control. Save the state
as a JSON file, send it to a collaborator, and they pick up where
you left off.

![status](https://img.shields.io/badge/tests-60%2B%20passed-brightgreen)
![python](https://img.shields.io/badge/python-3.11-blue)
![license](https://img.shields.io/badge/license-MIT-blue)

---

## Installation

vMPT is a local-only tool: it runs on your machine, files stay on
your disk, computation uses your local Python. There are three
install paths — pick whichever matches your environment.

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
pytest tests/    # 60+/60 should pass; ~7 seconds
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

2. **Load an example** from the **Image** tab (left sidebar):
   - **"Load Abell 370 example"** — a 42 MB three-band FITS at
     `example_a370/a370_f182m_f200w_f210m.fits` (JWST/NIRCam F182M +
     F200W + F210M). Includes a target catalog and an APT MPT plan
     from GTO-1208 you can load to see the full picking workflow.
   - **"Load RXCJ0600 example"** — a 17 MB JPG + WCS-sidecar pair
     (`example_r0600/*`). Demonstrates the JPG+sidecar workflow.
     Includes a ~28k-source target catalog.

   Both examples are committed in the repo (no extra download
   needed). A full-page spinner overlays the canvas during loads so
   you know the app is busy.

3. **Aim the MSA** in the **Aim** tab:
   - The pointing center auto-fills to the image center (unless a
     plan was loaded first — in that case the plan's RA/Dec is kept).
   - Drag the **V3 PA** slider or type into the **APA** box. The
     spinner appears during recomputation.
   - Optional: enter a date and click **"Compute allowed V3 PA"** —
     queries `jwst_gtvt` and reports the valid V3 PA window for that
     date.

4. **Pick shutters** in the **Pick** tab + on the image:
   - **Click on the image** → opens an **N-shutter slitlet** at the
     nearest operable shutter. Choose N in the **N-shutter slitlet**
     dropdown:
       - `N=1` → only the clicked shutter
       - `N=2` → click + one row lower (lower y on the detector)
       - `N=3` → centred 3-shutter slitlet (the standard for MOS)
       - `N=5` → centred 5-shutter slitlet
   - **Click an open shutter** → closes it AND its slitlet siblings.
   - **Double-click** → toggles a cyan highlight (a visual flag, not
     exported).
   - **Shift-click** → moves the pointing center to that location.
   - **Wheel** → zoom both axes equally.
   - **Drag** → pan.

   If a catalog is loaded, vMPT auto-tags the slitlet with the
   catalog source ID whose footprint falls inside any opened shutter.
   The status bar shows which source was matched.

5. **Save and export** in the **MPT** tab:
   - **Save session** → writes the bundle (see "Bundle output" below)
     to a single chosen directory. Share the directory with a
     collaborator to hand off the work.
   - **Load session** → restores any vMPT bundle (point at either
     `MPT_plan.json` or `vMPT_workspace.json` — both work; the
     sibling is auto-discovered).
   - **Export eMPT bundle** → same writer as Save session, but into
     a fresh timestamped subfolder of `exports/`.

---

## Color legend

What each color means on the figure:

| Color | Meaning |
|---|---|
| **Dodgerblue rectangles** | The four MSA quadrant outlines |
| **Gold polygons** | The 5 NIRSpec fixed slits (always visible) |
| **Lime cross** | Current pointing center (shift-click to move) |
| **Silver-edge boxes** (α=0.2) | Operable, **unaffected**, ready-to-pick shutters (toggle "Show operable shutters" in Pick → Layers) |
| **Dark-red thick outline** | Stuck-open shutter (always visible) |
| **Red-filled (#ff8888)** | User-opened shutter |
| **Cyan edge** | Highlighted shutter (double-click marker) |
| **Orange fill** (α=0.10, stackable) | Spectral-conflict warning — operable shutters whose spectra would overlap on the detector with an open or stuck-open shutter's. Darker orange = multiple opens contribute. ±1 row from each dispersion source; cross-quadrant via NRS1 (Q1↔Q3) and NRS2 (Q2↔Q4) detector pairing. |
| **Yellow circles** | Catalog targets (toggle "Show catalog targets" in Pick → Layers) |

Failed-closed shutters are not drawn at all — they don't exist for
the user's purposes.

---

## Loading your own data

### Image: FITS

In the **Image** tab, paste the absolute path into "FITS path (local)"
and hit Enter (or use the **Browse…** button). First HDU with image
data is auto-selected; the WCS comes from that HDU's header.

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

When you open a shutter that contains a catalog source, the slitlet
is auto-tagged with that source's ID. Matching follows APT's
*[Unconstrained](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-planner)*
Source Centering rule — a source matches the shutter whose **full
pitch cell** (≈0.27″×0.53″) contains its centre, so sources sitting
behind the MSA bars still get matched (to whichever neighbouring
shutter is nearest). The status bar names the matched source.
Slitlets with no real source get a synthesized entry at the slitlet's
centre at export time, tagged in the catalog's `Label` column as
`vMPT_synth`.

A small example catalog lives at `tests/fixtures/tiny_catalog.csv`.

---

## Loading an APT plan

The **MPT** tab also reads plans straight from APT exports:

- **Load plan from JSON** — paste the path to an APT MPT JSON
  (each `configs[]` entry becomes a selectable plan; pointing, V3 PA,
  disperser, and slitlets are applied on click).
- **Load shutter CSV** — the open-mask CSV that APT/MPT/eMPT write.
- **Fetch / open .aptx** — point at a local `.aptx` archive **or**
  enter a JWST program ID; vMPT downloads the latest archive from
  STScI's public proposal-info endpoint, lists the embedded MPT
  plans, and lets you load any one. (Some programs may not be
  publicly available yet — the fetch will report 404.)

If you load a plan first and then an image, vMPT keeps the plan's
pointing (instead of recentering on the image).

---

## Bundle output

When you click **Save session** or **Export eMPT bundle**, vMPT
writes a directory with six files. The prefixes telegraph the role:

| File | Role | Format |
|---|---|---|
| **`MPT_plan.json`** | APT MPT plan — load via APT MOS → MSA Planner | APT MPT JSON, matches the reference schema field-for-field |
| **`<catalog>.cat`** | APT-importable Target List — name matches the user's catalog (or `MPT_catalog.cat` if none was loaded) | ASCII, tab-separated, `#`-header with the JDox-recognized labels (`ID`, `RA`, `DEC`, `Weight`, `Primary`, `Label`). The `Label` column carries `real` or `vMPT_synth` so you can tell which rows came from your input catalog. |
| **`vMPT_workspace.json`** | vMPT-only state — per-shutter `target_id`+`role`, highlighted set, image / sidecar / catalog paths, slitlet height, exact V3 PA | vMPT-internal JSON |
| **`eMPT_observed_targets.cat`** | eMPT-style target list | eMPT format |
| **`eMPT_pointing_summary.txt`** | eMPT-style pointing summary | eMPT format |
| **`eMPT_shutter_mask.csv`** | 730×342 MSA mask, byte-compatible with eMPT's writer (`shutter_routines_new.f90`) | eMPT format |

### Loading the bundle into APT

1. **Import the target list**: APT → *Targets → Target Lists → Import…*
   → select the `<catalog>.cat` file. APT names the list after the
   file stem; that stem matches `catalog.name` inside `MPT_plan.json`.
2. **Load the plan**: APT → MOS template → *MSA Planner → Load Plan*
   → select `MPT_plan.json`. APT pairs the plan with the Target List
   imported in step 1.

### Loading the bundle back into vMPT

Point **Session load path** at either `MPT_plan.json` or
`vMPT_workspace.json`. vMPT auto-discovers the sibling and restores:
pointing, V3 PA, disperser/filter, slitlet height, every open
shutter with its target_id+role, the highlighted set, and (if the
image still exists on disk) the image + WCS sidecar.

---

## Collaborating on a target list

```
You                                       Collaborator
───                                       ────────────

1. Open vMPT, load image + catalog
2. Pick shutters
3. Save session  ──── bundle dir ───────> Load session
                                          (vMPT loads the same image +
                                           catalog + picks)
                                          Add / remove / adjust picks
                                          Save session
8.  Load session <──── bundle dir ────────
9.  Continue picking
...

When done:
   Export eMPT bundle  →  6 files, ready for APT + eMPT
```

The workspace JSON contains paths to the image and catalog. For
those to resolve on the collaborator's machine, use a shared mount
(Dropbox / Drive / network share / `git lfs`-tracked data folder),
or edit `image_path` / `catalog_path` / `wcs_sidecar_path` in
`vMPT_workspace.json` before each handoff. The MPT plan JSON itself
carries no file paths — it's safe to share standalone.

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
"Browse…" file picker — the file is read directly from disk, no
WebSocket size limit applies.

**APT can't find the catalog when loading my plan**
Import the matching `<catalog>.cat` file in APT first (Targets →
Target Lists → Import). The file stem and `MPT_plan.json`'s
`catalog.name` are aligned by the writer, but you still need to
load the target list once on the APT side.

**`jwst_gtvt` query takes forever the first time**
First call downloads JWST's ephemeris file (~30 MB). Subsequent calls
in the same session are fast.

**`session.json` from an old vMPT version doesn't load**
Pre-1.1 sessions (flat top-level `open_shutters`) still load on the
legacy path. Pre-1.4 sessions used filenames `session_MPT_plan.json`
+ `vmpt_workspace.json`; those are recognised as fallback names and
still work.

---

## Tool architecture

```
app/
├── main.py            Bokeh server entry; UI wiring; on_tap / on_export
├── coords.py          V2/V3 ↔ RA/Dec transforms (pysiaf-backed)
├── msa.py             MSA shutter grid + CRDS operability loader
├── wavelengths.py     Analytic per-grating dispersion + cutoffs
├── image_io.py        FITS + JPG-with-sidecar loaders
├── catalog.py         CSV/ASCII/FITS catalog reader
├── empt_io.py         eMPT-format + MPT-catalog writers
├── session_io.py      Bundle save/load (MPT plan + workspace sidecar)
├── mpt_io.py          APT MPT JSON parser + .aptx archive reader
├── static/            favicon.svg
└── templates/         index.html (injects the favicon as a data URI)

data/
└── nirspec_msa_v2v3.npz   Per-shutter V2/V3 coordinates (4×171×365)

tests/                 pytest suite (60+ tests, ~7 s)
example_a370/          Abell 370 cluster FITS (44 MB)
example_r0600/         RXCJ0600 JPG + sidecar (240 MB)
exports/               default output dir for bundles
```

### Performance

`refresh_overlays` runs in ~10 ms for the light path (MSA outline +
pointing handle only, during slider drags) and ~70 ms for the full
path with operable + spec-overlap layers on. The hot path is pure
numpy: precomputed V2/V3 offsets for all 249,660 shutters, a single
WCS inverse-Jacobian per refresh, and two matmuls (rotation by PA,
then sky→pixel).

The operable-shutter layer is filtered to *unaffected, ready-to-pick*
shutters only (excludes user-opens, stuck-opens, spec-overlap rows),
keeping the rendered polygon count well below the `MAX_OPERABLE_RENDER`
cap (10,000) at the typical "looking at one quadrant" zoom level.

---

## Known limitations

- **V2 dispersion calibration**: per-disperser spectrum extents are
  approximated; PRISM is calibrated against eMPT's `prism_sep.dat`
  (35″ V2 half-extent), M/H gratings are approximations
  (200″ and 500″ respectively). For research-quality numbers,
  replace with JDox-sourced or CRDS-derived constants.
- **`plannerSpecification` in MPT plan JSON**: written with sensible
  defaults (matching the reference G395H_F290LP plan schema field
  for field) but the dither / search-grid parameters don't reflect
  any vMPT internal state — APT uses them as starting values for
  re-planning if you choose to.
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
- **NIRSpec MOS / MPT**: [JDox MPT page](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template)
- **MPT catalog format**: [JDox MPT Catalogs](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-apt-templates/nirspec-multi-object-spectroscopy-apt-template/nirspec-mpt-catalogs)
- **MSA operability**: STScI CRDS `jwst_nirspec_msaoper_*.json` (auto-loaded if `CRDS_PATH` is set)
- **jwst_gtvt** (visibility): [GitHub](https://github.com/spacetelescope/jwst_gtvt)

## Example data — attribution

The two example fields shipped under `example_a370/` and
`example_r0600/` are JWST/NIRCam images of well-studied lensing
clusters. The image files were prepared for use as vMPT examples
(RGB stretches; the R0600 JPG was re-encoded at JPEG quality 85 to
keep the repo small — dimensions and WCS are unchanged from the
science-grade version). The accompanying target catalogs and APT
MPT plans are research products from real JWST programs.

If you use the example data for anything beyond trying vMPT itself,
please cite the originating program / data release directly — vMPT
just ships them as a starting point.

## License

MIT. See [LICENSE](LICENSE).

## Citation

If vMPT helps you plan an observation that ends up in a paper, a
mention is appreciated. The export-bundle format is calibrated
against eMPT; please cite Bonaventura et al. 2023 if you use the
`eMPT_*` files.
