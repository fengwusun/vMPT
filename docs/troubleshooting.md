# Troubleshooting

Common issues you may hit on a fresh install or when moving vMPT
between machines, with the fix in one line.

## "Load Abell 370 example" / "Load RXCJ0600 example" do nothing

The two example datasets are NOT bundled in the pip wheel (they'd
add ~70 MB of FITS/JPG data that most users never need). Run:

```bash
vmpt examples download
```

once per machine. Since v1.3.1 the tarball drops into
`~/.vmpt/examples/` (writable on all systems, persistent across
sessions). The in-app buttons search there first, then your dev
source-checkout, then your current directory — so you can click
either button straight after the download with no restart.

If the download itself fails, the message gives the GitHub release
URL it tried; you can `curl -L` the tarball by hand and extract
into `~/.vmpt/examples/`.

For a non-standard location (cluster home, shared mount), set
`VMPT_EXAMPLES_DIR=/some/path` before running the download AND
before launching the app.

## Port 5006 / 5099 already in use

```
ERROR: Cannot listen on port 5006 (in use).
```

Either kill the conflicting process or pick a different port:

```bash
vmpt --port 5010
```

The same flag works for `run.sh` from a source checkout.

## "Web socket origin error" in the browser console

Happens when you open vMPT on a different host/port than Bokeh was
told to accept. The `vmpt` console script auto-passes
`--allow-websocket-origin localhost:<port>`; if you're running on a
remote machine and tunneling, add the visible host yourself:

```bash
vmpt --port 5099 --allow-websocket-origin myserver.example.edu:5099
```

## The MSA quadrants disappear when I load a second image

Happens when the first image set a pointing centre and the second
image is in a completely different field. vMPT keeps your pointing
centre on image swap by design (so APT roll edits don't reset
you). Either:

- Click the **Load … example** button — these are HARD resets and
  always re-centre the MSA on the image.
- Or in the **Pointing** tab, clear the RA/Dec inputs and reload
  the image: an empty RA/Dec auto-centres.

## "FITS WCS not picked up" / image loads as a square greyscale

Two flavours:

1. Your image is a JPG/PNG without a sidecar `.fits` next to it →
   vMPT can't derive sky coordinates. Either use the **JPG + sidecar
   FITS** path (point the Sidecar input at a small FITS with the
   correct WCS), or convert the JPG to a real FITS first.

2. Your FITS lacks a `CTYPE1/CTYPE2 = RA---TAN/DEC--TAN` pair →
   `astropy.wcs` warns and rejects. Re-create the FITS with a
   conforming WCS header, or open it via Astropy and re-save.

## Catalog loads but no targets show on the image

- Open **Settings → Layers** and check **Show catalog targets**.
- Click `Edit catalog…` — the popover shows the parsed columns.
  vMPT uses case-insensitive substring matching for `ra`, `dec`,
  `id`, etc., but if your column header is something exotic like
  `RA_J2000_deg`, manual rename in the editor (or
  pre-rename in the CSV) fixes it.
- Targets outside the current image's WCS footprint are silently
  clipped — pan/zoom out to see them. If your image WCS covers a
  tiny area but the catalog spans the sky, that's a real mismatch.

## The orange / pink / purple overlay is missing for some dispersers

The detector-pixel collision check needs the per-shutter x/y data
shipped in `data/dispersion_cutoffs.npz`. The wheel includes all 9
supported (disperser, filter) combos by default. If you ship a
patched build and a combo is missing, the live overlay
transparently falls back to the V2-distance heuristic (and may
slightly over-flag at filter edges — see [Spec-overlap colours](spec_overlap.md)).

Check:

```python
python -c "import numpy as np, vmpt; \
  f=np.load(vmpt.__path__[0]+'/data/dispersion_cutoffs.npz'); \
  print('combos with xy data:', sorted({k.rsplit('_',2)[0] for k in f.files if k.endswith('_x_lo_nrs1')}))"
```

A complete wheel reports 9 combos; missing entries are the ones
that fall back.

## Operability map looks empty, or a shutter is failed in APT/MPT but operable in vMPT

vMPT resolves the failed-/stuck-shutter map in this order:

1. **Live CRDS reference** — the current MOS `jwst_nirspec_msaoper_*.json`
   in your CRDS cache (`$CRDS_PATH`, else `~/crds_cache/references/jwst/nirspec/`).
   vMPT also refreshes it from `https://jwst-crds.stsci.edu` on startup
   (best-effort, ~8 s budget). This is the same source APT/MPT uses.
2. **Bundled snapshot** — `vmpt/data/msaoper_fallback.npz`, shipped with vMPT.
   Used when no CRDS reference is reachable: a fresh `pip install`, no network,
   a firewall, an unwritable cache, or a slow first fetch. It always loads, but
   **may lag** the current operability.

The startup log tells you which is in use (and case 2 also shows a banner in
the app's status bar):

```
[msa] Loaded operability from .../jwst_nirspec_msaoper_0017.json      # live CRDS — current
[msa] Using the BUNDLED operability snapshot (…); it ships with vMPT…  # snapshot — may lag
```

If you're on the bundled snapshot (or APT/MPT marks a shutter failed that vMPT
shows operable), get the **current** reference and restart:

- Make sure the environment can reach CRDS — a **writable `CRDS_PATH`**
  (default `~/crds_cache`) and **network access** to `https://jwst-crds.stsci.edu`.
  vMPT auto-fetches the correct MOS reference on startup, so usually just fixing
  network/cache and **restarting** is enough. (Behind a proxy/offline cluster?
  copy a populated `~/crds_cache` from a networked machine, or run `crds sync`.)
- From a repo checkout you can force it now:
  `python scripts/update_operability.py` (maintainers: add `--bundle` to also
  refresh the shipped snapshot).

Operability is read once at startup, so **restart vMPT** to pick up a new
reference. (CRDS reference *numbers* are delivery order, not operability date;
vMPT resolves the correct one by USEAFTER for today, so you don't have to guess.)

## I see a single row "sticking out" of an N-shutter slitlet's orange band

That used to happen up to v1.3.0 because of how tilt was rendered
(see CHANGELOG v1.3.1). Since v1.3.1 the orange band is clamped to
exactly N+2 rows across the entire dispersion direction. If you
still see it on v1.3.1+, please open a GitHub issue with the
shutter coordinates (q, s, d) and the disperser/filter.

## The optimizer runs forever / never returns

Two usual causes:

- **Search radius too wide**: `ΔRA = 60′` × `ΔDec = 60′`
  × `ΔPA = 360°` × grid step = > 10⁶ points. Tighten to ±5′ on
  pointing and ±5° on PA for first runs.
- **Collision protection on with very few targets**: the constraint
  rejects nearly every candidate. Either expand the catalog, lower
  the priority threshold, or temporarily uncheck "Protect spectra
  from collisions" to confirm the optimizer is finding ANY valid
  pointing.

Cancel via the modal's **×** corner button (history is preserved
even after cancel).

## "Cannot find aperture NRS_FULL_MSA" / pysiaf error on import

```
KeyError: 'NRS_FULL_MSA'
```

You're on a pysiaf release that predates the JWST PRD update that
introduced this aperture. Update:

```bash
conda update pysiaf      # or `pip install --upgrade pysiaf`
```

vMPT pins `pysiaf >= 0.21` in `requirements.txt`; older venvs may
have an unconstrained install.

## "Save session" writes nothing / overwrite dialog disappears

Most often a filesystem permission issue (clicking Save with a
path inside `/usr/local` etc.). vMPT writes to the directory you
pick — make sure your user owns it. The Browse popup defaults to
your home directory.

## I want to undo a long picking session

**Settings → Actions → Undo last** reverts the most recent
open/close. History is 50 deep. The deletion is destructive past
the 50th step — save the session beforehand if you'll need to
return to an earlier state.

## Something broke on a new vMPT release

Save a regression test case before you upgrade:

```bash
pip install jwst-vmpt==1.3.0    # downgrade to known-good
```

…then open a GitHub issue with the saved `vMPT_workspace.json`
and the (disperser, filter, V3 PA) values. The session file is
self-contained and reproduces the exact state.

## A huge FITS looks coarse when zoomed out

That's by design, and it sharpens when you zoom in. Since v1.8.0 vMPT
**memory-maps** FITS files and shows a downsampled base view (~4000 px on the
long axis) so GB-scale mosaics open in well under a second instead of hanging the
tab. **As you zoom in, vMPT automatically re-reads just the visible region at
higher resolution — up to the file's exact native pixels at full zoom.** The
refresh is debounced, so it sharpens a moment after you stop panning/zooming
(you may see a brief coarse-then-sharp step on deep zoom — that's the on-demand
crop loading). No setting to tune; the WCS is rescaled at every level so catalog
targets and overlays stay aligned.

## "Load DS9 regions" / "Load contours" shows nothing

A few things to check:

- **An image must be loaded first.** Region/contour shapes are projected through
  the image WCS, so they're queued (with a status note) until an image is
  present, then drawn when it loads.
- **The layer is toggled off.** Tick *Show DS9 regions* / *Show contours* under
  **Settings → Layers**.
- **Contour coordinate frame.** A DS9 `.ctr` file records its frame (e.g.
  `icrs` or `image`) and vMPT reads it automatically. A bare `.con` with no
  frame line defaults to **sky**; if such a file is actually in image pixels its
  contours will land far off — re-export it from DS9 as `.ctr` (which carries
  the frame) so the coordinate system is detected correctly.
- **`regions` package missing.** DS9 `.reg` parsing needs the
  [`regions`](https://astropy-regions.readthedocs.io/) package (a vMPT
  dependency since v1.8.0). On an old environment, `pip install regions`.
- **Off-field overlay.** If the shapes load but you can't see them, they may be
  outside the current view — the region coordinates may not match your image's
  field. Double-check the file is for this field.
