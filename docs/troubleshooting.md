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

## A shutter is failed-closed in APT/MPT but operable in vMPT (stale operability)

vMPT reads the failed-/stuck-shutter map from the newest
`jwst_nirspec_msaoper_*.json` in your local **CRDS cache**
(`~/crds_cache/references/jwst/nirspec/`). That file only updates when
something populates the cache — running the JWST pipeline, `crds sync`,
or the helper below. If your cache is behind the current operability,
vMPT will show shutters as operable that APT/MPT (which uses the live
operability) marks as failed-closed.

The startup log line tells you which file is in use:

```
[msa] Loaded operability from .../jwst_nirspec_msaoper_0017.json
```

To fetch the **current** MOS operability reference (the one the
operational CRDS context selects for `EXP_TYPE = NRS_MSASPEC` — the same
source APT/MPT uses) into your cache:

```bash
python scripts/update_operability.py    # needs network to jwst-crds.stsci.edu
```

Then **restart vMPT** — operability is read once at startup, so a running
server keeps the old map until relaunched. (CRDS reference *numbers* are
delivery order, not operability date; the script resolves the correct one
by USEAFTER for today, so you don't have to guess.)

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
