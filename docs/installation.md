# Installation

vMPT is a local-only tool: it runs on your machine, files stay on
your disk, computation uses your local Python.

## With pip (recommended, **v1.2.2+**)

```bash
pip install jwst-vmpt
vmpt                                                             # opens the app
```

The console script `vmpt` accepts the same flags as the legacy
`run.sh`:

```bash
vmpt --port 5010                                                 # different port
vmpt --fits img.fits --catalog a.csv --catalog b.csv             # stack catalogs
vmpt --jpg img.jpg --wcs wcs.fits --catalog targets.csv          # JPG + WCS pair
vmpt examples download                                           # grab example_a370 + example_r0600
```

The wheel itself is ~20 MB (the required MSA grid + per-shutter
dispersion table). The two example datasets (~64 MB combined) are
fetched on demand via `vmpt examples download` and dropped into the
current directory; they are **not** bundled in the wheel.

## From source (developers)

```bash
git clone https://github.com/fengwusun/vMPT.git
cd vMPT
pip install -e .            # editable install, picks up local edits
./run.sh                    # same as `vmpt` after install
pytest tests/               # 139 passed, 4 skipped
```

### STScI's `stenv`

If you already use the [STScI JWST/HST pipeline environment](https://stenv.readthedocs.io/),
the heavy deps (`astropy`, `jwst`, `pysiaf`, …) are already there:

```bash
conda activate stenv
pip install jwst-vmpt       # or `pip install -e .` from a checkout
vmpt
```

## What's installed

After `pip install jwst-vmpt`:

| Thing | Where | Purpose |
|---|---|---|
| `vmpt/` package | site-packages | Bokeh app + optimizer + IO |
| `vmpt/data/nirspec_msa_v2v3.npz` | site-packages | MSA shutter centres in V2/V3 |
| `vmpt/data/dispersion_cutoffs.npz` | site-packages | per-shutter λ table |
| `vmpt` console script | `$PATH` | starts the Bokeh server |

## Verifying

```bash
vmpt --help                            # should print usage
python -c "import vmpt; print(vmpt.__file__)"
```

If `vmpt` is missing from `$PATH`, your venv's `bin/` directory
isn't being sourced — activate the venv or run with the absolute
path printed by `pip show -f jwst-vmpt | grep vmpt$`.
