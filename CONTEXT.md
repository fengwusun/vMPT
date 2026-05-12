# NIRSpec MSA Overlay & Planner — Project Context

Cold-start reference for the interactive NIRSpec-MSA-on-image **planning** tool. Think of it as a hand-driven MPT/eMPT: load an image, drop a target catalog on top, pick a pointing (RA, Dec, APA_V3), hand-pick which shutters to open, then export a configuration that can be loaded into APT.

---

## Goal

A local app that:

1. **Loads an image**:
   - FITS with embedded WCS (primary path), or
   - JPG/PNG paired with a **sidecar FITS** whose header supplies the WCS. The image array comes from the JPG, the WCS comes from the FITS header — the JPG and FITS must share the same pixel grid (or at least the same WCS solution).
2. **Loads a target catalog** (CSV/FITS with at least `ID, RA, DEC`; optional `priority, mag, z, label`) and overlays target markers on the image.
3. **Overlays the NIRSpec MSA** at a user-chosen pointing **(RA₀, Dec₀)** and **APA_V3** (= "NIRSpec PA"; aperture PA of the V3 axis, in degrees). One PA per session — no multi-PA stacking.
4. **Toggles** overlay layers:
   - Full MSA outline (4 quadrants).
   - All shutters (LOD-aware: only those inside view + margin).
   - **Operable** shutters only (failed-closed/failed-open masked out and colored differently).
   - **Open shutters** the user has selected (highlighted).
   - **Spectral trace** of each open slitlet for the chosen disperser/filter (shows the dispersed-light footprint on the sky so the user can avoid overlap and bad-pixel rows).
5. **Hand-picks shutters / slitlets**:
   - Click a target → app proposes the nearest operable shutter and a 3-shutter slitlet centered on it (NIRSpec MOS standard, with the target in the center shutter; the other two are used for nod-and-shuffle sky).
   - Click any shutter directly to add/remove it from the open set.
   - Detect and warn on spectral-trace conflicts between open slitlets.
6. For any shutter the user picks (or hovers), compute the **λ_blue, λ_red, λ_gap_lo, λ_gap_hi** of the dispersed spectrum given the chosen disperser/filter (G140M/F070LP, G140M/F100LP, G235M/F170LP, G395M/F290LP, G140H/F070LP, G140H/F100LP, G235H/F170LP, G395H/F290LP, PRISM/CLEAR).
7. **Exports a configuration** that can be ingested back by APT (or at least by hand-translation into APT's MPT). Concretely: pointing (RA, Dec, APA_V3), disperser/filter, and the full list of open shutters with their `(q, d, s)` indices plus the host target ID for each. See **MPT-style JSON** below.

---

## Coordinate plumbing (do not re-derive, copy from `footprint_emerald.ipynb`)

### Data file: `nirspec_msa_v2v3.npz`

Stored at `/Users/sunfengwu/jwst_cycle4/nirspec_msa_v2v3.npz` (copy to `data/` here on first use; ~4 MB).

```python
d = np.load("nirspec_msa_v2v3.npz")
v2_msa = d["v2_msa"]   # (4, 171, 365), float64, arcsec in V2
v3_msa = d["v3_msa"]   # (4, 171, 365), float64, arcsec in V3
```

Indexing convention (matches MPT JSON `slitlets` and the notebook):
- axis 0: **quadrant q − 1**  (q ∈ {1,2,3,4})
- axis 1: **shutter row s − 1**  (s ∈ {1…171})  — the "horizontal stripe" index
- axis 2: **shutter column d − 1**  (d ∈ {1…365}) — the "dispersion-direction" index along V2

Quadrant V2/V3 bounding boxes (sanity check):
- Q1: V2 ∈ [399.55, 533.73],  V3 ∈ [−503.04, −370.62]
- Q2: V2 ∈ [316.00, 448.20],  V3 ∈ [−406.62, −275.99]
- Q3: V2 ∈ [309.20, 442.30],  V3 ∈ [−583.70, −450.22]
- Q4: V2 ∈ [226.02, 357.41],  V3 ∈ [−486.74, −355.33]

Shutter pitch: **0.20″ along V2** (within a row), **~0.40″ along V3** (between rows). The shutter open area is **0.20″ × 0.46″** with a small bar between adjacent shutters in the row direction.

The MSA is rotated by **138.5°** within the V2/V3 frame (this is the constant from the notebook's `rot_matrix(138.5)` in `shutter_corners_v2v3`). That rotation maps "shutter-local x/y" → V2/V3 displacement so the corners come out aligned with the MSA grid.

### Apertures from pysiaf

```python
import pysiaf
siaf = pysiaf.Siaf('NIRSpec')
msa_ap = siaf['NRS_FULL_MSA']         # full MSA, gives V2Ref, V3Ref
msa1, msa2, msa3, msa4 = (siaf[f'NRS_FULL_MSA{i}'] for i in (1,2,3,4))
# msa_ap.V2Ref, msa_ap.V3Ref ≈ (378.563, -428.403) arcsec
# msaN.closed_polygon_points(to_frame='tel') → quadrant outline in (V2,V3)
```

### V2/V3 → RA/Dec at a given fiducial and PA_V3

Verbatim from the notebook (do not change the sign conventions):

```python
def rot_matrix(rotation=30):
    th = np.radians(rotation)
    c, s = np.cos(th), np.sin(th)
    return np.array(((c, -s), (s, c)))

def shutter_corners_v2v3(v2c, v3c, w=0.20, h=0.46):
    """4 corners of one shutter in V2/V3 (arcsec). Order: LL, LR, UR, UL."""
    return np.array([v2c, v3c]) + np.dot(
        np.array([[-w/2,  w/2, w/2, -w/2],
                  [-h/2, -h/2, h/2,  h/2]]).T,
        rot_matrix(138.5),
    )

def v2v3_to_radec(coord_c, pa_v3, corners_v2v3):
    """coord_c: SkyCoord fiducial (typically corresponding to msa_ap V2Ref,V3Ref).
       pa_v3:   APA_V3 in degrees. corners_v2v3: (N,2) in arcsec."""
    offsets = corners_v2v3 - np.array([msa_ap.V2Ref, msa_ap.V3Ref])
    offsets = np.dot(offsets, rot_matrix(pa_v3))      # rotate into sky frame
    return coord_c.spherical_offsets_by(
        offsets.T[0]*u.arcsec, offsets.T[1]*u.arcsec
    )
```

`spherical_offsets_by` takes (Δlon east, Δlat north). The sign convention is right because of the 138.5° pre-rotation: tangent-plane "east" lines up with rotated-V2.

### RA/Dec → image pixel

Use the image's `astropy.wcs.WCS` plus `astropy.wcs.utils.skycoord_to_pixel`. For JPG/PNG inputs, the app builds a WCS from user-supplied parameters (see "JPG mode" in `PLAN.md`).

---

## Export formats (APT ingestion)

We will emit **the same three artifacts eMPT (Bonaventura et al. 2023) writes**, because that pipeline is documented, has been used by us before, and gives us a known-working path into APT/MPT. Reference: arxiv.org/abs/2302.10957 and github.com/esdc-esac-esa-int/eMPT_v1 (`eMPT_user_guide_release_v1_doc_v1.1.pdf` §4).

The APT ingestion workflow (which we'll mirror in the export panel's instructions):

1. In APT **Form Editor → Targets**, "Import MSA Source Catalog" → load our `observed_targets.cat`.
2. In APT **MSA Planner** tab, paste `PA_AP`, RA, Dec from our `pointing_summary.txt` into the Search Grid panel; set search box width/height to 0 and Number of configurations to 1; **Generate Plan**.
3. For each nod row in the generated plan, **Edit Configuration → Edit → Import CSV** → load our `shutter_mask.csv`. This overwrites APT/MPT's auto-generated mask with our exact one.

### 1. `observed_targets.cat` — source catalog

Plain text, whitespace-separated, eMPT column convention:

```
# No   No_sub      No_cat    Pr    RA[deg]     Dec[deg]
   1        1        14170   1   53.1633910  -27.7756740
   2        1         8821   2   53.1641205  -27.7748813
   ...
```

- `No` — running index over rows actually placed in this configuration.
- `No_sub` — sub-index for multi-shutter sources; `1` for ordinary point sources.
- `No_cat` — original target ID from the user-supplied catalog.
- `Pr` — priority class (1 = highest; integer).
- `RA[deg]`, `Dec[deg]` — decimal degrees.

Only targets whose host slitlet is in the open-shutter set get written here.

### 2. `pointing_summary.txt` — pointings & PA

Free-form text, human-read, **copy-pasted manually into APT** (no auto-import). Must contain at minimum:

```
RA, Dec of Central Pointing:
 Nod 0:    189.1234567   62.2109876
Official Assigned APT/MPT roll angle:
 PA_AP:    273.000000
PA_V3:     <derived from APA_V3 = APA_aperture + V3IdlYAngle of the aperture used>
```

Even though our tool works in a single PA per session, we still emit a `Nod 1`/`Nod 2` block (with the same RA/Dec as `Nod 0`) so the file parses cleanly in APT's expected layout. The "PA_AP vs PA_V3" distinction matters for APT — confirm at M6 time which our session-level `apa_v3_deg` maps to.

### 3. `shutter_mask.csv` — MSA shutter configuration (the APT-loadable one)

CSV grid; one cell per shutter. Cell alphabet (from eMPT's `reference_files/shutter_routines_new.f90`):

- `x` — failed-closed (operability)
- `s` — failed-open (operability)
- `1` — functional, commanded closed
- `0` — commanded open ← the user's picks

The grid layout tiles the four quadrants into one matrix. Documented dimensions per the eMPT writer: **730 rows × 365 cols**, with rows 1–365 carrying the left half (Q1 + Q2 stacked) and rows 366–730 carrying the right half (Q3 + Q4 stacked). **Verify by reading `trial_00/m_pick_output/pointing_100/shutter_mask.csv` from the eMPT repo before finalizing the writer at M6** — exact tiling and row/column orientation must match byte-for-byte or APT's import will silently mis-place shutters.

Header line is literal:

```
# This CSV indicates which shutters should be open/closed on the MSA - created by ESA NIRSpec Team
```

We'll keep that header (drop-in compatibility) and add a second `#` comment line with our tool name and the source pointing.

### Internal session format (separate concern)

For session save/load inside our tool (round-tripping the user's picks without going through APT), we'll still use a JSON like:

```json
{
  "pointing":   {"ra_deg": 189.12, "dec_deg": 62.21, "apa_v3_deg": 273.0},
  "instrument": {"disperser": "G395M", "filter": "F290LP"},
  "open_shutters": [
    {"q": 2, "d": 200, "s": 86, "target_id": 123456, "role": "target"},
    {"q": 2, "d": 200, "s": 85, "target_id": 123456, "role": "sky"},
    {"q": 2, "d": 200, "s": 87, "target_id": 123456, "role": "sky"}
  ],
  "image_path": "/path/to/loaded.fits",
  "catalog_path": "/path/to/catalog.fits"
}
```

This is **internal only**; it doesn't go to APT. The three eMPT files above do.

---

## Target-catalog format

Minimum: `ID, RA, DEC` (decimal degrees). Optional columns the UI will use if present:

- `priority` / `Priority` — numeric, higher = more important; controls marker size/color.
- `mag_F444W` (or any single `mag_*` column) — numeric, controls marker color.
- `zspec`, `z` — numeric; if a line is configured, the app can show where it falls in the chosen disperser.
- `label` / `name` — string for hover/label text.

Accepted file types: CSV, ASCII (whitespace-separated with `astropy.io.ascii`), and FITS table.

---

## Operability mask

NIRSpec MSA shutters are individually failed-open / failed-closed / operable. STScI maintains the operability state per epoch:

- Reference files are distributed via **CRDS** as `jwst_nirspec_msaoper_*.json` (see jwst-docs MSA operability page).
- The JSON is a list of failure entries with fields `Q`, `x` (= column d), `y` (= row s), `state` (e.g. "open", "stuck closed", "stuck open").

Action item: fetch the latest `msaoper_*.json` once into `data/` and provide a loader that returns a `(4,171,365)` boolean array `operable[q-1, s-1, d-1]`. Keep both "failed-open" and "failed-closed" reasons in a parallel array so the UI can color them differently.

Until then, treat all shutters as operable.

---

## Wavelength cutoffs per shutter

NIRSpec is a slit spectrograph: the location of a shutter on the MSA determines the wavelength range that lands on each detector (NRS1, NRS2), with a physical gap between them. Computing **λ_blue, λ_gap_lo, λ_gap_hi, λ_red** for each (shutter, disperser, filter) requires the JWST WCS pipeline.

Two viable routes:

1. **`jwst.assign_wcs.nirspec`** — given an MSA shutter (q, s, d), grating, and filter, produce a `gwcs` object that maps shutter coordinates → detector (x_det, y_det) → λ. Walk along the spectral trace and find where it leaves NRS1 (start of gap), enters NRS2, and leaves NRS2. This requires the `jwst` calibration package and the relevant CRDS reference files (camera, collimator, fpa, msa, disperser, filteroffset, wavelengthrange).
2. **Precomputed lookup tables** — for each (disperser, filter) pair, run route 1 once over a (q, s, d) grid (or its V2/V3 equivalent) at app-build time and cache as a NetCDF/Parquet file: columns `q, s, d, disp, filt, lam_blue, lam_gap_lo, lam_gap_hi, lam_red`. At runtime, just look up.

Recommend **route 2** for the app (instant, no CRDS dependency at runtime). Route 1 only as a one-shot offline build step.

Dispersers/filters to support initially: G140M/F070LP, G140M/F100LP, G235M/F170LP, G395M/F290LP, G140H/F070LP, G140H/F100LP, G235H/F170LP, G395H/F290LP, PRISM/CLEAR.

---

## Environment

- Activate with `conda activate stenv` (path: `/Users/sunfengwu/anaconda3`). `stenv` has `astropy`, `pysiaf`, `jwst`, `numpy`, `matplotlib`. **Do not** rely on system `python3` — it lacks numpy.
- `pysiaf` version note: stenv reports "PRDOPSSOC-068 doesn't match online PRDOPSSOC-072". The 138.5° MSA tilt and V2/V3 reference values are stable across these PRDs; safe to ignore unless precision <0.05″ matters. Upgrade only if a real disagreement shows up.

---

## Files in this directory

- `CONTEXT.md` — this file (what / where / why).
- `PLAN.md` — milestone-level implementation plan.
- (later) `app/` — application code.
- (later) `data/` — copy of `nirspec_msa_v2v3.npz`, MSA operability JSON, wavelength-cutoff lookup tables.
- (later) `notebooks/` — prototypes & one-shot builders.

## External references

- Source of MSA file & all coordinate code: `/Users/sunfengwu/jwst_cycle4/footprint_emerald.ipynb` (cells around In[54], In[151], In[164], In[166], In[125]).
- MSA operability docs: jwst-docs.stsci.edu → NIRSpec → MSA → Operability.
- pysiaf NIRSpec apertures: `NRS_FULL_MSA`, `NRS_FULL_MSA1..4`, `NRS_S200A1_SLIT`, fixed-slit names live in the same SIAF.
