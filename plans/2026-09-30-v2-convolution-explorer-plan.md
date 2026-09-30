# v2: convolution-based meltwater anomalies (app first, reusable pipeline)

## Context

v1 of the explorer showed discharge on a globe with hand-drawn region outlines. The site anomalies in the manuscript used an *equilibrium* approximation (the 500-yr weights in `data/meltmodel_site_weights.csv`). v2 uses the real forward model:

- **Method.** Take the source-region δ¹⁸O forcing, convolve it with the HadCM3 pulse (impulse-response) fields, and extract the result at the proxy sites. This follows the old convolution notebook (`old script/[2025.05]_Convolution.ipynb`).
- **Site extraction.** Ocean sites use a ±2° box around the site. Land sites are weighted by the moisture-uptake field from the trajectories, as in the meltmodel repo.
- **Code layout.** Follow the meltmodel repo (`myconfig/`, `mymodules/`, `scripts/`). The pipeline must also feed other Quarto reports and handle more forcings than GLAC-1D and ICE-6G.

**Blocker: the kernels are not on this Mac.** The notebook reads `regions.pkl`, `NA_pulses_dict.pkl` and `dyepulse_series_*.nc` from `/nfs/mary/Users/eelse/...`, and none of them is here. `NA_pulses_dict.pkl` would not be enough anyway, because it only covers 80°W–0°, 20–64°N. That box misses PS2644-5, MD95-2010, NGRIP and most of the Cave Without a Name footprint. Laura will supply the pulse experiment IDs (cold, zonal, merid), and step 1 exports a compact full-grid bundle on the Leeds server.

Decisions made:
- **Sites:** the 9 sites of Fig. 7 (6 marine cores; NISA, Cave Without a Name and NGRIP as land sites), plus a 10th option "Discharge" that shows the forcing itself. Sites are driven by config.
- **Pathways:** `cold`, `zonal` and `mixed`. The mixed schedule comes from notebook cell 73: zonal 21.5–18.0 ka, cold 18.0–14.7, merid 14.7–13.7, zonal 13.7–12.6, cold 12.6–11.3, zonal 11.3–10.0. Each segment is convolved with its own mode's kernel, and the tails are summed across segment boundaries (cell 74). So the merid kernel is needed too.

## Step 1: export script for Leeds (write, test locally on fake data, hand to Laura)

`scripts/export_pulse_kernels.py` has no dependency on this repo beyond xarray, numpy and netCDF4.

- **Config block at the top.** It is pre-filled from `old script/Create_impulsets.ipynb`:
  - pulse runs xpran (cold), xprao (zonal) and xpujc (merid);
  - each pulse run starts after its 10-yr pulse, so the first 10 yr of the parent constant run (xpraj, xprak, xpram) are prepended;
  - the merid parent is xpram (4511–4520), confirmed from the time axes on 2026-09-30;
  - `BASE_DIR = /nfs/annie/earpal/database/experiments`;
  - `PULSE_YEARS = 10`.
  - A time-axis continuity check warns if a pulse run does not start the year after its prepended years.
- **Reading:** open `{exp}/time_series/{exp}.dye0?.annual.nc` with the same pattern as `open_experiment()` in the meltmodel `mymodules/dyefield_computation.py`, renamed to dye00–08.
- **NISA in the notebook** used Miguel's `moisturesource.nc` recharge weights. v2 uses the trajectory uptake masks instead, so the NISA cross-check will differ somewhat.
- **Processing:**
  - Keep the surface level only.
  - Keep latitudes ≥ 0°N.
  - Keep 500 years from the pulse start (10 parent years plus 490 pulse-run years), averaged to decadal means (50 steps).
  - Store as float32 with zlib compression.
- **Output:** `pulse_kernels_surface.nc`, with dimensions (mode, dye, t_decade, lat, lon) and provenance attributes. Estimated size is about 110 MB uncompressed and about 30–50 MB compressed.
- **Optional:** if `NA_pulses_dict.pkl` and `regions.pkl` are found, also write their NISA, PS and NA point kernels to `notebook_point_kernels.nc`. They serve only as a cross-check.
- **At the end:** print file sizes and one `tar czf` command to scp home. The bundle goes to `data/kernels/` here. It is gitignored, and a committed `data/kernels/README.md` records its provenance.
- **Before handing over:** run it here against a small synthetic NetCDF made in the scratchpad, to prove it runs end to end.

## Step 2: reusable pipeline (in this repo, meltmodel style)

```
myconfig/  SITES.py      9 sites (copied from meltmodel PROXYSITES) + domain + uptake key
           DYES.py       dye00–08 ↔ code (Med…GulofMex) ↔ paper name ↔ colour
           FORCINGS.py   registry: name → file, time column, sign, column template, unit
           PATHWAYS.py   cold / zonal / mixed schedule (list of (start, end, mode))
           PATHS.py      meltmodel repo path, kernel bundle path, output dir
mymodules/ forcing.py    load any registered forcing → 10-yr grid, yr BP, 9 region columns
           kernels.py    load bundle; point kernel at a site = site extraction of the field
           sites.py      extract_ocean_proxy / extract_land_proxy / cell_area_2d
                         (ported from meltmodel scripts/precompute_proxymag.py)
           convolution.py convolve(forcing, kernel, pathway) for 1-D series and 2-D fields
scripts/   make_region_geojson.py   HadCM3 input cells (dye_regions_norm.nc) → merged polygons
           make_uptake_weights.py   land_uptakemasks.pkl → data/meltmodel/uptake_weights.nc
           run_convolution.py       CLI: --forcing all|NAME --frame-step 10
R/         site_metrics.R           five metrics + spider plot, sourced by any .qmd (written in step 4)
```

**Forcings.** Adding a forcing means one entry in `FORCINGS.py`. As a fallback, any CSV with a `time_bp` column plus the 9 region-code columns also works. GLAC-1D and ICE-6G are configured as the current files (negate `time` or `t_adj`, `<code> d18O (-35.0)`).

**Convolution maths.** `np.convolve(input_decadal, kernel_decadal, "full")`. The notebook's ×10 was dropped after the scaling check on 2026-09-30 (decision: Laura). The output is extended by 500 yr for the tails.

**Speed.** Convolve the point kernels for the site series. For the maps, use `scipy.signal.fftconvolve` along time on all cells at once, instead of the notebook's per-cell loop that took 30 minutes.

**Outputs** (small, committed, in `outputs/convolution/`):
- `site_anomaly.csv`, tidy: `forcing, pathway, site, region, time_bp, anomaly`.
- `field_decadal_<forcing>_<pathway>.nc` (not committed if too large), which holds the total surface anomaly at every decade, cropped to the North Atlantic and Arctic. FFT convolution along time per cell makes this fast.
- `app/frames/<forcing>_<pathway>.bin` plus `frames_meta.json`, the 8-bit frames for the app.
- `run_info.json`, which records the kernel bundle hash, the settings and the git commit.

**Scaling check (validation).** Convolve a constant input over 500 yr and compare the site value with the equilibrium weight in `data/meltmodel_site_weights.csv` and the meltmodel `mean_dye_{cold,zonal}.nc`. The ratio should be about 1 per site, region and mode. If it isn't, the `*10` convention or `PULSE_AMPLITUDE` is wrong, and Laura has to confirm which. A second check against `notebook_point_kernels.nc` covers NISA and PS, if they were exported.

**Before the data arrives.** The pipeline can use a development-only placeholder kernel, the equilibrium field × (1 − e^(−t/τ)) differenced. It is flagged in `run_info.json` and in the app banner, and never published. This lets the app be built while waiting for the bundle.

## Step 3: explorer v2 (`app/`)

- **`build_app.py`.** It stays standard-library only and reads the pipeline outputs. The fields are exported once by `run_convolution.py` as a compact JSON of rounded values on ocean cells, so the app build needs no xarray. It also takes `--publish-out` as before.
- **Left: globe.**
  - Draw the real HadCM3 input regions from `app/regions_hadcm3.geojson`, using the same `code` property. This replaces `regions_approx.geojson`.
  - Overlay the surface anomaly field as canvas cells, **quasi-continuous at the decadal model resolution** (about 1150 frames, 21.5–10 ka), with play, speed and scrub controls.
  - Frames are quantised to 8 bit (one scale per forcing and pathway) and stored as one binary file per forcing and pathway, about 6 MB each at 1.25° on North Atlantic and Arctic ocean cells. They are published through the artifact's `files` and loaded when that combination is chosen. A preflight keeps each file under 16 MB and the page small.
  - Toggles for forcing (GLAC-1D or ICE-6G) and pathway (cold, zonal or mixed).
  - A diverging or sequential colour scale checked with the dataviz skill in both themes.
  - Site markers are clickable and select the site on the right.
- **Right: time series.**
  - A site selector with the 9 sites and "Discharge". The Glas and speleothem data are removed.
  - A stacked regional contribution chart for the selected forcing and pathway, keeping the existing region colours and stacking order.
  - Thin total lines for the other pathways, for comparison.
  - The globe's time is shown as a cursor.
  - With "Discharge" selected, the view is the v1 stacked discharge.
- **Keep:** the colour tokens in `template.html` stay in sync with `reg_col`. Update `app/README.md` and CLAUDE.md (layout, commands, the new pipeline, and the kernel bundle being untracked).
- **Size budget.** About 1150 decadal frames × about 5k cells, at 1 byte each, is about 6 MB per forcing and pathway, or about 35 MB for all 6 files. `--frame-step` (default 10 yr) and `--coarsen` are the fallbacks if that is too heavy.

## Step 4 (after the app works): spider figures + metrics explained

- **New report:** `site_anomalies_v2.qmd`, with the YAML of `glas_age_range.qmd`. It sources `R/site_metrics.R` and reads `outputs/convolution/site_anomaly.csv`.
- **Figures:**
  - One spider figure with 9 small multiples, one per site. Each compares GLAC-1D and ICE-6G for the mixed pathway.
  - A supplementary figure for cold and zonal.
  - The metrics use the common 500-yr grid, the same rule as v1.
- **"What the five numbers mean" section.** A plain-language explanation of each metric and how it is computed:
  - **Share:** the region's part of the summed anomaly.
  - **Dominant:** how often the region is the largest contributor.
  - **SD:** how much the region's contribution swings over time.
  - **Var. share:** Cov(region, total) / Var(total), meaning how much of the ups and downs of the total the region causes. It sums to 100 %.
  - **SD Δ:** the typical size of a 500-yr jump, which measures abruptness.
  - The definitions come from `nisa_meltwater_sources.qmd` lines 186–200 and 301–305.
- **Numbers:** every number is inline R, and the text is read back with pandoc after the render.

## Verification

1. The export script runs on a synthetic dataset. Laura then runs it on Leeds, and the bundle's dimensions and attributes are checked on arrival.
2. The scaling check passes within a few percent for cold and zonal. If not, stop and ask. Merid has no paper weights. Compare it with `mean_dye_merid.nc` and report the ratio, but a mismatch may be real, because its parent xpram is not the paper's merid run xpral.
3. The mixed pathway is continuous at the boundaries, with no step caused by dropped tails. MWP-1A timing in the Discharge view still peaks at 0.260 Sv at 14.4 ka and 0.243 Sv at 14.0 ka.
4. The app is served with `python3 -m http.server 8765 --bind 127.0.0.1` from `app/` and tested in Chrome:
   - all toggles work and play animates the slices;
   - site clicks sync with the selector;
   - no console errors;
   - no horizontal scroll at 390 px, in light and dark themes.
   Then republish to https://claude.ai/artifact/HMQMC2ACsbpvvDQNyRgcPG.
5. Render the Quarto report and read it back with `pandoc … -t plain`, checking every inline value and qualitative claim.
6. Work on `dev` and commit per step with `ghe-skills:commit`. Don't push unless asked.
