# Pulse (impulse-response) kernels

This folder holds the HadCM3 dye-pulse fields that the convolution pipeline
needs. They are **not in git**, because they are too large. Create them on the
Leeds server with `scripts/export_pulse_kernels.py`, then copy them here.

## How to make the bundle

1. Copy `scripts/export_pulse_kernels.py` to the Leeds server.
2. Check `PULSE_EXPERIMENTS` in its configuration block. It is pre-filled
   from `Create_impulsets.ipynb`:

   | Mode | Pulse run | Pulse years from (parent) |
   |---|---|---|
   | cold | xpran | xpraj, years 5681–5690 |
   | zonal | xprao | xprak, years 5101–5110 |
   | merid | xpujc | xpram, years 4511–4520 (checked from the time axes; the notebook named xpram but never loaded it) |

   Each pulse run starts *after* its 10-yr pulse, so the script prepends the
   first 10 years of the parent constant-input run, as the notebook did.
3. Run `python export_pulse_kernels.py --dry-run` with **Python 3.8 or newer**.
   The system `python` on the server may be Python 2, which fails with
   `SyntaxError` at the first f-string. Use the notebooks' environment
   instead, `/nfs/annie/eelse/conda/envs/py3/bin/python`, or run
   `conda activate py3` first. It should find 9 of 9 dye files
   (dye00–dye08) per experiment. Some runs also carry a 10th tracer, `dye09`,
   which is not one of the nine input regions; it is listed as "not used". It also prints each time axis and warns if a pulse run does
   not start the year after its prepended pulse years.
4. Run `python export_pulse_kernels.py` in the same environment. It needs
   numpy, xarray, dask and netCDF4.
5. Copy the `pulse_kernels.tar.gz` it prints home and unpack it in this folder.

## Files

| File | Content |
|---|---|
| `pulse_kernels_surface.nc` | `kernel(mode, dye, lag, latitude, longitude)`: surface dye after a pulse, as 10-yr means of the first 500 yr, north of 0°N |
| `notebook_point_kernels.nc` | Optional. The NA, NISA and PS series from the 2025.05 convolution notebook's pickles, used only to cross-check the scaling |

The file attributes record the experiment IDs, the pulse design and the
creation date.

## Caveat: the merid pulse has a different parent run

The cold and zonal pulse runs branch from the paper's constant-input runs
for 17.8 ka (xpraj and xprak, as listed in the meltmodel `EXPERIMENTS.py`).
The merid pulse run xpujc does **not**. It branches from **xpram**, which
starts in 4511, while xpujc starts in 4521. The paper's merid 17.8 ka run is
**xpral**, which starts in 6381. xpram does not appear in the meltmodel
`EXPERIMENTS.py`. This was established from the time axes on the Leeds
server on 2026-09-30.

Consequences:

- The merid kernel's background climate is xpram's, not that of the merid
  equilibrium field used in Endres et al. (2026b). Before relying on the
  merid kernel, check which boundary conditions and AMOC state xpram has.
- The scaling check cannot test merid against the paper's weights:
  `data/meltmodel_site_weights.csv` has only the cold and zonal modes. The
  nearest reference is the meltmodel `mean_dye_merid.nc`, which comes from
  the paper's merid runs, not from xpram. A mismatch there may be real.
- The merid kernel is used only in the **mixed** pathway, for 14.7–13.7 ka.
  The cold and zonal pathways do not depend on it.

## Checks on the bundle (2026-09-30)

- **Contents.** `kernel(mode, dye, lag, latitude, longitude)` with 3 modes, 9 dyes, 50 decadal lags and 72 × 288 cells (0–90°N). About 48 % of cells are land (NaN). Small negative values, down to −0.037, hold 0.1 % of the dye mass and are left as they are.
- **Pulse continuity.** In the source regions, the second decade (the first of the pulse run) is 0.2–0.35× the first decade (the parent's pulse years), and the dye keeps decaying after that. So the pulse run continues a pulse of 1 unit per year, not 10.
- **Scaling.** The 50 decadal responses add up to the constant-input equilibrium:
  - Cold and zonal, per site (weighted over regions): 0.85–1.16 against `data/meltmodel_site_weights.csv`.
  - Over the North Atlantic box: 1.00 (cold) and 1.02 (zonal).
  - So with no extra factor, a steady forcing reaches the paper's equilibrium weights. The 2025.05 notebook's deglacial cells multiplied the forcing by 10. The pipeline runs at ×1, and the explorer can switch to ×10. **Which factor is right is under review** (see "Open questions" in CLAUDE.md).
- **Merid.** Against the paper's `mean_dye_merid.nc`, the ratio is 1.08 over the North Atlantic box but 0.49 (PS2644-5) to 1.16 at single sites. This fits the xpram parent (see the caveat above).
- **Notebook cross-check** (`notebook_point_kernels.nc`):
  - PS2644-5, cold and zonal: identical to the notebook (r = 1.000, amplitude 1.00).
  - NISA: same shape (r ≈ 0.998), but 1.4× (cold) and 1.7× (zonal) larger. The notebook weighted NISA with `moisturesource.nc`, while v2 uses the trajectory uptake masks.
