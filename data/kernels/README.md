# Pulse (impulse-response) kernels

This folder holds the HadCM3 dye-pulse fields that the convolution pipeline
needs. They are **not in git**, because they are too large. Create them on the
Leeds server with `scripts/export_pulse_kernels.py`, then copy them here.

## How to make the bundle

1. Copy `scripts/export_pulse_kernels.py` to the Leeds server.
2. In its configuration block, fill in:
   - `PULSE_EXPERIMENTS`: the pulse experiment ID for each of the `cold`,
     `zonal` and `merid` AMOC modes;
   - `PULSE_YEARS` and `PULSE_AMPLITUDE`: the pulse design.
3. Check the input files with `python export_pulse_kernels.py --dry-run`. It
   should find 9 dye files per experiment.
4. Run `python export_pulse_kernels.py`. It needs numpy, xarray and netCDF4;
   the `dyetracer` environment has them.
5. Copy the `pulse_kernels.tar.gz` it prints home and unpack it in this folder.

## Files

| File | Content |
|---|---|
| `pulse_kernels_surface.nc` | `kernel(mode, dye, lag, latitude, longitude)`: surface dye after a pulse, as 10-yr means of the first 500 yr, north of 0°N |
| `notebook_point_kernels.nc` | Optional. The NA, NISA and PS series from the 2025.05 convolution notebook's pickles, used only to cross-check the scaling |

The file attributes record the experiment IDs, the pulse design and the
creation date.
