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
   | merid | xpujc | xpral, **to confirm**: the notebook names `xpram`, which it never loads |

   Each pulse run starts *after* its 10-yr pulse, so the script prepends the
   first 10 years of the parent constant-input run, as the notebook did.
3. Run `python export_pulse_kernels.py --dry-run`. It should find 9 dye files
   per experiment. It also prints each time axis and warns if a pulse run does
   not start the year after its prepended pulse years.
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
