#!/usr/bin/env python3
"""
Export the two remaining inputs of the forward-model paper, on the Leeds server.

1. constant_runs_surface.nc
   field(run, dye, lag, latitude, longitude): surface dye of the 500-yr
   constant-input runs, as decadal means, north of LAT_MIN. Used for the
   spatial validation (RMSE maps): the convolution of a constant forcing with
   the pulse kernels against these runs.
2. amoc_poeppelmeier_25N.csv
   AMOC strength (Sv) at 25 N, averaged over 1000-1500 m, from the deglacial
   best-fit simulation of Poeppelmeier et al. (2023), as used in the 2025.05
   convolution notebook.

Copy this file next to export_pulse_kernels.py (it reuses its readers) and run
it with Python >= 3.8 and xarray, e.g.

    /nfs/annie/eelse/conda/envs/py3/bin/python export_validation_inputs.py --dry-run
    /nfs/annie/eelse/conda/envs/py3/bin/python export_validation_inputs.py

Author: Laura Endres
"""

import argparse
import datetime
import os
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import export_pulse_kernels as epk  # noqa: E402  (same folder; shared readers and settings)

# =============================================================================
# Configuration: edit here
# =============================================================================

# Constant-input runs (17.8 ka). xpral is the paper's meridional run; xpram is
# the parent of the meridional pulse run, kept to check that kernel too.
CONSTANT_EXPERIMENTS = {"cold": "xpraj", "zonal": "xprak", "merid": "xpral", "merid_xpram": "xpram"}

AMOC_FILE = Path("/nfs/mary/Users/eelse/work/Poeppelmeier2022_DeglacialBestFit.nc")
AMOC_VAR = "AtlanticStreamfunction"
AMOC_LAT = 25.0
AMOC_DEPTH = (1000.0, 1500.0)

OUTDIR = Path("./validation_export")


def export_constant(dry_run=False):
    runs = list(CONSTANT_EXPERIMENTS)
    for run, exp in CONSTANT_EXPERIMENTS.items():
        files = epk.experiment_files(exp)
        years = epk.time_years(exp) if files else []
        print(f"  {run:12s} {exp}: {len(files)} of {len(epk.DYES)} dye files, "
              f"{len(years)} yr" + (f" ({years[0]}-{years[-1]})" if years else ""), flush=True)
    if dry_run:
        return None
    field = None
    for r, run in enumerate(runs):
        ds = epk.open_experiment(CONSTANT_EXPERIMENTS[run])
        for d, dye in enumerate(epk.DYES):
            vals, lat, lon = epk.surface(ds[dye], epk.N_YEARS)
            vals = epk.decadal_means(vals)
            if field is None:
                field = np.full((len(runs), len(epk.DYES), *vals.shape), np.nan, dtype="float32")
            field[r, d] = vals
        ds.close()
        print(f"  {run} done", flush=True)
    n_dec = field.shape[2]
    out = xr.Dataset(
        {"field": (("run", "dye", "lag", "latitude", "longitude"), field)},
        coords={"run": runs, "experiment": ("run", [CONSTANT_EXPERIMENTS[r] for r in runs]),
                "dye": epk.DYES, "region_code": ("dye", epk.REGION_CODES),
                "lag": ("lag", np.arange(n_dec) * epk.DECADE), "latitude": lat, "longitude": lon},
    )
    out["lag"].attrs = {"units": "years", "long_name": "start of the 10-yr averaging window after the start of input"}
    out["field"].attrs = {"long_name": "surface dye under constant input", "note": "10-yr means of annual output"}
    out.attrs = {"title": "HadCM3 constant-input dye runs (surface), for validating the pulse kernels",
                 "experiments": ", ".join(f"{r}={e}" for r, e in CONSTANT_EXPERIMENTS.items()),
                 "n_years": epk.N_YEARS, "lat_min": epk.LAT_MIN,
                 "created": datetime.datetime.now().isoformat(timespec="seconds"),
                 "created_by": "export_validation_inputs.py"}
    path = OUTDIR / "constant_runs_surface.nc"
    out.to_netcdf(path, encoding={"field": {"zlib": True, "complevel": 4,
                                            "chunksizes": (1, 1, n_dec, len(lat), len(lon))}})
    return path


def export_amoc(dry_run=False):
    print(f"  {AMOC_FILE}: {'found' if AMOC_FILE.exists() else 'NOT FOUND'}", flush=True)
    if dry_run or not AMOC_FILE.exists():
        return None
    ds = xr.open_dataset(AMOC_FILE)
    sf = ds[AMOC_VAR].sel(lat_u=AMOC_LAT, method="nearest").sel(z_w=slice(*AMOC_DEPTH)).mean("z_w")
    tname = [d for d in sf.dims][0]
    out = sf.to_dataframe(name="amoc_sv").reset_index()[[tname, "amoc_sv"]]
    path = OUTDIR / "amoc_poeppelmeier_25N.csv"
    with open(path, "w") as f:
        f.write(f"# AMOC strength: {AMOC_VAR} at lat_u={float(sf.lat_u):g}, mean over z_w {AMOC_DEPTH[0]:g}-{AMOC_DEPTH[1]:g} m\n")
        f.write(f"# source: {AMOC_FILE} (Poeppelmeier et al. 2023); time units: {ds[tname].attrs.get('units', 'see source')}\n")
        out.to_csv(f, index=False)
    return path


def main():
    global OUTDIR, AMOC_FILE
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="only check the inputs")
    ap.add_argument("--base-dir", type=Path, help=f"override BASE_DIR ({epk.BASE_DIR})")
    ap.add_argument("--amoc-file", type=Path, help=f"override AMOC_FILE ({AMOC_FILE})")
    ap.add_argument("--outdir", type=Path, help=f"override OUTDIR ({OUTDIR})")
    args = ap.parse_args()
    if args.base_dir:
        epk.BASE_DIR = args.base_dir
    if args.amoc_file:
        AMOC_FILE = args.amoc_file
    if args.outdir:
        OUTDIR = args.outdir
    if not args.dry_run:
        OUTDIR.mkdir(parents=True, exist_ok=True)
    print("Constant-input runs:")
    written = [export_constant(args.dry_run)]
    print("AMOC series:")
    written.append(export_amoc(args.dry_run))
    written = [p for p in written if p is not None]
    if args.dry_run:
        return
    print("\nWritten:")
    for p in written:
        print(f"  {p}  ({os.path.getsize(p) / 1e6:.1f} MB)")
    tar = OUTDIR.resolve().parent / "validation_inputs.tar.gz"
    print(f"\nPack and copy home:\n  tar czf {tar} -C {OUTDIR.resolve()} {' '.join(p.name for p in written)}")


if __name__ == "__main__":
    main()
