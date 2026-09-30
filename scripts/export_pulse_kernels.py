#!/usr/bin/env python
"""
Export the HadCM3 dye-pulse (impulse-response) fields as one small bundle.

Run this ON THE LEEDS SERVER, where the raw experiment output lives. It needs
only numpy, xarray and netCDF4 (the ``dyetracer`` environment is fine), and
nothing else from this repository.

What it writes (to OUTDIR)
--------------------------
pulse_kernels_surface.nc
    kernel(mode, dye, lag, latitude, longitude), float32, zlib-compressed.
    Surface level only, latitudes >= LAT_MIN, decadal means of the first
    N_YEARS years after the pulse starts (lag 0 = first decade).
notebook_point_kernels.nc   (only if the notebook pickles are found)
    point(site, mode, dye, lag): the NA / NISA / PS pulse responses of
    regions.pkl, and const(site, mode, dye, lag): the 500-yr constant-input
    responses of chp5_constdye.pkl. Used only to cross-check the scaling.

Usage
-----
1. Fill in PULSE_EXPERIMENTS (and check the other settings) below.
2. python export_pulse_kernels.py            # or: --dry-run to list files only
3. Copy the printed tar file home and unpack it into data/kernels/.

Author: Laura Endres
"""

import argparse
import datetime
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import xarray as xr


# =============================================================================
# Configuration: edit here
# =============================================================================

# Pulse (impulse-response) experiment ID per AMOC mode.
PULSE_EXPERIMENTS = {
    "cold": "TODO",
    "zonal": "TODO",
    "merid": "TODO",
}

# Same layout as mymodules/dyefield_computation.py in the meltmodel repo:
# BASE_DIR/<exp>/time_series/<exp>.dye0?.annual.nc
BASE_DIR = Path("/nfs/annie/earpal/database/experiments")
FILE_PATTERN = "{exp}/time_series/{exp}.dye0?.annual.nc"

# Pulse design, written into the file attributes (please fill in).
PULSE_YEARS = None          # length of the dye pulse in years, e.g. 10
PULSE_AMPLITUDE = None      # dye input during the pulse (units as in the model)

# Index (in years from the first file time step) at which the pulse starts.
PULSE_START_INDEX = 0
N_YEARS = 500               # length of the response kept
DECADE = 10                 # averaging window (years)
LAT_MIN = 0.0               # keep latitudes >= LAT_MIN (degrees N)

# Optional notebook pickles for the cross-check (skipped if missing).
NOTEBOOK_DIR = Path("/nfs/mary/Users/eelse/scripts/jupyter")
NOTEBOOK_PICKLES = {
    "point": NOTEBOOK_DIR / "regions.pkl",
    "const": NOTEBOOK_DIR / "myintermediates/dyestuff_modelpaper/chp5_constdye.pkl",
}

OUTDIR = Path("./pulse_kernel_export")

DYES = [f"dye{i:02d}" for i in range(9)]
REGION_CODES = ["Med", "Bri", "Fen", "EurArc", "AmeArc", "GIS", "NLau", "SLau", "GulofMex"]


# =============================================================================
# Helpers
# =============================================================================

def experiment_files(exp):
    return sorted(BASE_DIR.glob(FILE_PATTERN.format(exp=exp)))


def open_experiment(exp):
    """Open one experiment, rename dye variables to dye00...dye08."""
    files = experiment_files(exp)
    if len(files) != 9:
        raise FileNotFoundError(
            f"{exp}: expected 9 dye files under {BASE_DIR / exp}, found {len(files)}")
    ds = xr.open_mfdataset([str(f) for f in files], combine="by_coords", chunks={})
    return ds.rename_vars({old: DYES[i] for i, old in enumerate(ds.data_vars)})


def dim_like(da, *names):
    """Return the first dimension of ``da`` whose name starts with one of ``names``."""
    for d in da.dims:
        if any(d.lower().startswith(n) for n in names):
            return d
    return None


def surface_decadal(da):
    """Surface level, lat >= LAT_MIN, decadal means of the first N_YEARS years."""
    tdim = dim_like(da, "t")
    zdim = dim_like(da, "depth", "lev", "z")
    ydim = dim_like(da, "lat")
    xdim = dim_like(da, "lon")
    if zdim is not None:
        da = da.isel({zdim: 0}, drop=True)
    da = da.isel({tdim: slice(PULSE_START_INDEX, PULSE_START_INDEX + N_YEARS)})
    if da.sizes[tdim] < N_YEARS:
        raise ValueError(f"only {da.sizes[tdim]} years after the pulse start, need {N_YEARS}")
    da = da.where(da[ydim] >= LAT_MIN, drop=True)
    da = da.rename({ydim: "latitude", xdim: "longitude"})
    vals = da.transpose(tdim, "latitude", "longitude").values.astype("float32")
    n_dec = N_YEARS // DECADE
    vals = vals[: n_dec * DECADE].reshape(n_dec, DECADE, *vals.shape[1:]).mean(axis=1)
    return vals, da["latitude"].values, da["longitude"].values


def export_fields(dry_run=False):
    modes = list(PULSE_EXPERIMENTS)
    for mode, exp in PULSE_EXPERIMENTS.items():
        files = experiment_files(exp)
        print(f"  {mode:6s} {exp}: {len(files)} dye files")
        for f in files:
            print(f"         {f}")
    if dry_run:
        return None
    if any(e == "TODO" for e in PULSE_EXPERIMENTS.values()):
        sys.exit("Fill in PULSE_EXPERIMENTS first.")

    kernel = None
    for m, mode in enumerate(modes):
        ds = open_experiment(PULSE_EXPERIMENTS[mode])
        for d, dye in enumerate(DYES):
            vals, lat, lon = surface_decadal(ds[dye])
            if kernel is None:
                kernel = np.full((len(modes), len(DYES), *vals.shape), np.nan, dtype="float32")
            kernel[m, d] = vals
            print(f"  {mode} {dye}: max {np.nanmax(vals):.4g}")
        ds.close()

    n_dec = kernel.shape[2]
    out = xr.Dataset(
        {"kernel": (("mode", "dye", "lag", "latitude", "longitude"), kernel)},
        coords={
            "mode": modes,
            "dye": DYES,
            "region_code": ("dye", REGION_CODES),
            "lag": ("lag", np.arange(n_dec) * DECADE),
            "latitude": lat,
            "longitude": lon,
        },
    )
    out["lag"].attrs = {"units": "years", "long_name": f"start of {DECADE}-yr averaging window after pulse start"}
    out["kernel"].attrs = {"long_name": "surface dye concentration after a dye pulse",
                           "note": f"{DECADE}-yr means of annual output"}
    out.attrs = {
        "title": "HadCM3 dye-pulse impulse-response fields (surface)",
        "experiments": ", ".join(f"{m}={e}" for m, e in PULSE_EXPERIMENTS.items()),
        "pulse_years": str(PULSE_YEARS),
        "pulse_amplitude": str(PULSE_AMPLITUDE),
        "pulse_start_index": PULSE_START_INDEX,
        "n_years": N_YEARS,
        "lat_min": LAT_MIN,
        "source_pattern": str(BASE_DIR / FILE_PATTERN),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "created_by": "export_pulse_kernels.py",
    }
    path = OUTDIR / "pulse_kernels_surface.nc"
    enc = {"kernel": {"zlib": True, "complevel": 4,
                      "chunksizes": (1, 1, n_dec, len(lat), len(lon))}}
    out.to_netcdf(path, encoding=enc)
    return path


def as_series(obj, dye):
    """Pull one dye's 1-D series out of a Dataset / dict, as a float array."""
    da = obj[dye]
    tdim = dim_like(da, "t")
    da = da.squeeze(drop=True)
    return np.asarray(da.transpose(tdim).values, dtype="float32")


def decadal(series):
    n = min(len(series), N_YEARS) // DECADE
    return series[: n * DECADE].reshape(n, DECADE).mean(axis=1)


def export_notebook_points():
    found = {k: p for k, p in NOTEBOOK_PICKLES.items() if p.exists()}
    if not found:
        print("  notebook pickles not found, skipping cross-check file")
        return None
    arrays, sites, modes = {}, None, ["cold", "zonal", "merid"]
    for key, p in found.items():
        with open(p, "rb") as f:
            obj = pickle.load(f)
        sites = sorted(obj) if sites is None else sites
        arr = np.full((len(sites), len(modes), len(DYES), N_YEARS // DECADE), np.nan, dtype="float32")
        for s, site in enumerate(sites):
            for m, mode in enumerate(modes):
                if site not in obj or mode not in obj[site]:
                    continue
                for d, dye in enumerate(DYES):
                    try:
                        v = decadal(as_series(obj[site][mode], dye))
                    except Exception as err:  # keep going; report what failed
                        print(f"  {key} {site} {mode} {dye}: skipped ({err})")
                        continue
                    arr[s, m, d, : len(v)] = v
        arrays[key] = arr
        print(f"  {key}: {p}")
    out = xr.Dataset(
        {k: (("site", "mode", "dye", "lag"), v) for k, v in arrays.items()},
        coords={"site": sites, "mode": modes, "dye": DYES,
                "lag": np.arange(N_YEARS // DECADE) * DECADE},
    )
    out.attrs = {"note": "Notebook (2025.05 convolution) point series, decadal means; cross-check only",
                 **{f"source_{k}": str(p) for k, p in found.items()}}
    path = OUTDIR / "notebook_point_kernels.nc"
    out.to_netcdf(path)
    return path


def main():
    global BASE_DIR, OUTDIR, NOTEBOOK_DIR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="only list the input files")
    ap.add_argument("--base-dir", type=Path, help=f"override BASE_DIR ({BASE_DIR})")
    ap.add_argument("--outdir", type=Path, help=f"override OUTDIR ({OUTDIR})")
    ap.add_argument("--notebook-dir", type=Path, help=f"override NOTEBOOK_DIR ({NOTEBOOK_DIR})")
    args = ap.parse_args()
    if args.base_dir:
        BASE_DIR = args.base_dir
    if args.outdir:
        OUTDIR = args.outdir
    if args.notebook_dir:
        for k, p in NOTEBOOK_PICKLES.items():
            NOTEBOOK_PICKLES[k] = args.notebook_dir / p.relative_to(NOTEBOOK_DIR)

    print("Pulse fields:")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    written = [export_fields(dry_run=args.dry_run)]
    if args.dry_run:
        return
    print("Notebook point series:")
    written.append(export_notebook_points())

    written = [p for p in written if p is not None]
    print("\nWritten:")
    for p in written:
        print(f"  {p}  ({os.path.getsize(p) / 1e6:.1f} MB)")
    tar = OUTDIR.resolve().parent / "pulse_kernels.tar.gz"
    names = " ".join(p.name for p in written)
    print(f"\nPack and copy home:\n  tar czf {tar} -C {OUTDIR.resolve()} {names}")
    print(f"  scp <user>@<host>:{tar} <repo>/data/kernels/   # then: tar xzf pulse_kernels.tar.gz")


if __name__ == "__main__":
    main()
