#!/usr/bin/env python3
"""
Export the HadCM3 dye-pulse (impulse-response) fields as one small bundle.

Run this ON THE LEEDS SERVER, where the raw experiment output lives. It needs
Python >= 3.8 with numpy, xarray, dask and netCDF4, and nothing else from this
repository. The system `python` there may be Python 2; use the conda env the
notebooks ran in, e.g.

    /nfs/annie/eelse/conda/envs/py3/bin/python export_pulse_kernels.py --dry-run

What it writes (to OUTDIR)
--------------------------
pulse_kernels_surface.nc
    kernel(mode, dye, lag, latitude, longitude), float32, zlib-compressed.
    Surface level only, latitudes >= LAT_MIN, decadal means of the first
    N_YEARS years after the pulse starts (lag 0 = first decade). The pulse
    years are the first PULSE_YEARS years of the parent constant-input run,
    followed by the pulse run (as in Create_impulsets.ipynb).
notebook_point_kernels.nc   (only if the notebook pickles are found)
    point(site, mode, dye, lag): the NA / NISA / PS pulse responses of
    regions.pkl, and const(site, mode, dye, lag): the 500-yr constant-input
    responses of chp5_constdye.pkl. Used only to cross-check the scaling.

Usage
-----
1. Check PULSE_EXPERIMENTS below (the merid parent needs confirming).
2. python export_pulse_kernels.py --dry-run  # files + time-axis check only
3. python export_pulse_kernels.py
4. Copy the printed tar file home and unpack it into data/kernels/.

Author: Laura Endres
"""

import argparse
import datetime
import os
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

# Model years (5000s) are outside the numpy datetime range; cftime is fine.
warnings.filterwarnings("ignore", category=xr.SerializationWarning)


# =============================================================================
# Configuration: edit here
# =============================================================================

# Pulse (impulse-response) experiments per AMOC mode, as in the notebook
# Create_impulsets.ipynb. Each pulse run starts AFTER its 10-yr dye pulse;
# the pulse years themselves are the first PULSE_YEARS years of the parent
# constant-input run, and are prepended (notebook cell 2).
#   cold : xpran after xpraj 5681-5690
#   zonal: xprao after xprak 5101-5110
#   merid: xpujc after xpram 4511-4520. The notebook names xpram but never
#          loads it (and slices the wrong years), so its merid pulse was
#          skipped. Confirmed from the time axes: xpram runs 4511-5010,
#          xpujc starts 4521. xpram is not the paper's merid run (xpral).
PULSE_EXPERIMENTS = {
    "cold": {"pulse": "xpran", "parent": "xpraj"},
    "zonal": {"pulse": "xprao", "parent": "xprak"},
    "merid": {"pulse": "xpujc", "parent": "xpram"},
}

# Same layout as mymodules/dyefield_computation.py in the meltmodel repo:
# BASE_DIR/<exp>/time_series/<exp>.dyeNN.annual.nc, one file per dye.
# Some runs carry a 10th tracer (dye09) that is not one of the nine input
# regions; only dye00-dye08 are read (the notebook also left dye09 out).
BASE_DIR = Path("/nfs/annie/earpal/database/experiments")
FILE_PATTERN = "{exp}/time_series/{exp}.{dye}.annual.nc"

# Pulse design: 10 yr of the constant-input run's dye flux, then zero.
PULSE_YEARS = 10
PULSE_AMPLITUDE = "same dye flux as the parent constant-input run"

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
    """The nine region-dye files of an experiment (missing ones left out)."""
    files = [BASE_DIR / FILE_PATTERN.format(exp=exp, dye=d) for d in DYES]
    return [f for f in files if f.exists()]


def extra_files(exp):
    """Dye files present but not used (e.g. dye09)."""
    used = set(experiment_files(exp))
    pattern = FILE_PATTERN.format(exp=exp, dye="dye*")
    return sorted(f for f in BASE_DIR.glob(pattern) if f not in used)


def open_experiment(exp):
    """Open dye00...dye08 of one experiment, one file per dye.

    Each file is named after its dye, so the variable inside is renamed by
    file, not by position in a merged dataset.
    """
    files = experiment_files(exp)
    if len(files) != len(DYES):
        raise FileNotFoundError(
            f"{exp}: expected {len(DYES)} dye files (dye00-dye08) under {BASE_DIR / exp}, "
            f"found {len(files)}")
    parts = {}
    for dye, f in zip(DYES, files):
        ds = xr.open_dataset(f, chunks={})
        vars4d = [v for v in ds.data_vars if ds[v].ndim >= 3]
        if len(vars4d) != 1:
            raise ValueError(f"{f}: expected one gridded variable, found {vars4d}")
        parts[dye] = ds[vars4d[0]]
    return xr.Dataset(parts)


def dim_like(da, *names):
    """Return the first dimension of ``da`` whose name starts with one of ``names``."""
    for d in da.dims:
        if any(d.lower().startswith(n) for n in names):
            return d
    return None


def surface(da, n_years):
    """First `n_years` of the surface level at lat >= LAT_MIN, as values.

    Time, depth and latitude are cut by index before anything is read, so
    only the needed part of each file comes over NFS.
    """
    tdim = dim_like(da, "t")
    zdim = dim_like(da, "depth", "lev", "z")
    ydim = dim_like(da, "lat")
    xdim = dim_like(da, "lon")
    da = da.isel({tdim: slice(0, n_years)})
    if zdim is not None:
        da = da.isel({zdim: 0}, drop=True)
    da = da.isel({ydim: np.asarray(da[ydim].values >= LAT_MIN)})
    da = da.rename({ydim: "latitude", xdim: "longitude"}).transpose(tdim, "latitude", "longitude")
    return da.values.astype("float32"), da["latitude"].values, da["longitude"].values


def decadal_means(vals):
    """Average the first N_YEARS annual steps into DECADE-yr means."""
    if vals.shape[0] < N_YEARS:
        raise ValueError(f"only {vals.shape[0]} years after the pulse start, need {N_YEARS}")
    n_dec = N_YEARS // DECADE
    return vals[: n_dec * DECADE].reshape(n_dec, DECADE, *vals.shape[1:]).mean(axis=1)


def time_years(exp):
    """Model years of an experiment, read from the time variable of dye00 only.

    Uses netCDF4 directly: no data and no other coordinates are read, which
    keeps the check fast on NFS.
    """
    import netCDF4
    f = BASE_DIR / FILE_PATTERN.format(exp=exp, dye=DYES[0])
    with netCDF4.Dataset(f) as nc:
        grid = [v for v in nc.variables.values() if v.ndim >= 3][0]
        tname = [d for d in grid.dimensions if d.lower().startswith("t")][0]
        tv = nc.variables[tname]
        dates = netCDF4.num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))
    return [d.year for d in dates]


def check_continuity(mode, exps):
    """The pulse run should start the year after the prepended parent years."""
    tp, tq = time_years(exps["parent"]), time_years(exps["pulse"])
    first, last, nxt = tp[0], tp[PULSE_YEARS - 1], tq[0]
    n_total = PULSE_YEARS + len(tq)
    msg = (f"  {mode:6s} pulse years {first}-{last} from {exps['parent']}, "
           f"{exps['pulse']} starts {nxt}, {n_total} yr in total")
    if nxt != last + 1:
        msg += f"  <-- WARNING: expected {exps['pulse']} to start in {last + 1}"
    if n_total < N_YEARS:
        msg += f"  <-- WARNING: fewer than N_YEARS={N_YEARS}"
    print(msg, flush=True)


def export_fields(dry_run=False):
    modes = list(PULSE_EXPERIMENTS)
    for mode, exps in PULSE_EXPERIMENTS.items():
        for role in ("parent", "pulse"):
            files = experiment_files(exps[role])
            print(f"  {mode:6s} {role:6s} {exps[role]}: {len(files)} of {len(DYES)} dye files")
            for f in files:
                print(f"         {f}")
            for f in extra_files(exps[role]):
                print(f"         (not used) {f.name}")
    print("Time axes:", flush=True)
    ok = True
    for mode, exps in PULSE_EXPERIMENTS.items():
        try:
            check_continuity(mode, exps)
        except (OSError, IndexError, KeyError) as err:
            print(f"  {mode:6s} cannot read time axis: {err}")
            ok = False
    if dry_run:
        return None
    if not ok:
        sys.exit("Cannot export: fix the experiments above first.")

    kernel = None
    for m, mode in enumerate(modes):
        parent_ds = open_experiment(PULSE_EXPERIMENTS[mode]["parent"])
        pulse_ds = open_experiment(PULSE_EXPERIMENTS[mode]["pulse"])
        for d, dye in enumerate(DYES):
            head, lat, lon = surface(parent_ds[dye], PULSE_YEARS)
            tail, lat2, lon2 = surface(pulse_ds[dye], N_YEARS - PULSE_YEARS)
            if not (np.array_equal(lat, lat2) and np.array_equal(lon, lon2)):
                raise ValueError(f"{mode} {dye}: parent and pulse grids differ")
            vals = decadal_means(np.concatenate([head, tail], axis=0))
            if kernel is None:
                kernel = np.full((len(modes), len(DYES), *vals.shape), np.nan, dtype="float32")
            kernel[m, d] = vals
            print(f"  {mode} {dye}: max {np.nanmax(vals):.4g}", flush=True)
        parent_ds.close()
        pulse_ds.close()

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
        "experiments": ", ".join(f"{m}={e['pulse']} (pulse years from {e['parent']})"
                                 for m, e in PULSE_EXPERIMENTS.items()),
        "pulse_years": str(PULSE_YEARS),
        "pulse_amplitude": str(PULSE_AMPLITUDE),
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
