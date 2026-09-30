#!/usr/bin/env python3
"""
Moisture-uptake weights for the land sites, from the trajectory files.

Reads the meltmodel repo's precipitation-weighted uptake fields
(data/trajectories/<loc>_<run>_th00_UTOT_weighted.nc, 0.5 deg) and regrids
them to the HadCM3 1.25 deg dye grid exactly as the meltmodel
scripts/fig6_trajectories.py does (lon wrapped to 0-360, duplicates dropped,
linear interpolation). Writes data/meltmodel/uptake_weights.nc with
uptake(site, run, latitude, longitude).

Run from the repository root:
    python3 scripts/make_uptake_weights.py
"""

import pickle
import sys
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from myconfig.PATHS import MEAN_DYE, MELTMODEL_DATA, TRAJECTORIES, UPTAKE_WEIGHTS
from myconfig.SITES import SITES, UPTAKE_RUN

RUNS = sorted(set(UPTAKE_RUN.values()))


def regrid(path, lat_new, lon_new):
    da = xr.open_dataarray(path).squeeze()
    da = da.assign_coords(lon=("dimx_N", np.linspace(-180, 180, 721)),
                          lat=("dimy_N", np.linspace(-90, 90, 361)))
    da = da.swap_dims({"dimx_N": "lon", "dimy_N": "lat"})
    da = da.assign_coords(lon=da.lon % 360).sortby("lon")
    _, idx = np.unique(da.lon, return_index=True)
    da = da.isel(lon=np.sort(idx))
    return da.interp(lat=lat_new, lon=lon_new, method="linear")


def main():
    grid = xr.open_dataset(str(MEAN_DYE).format(mode="zonal"))
    locs = sorted({s["uptake"] for s in SITES.values() if s["domain"] == "land"})
    arr = []
    for loc in locs:
        arr.append(xr.concat(
            [regrid(str(TRAJECTORIES).format(loc=loc, exp=run), grid.lat, grid.lon) for run in RUNS],
            dim="run"))
    out = xr.concat(arr, dim="site").assign_coords(site=locs, run=RUNS)
    out = out.drop_vars([c for c in out.coords if c not in ("site", "run", "lat", "lon")])
    out = out.rename({"lat": "latitude", "lon": "longitude"}).astype("float32")
    out.name = "uptake"
    out.attrs = {
        "long_name": "precipitation-weighted moisture uptake, regridded to the HadCM3 1.25 deg grid",
        "units": "% / 1e5 km2",
        "source": str(TRAJECTORIES),
        "method": "as meltmodel scripts/fig6_trajectories.py (threshold th00, linear interpolation)",
    }
    UPTAKE_WEIGHTS.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(UPTAKE_WEIGHTS, encoding={"uptake": {"zlib": True, "complevel": 4}})
    print(f"wrote {UPTAKE_WEIGHTS}")

    # Cross-check against the meltmodel cache made by fig6_trajectories.py.
    pkl = MELTMODEL_DATA / "intermediates/dyestuff_modelpaper/land_uptakemasks.pkl"
    if pkl.exists():
        ref = pickle.load(open(pkl, "rb"))
        for loc in locs:
            for run in RUNS:
                a = out.sel(site=loc, run=run).values
                b = ref[loc][run].values
                d = np.nanmax(np.abs(a - b)) / np.nanmax(np.abs(b))
                print(f"  {loc:15s} {run}: max rel. difference to land_uptakemasks.pkl {d:.1e}")


if __name__ == "__main__":
    main()
