"""Extract a surface field at a proxy site.

Same maths as the meltmodel repo (scripts/precompute_proxymag.py):
ocean sites are a latitude-weighted mean over a small box, land sites a sum
weighted by moisture uptake times grid-cell area. Works on fields with any
leading dimensions (for example mode, dye, lag), because the reductions only
touch latitude and longitude.
"""

import numpy as np
import xarray as xr

R_EARTH = 6.371e6


def cell_area(lat, dlat=1.25, dlon=1.25):
    """Grid-cell area (m^2) on a regular grid, as a latitude DataArray."""
    lat_rad = np.deg2rad(np.asarray(lat, float))
    a = R_EARTH ** 2 * np.abs(np.sin(lat_rad + np.deg2rad(dlat) / 2)
                              - np.sin(lat_rad - np.deg2rad(dlat) / 2)) * np.deg2rad(dlon)
    return xr.DataArray(a, coords={"latitude": lat}, dims="latitude")


def extract_ocean(field, lat, lon, box=2):
    """Latitude-weighted mean of `field` in a +/- `box` degree box."""
    lon = lon % 360
    sub = field.sel(latitude=slice(lat - box, lat + box),
                    longitude=slice(lon - box, lon + box))
    if sub.sizes["latitude"] == 0 or sub.sizes["longitude"] == 0:
        raise ValueError(f"site ({lat}, {lon}) is outside the field")
    w = np.cos(np.deg2rad(sub.latitude))
    return sub.weighted(w).mean(("latitude", "longitude"), skipna=True)


def land_weights(uptake):
    """Normalised uptake x cell-area weights on the uptake field's grid."""
    w = uptake * cell_area(uptake.latitude)
    return w / w.sum(skipna=True)


def extract_land(field, weights):
    """Uptake-weighted sum of `field`. NaN (land) cells count as zero.

    `weights` are normalised on the full grid; cells outside `field`
    (for example south of the kernel domain) drop out, as a missing
    contribution, not by renormalising.
    """
    w = weights.sel(latitude=field.latitude, longitude=field.longitude,
                    method="nearest", tolerance=0.01)
    w = w.assign_coords(latitude=field.latitude, longitude=field.longitude)
    return (field * w).sum(("latitude", "longitude"), skipna=True)
