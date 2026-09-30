"""Impulse-response kernels: the exported HadCM3 pulse bundle or a placeholder.

A kernel is kernel(mode, dye, lag, latitude, longitude): the surface dye
response, in 10-yr means, to a 10-yr dye pulse in one input region. The
convolution (mymodules/convolution.py) multiplies the decadal forcing by 10,
as the 2025.05 notebook did.

The placeholder exists only to build and test the pipeline before the real
bundle is available. It spreads each mode's 500-yr equilibrium field
(meltmodel mean_dye_<mode>.nc) over time with a single first-order response
time TAU, so that a constant forcing reaches exactly the equilibrium weights
after 500 yr. Its timing is invented. Never publish results made with it.
"""

import hashlib
from pathlib import Path

import numpy as np
import xarray as xr

from myconfig.DYES import DYES
from myconfig.PATHS import KERNELS, MEAN_DYE, UPTAKE_WEIGHTS
from myconfig.SITES import SITES, UPTAKE_RUN
from mymodules import sites as st

MODES = ["cold", "zonal", "merid"]
N_LAG = 50
DECADE = 10
PLACEHOLDER_TAU = 100.0  # years


class Kernels:
    """Holds kernel fields and extracts point kernels at the proxy sites."""

    def __init__(self, field, source, placeholder, info=None):
        self.field = field            # DataArray (mode, dye, lag, latitude, longitude)
        self.source = source
        self.placeholder = placeholder
        self.info = info or {}
        self._uptake = None

    # ---- constructors ------------------------------------------------------
    @classmethod
    def from_bundle(cls, path=KERNELS):
        ds = xr.open_dataset(path)
        da = ds["kernel"].load()
        da = da.assign_coords(dye=[str(d) for d in da.dye.values],
                              mode=[str(m) for m in da.mode.values])
        info = dict(ds.attrs)
        info["sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        return cls(da, source=str(path), placeholder=False, info=info)

    @classmethod
    def placeholder_from_equilibrium(cls, tau=PLACEHOLDER_TAU, lat_min=0.0):
        t_end = np.arange(1, N_LAG + 1) * DECADE
        step = (1 - np.exp(-t_end / tau)) / (1 - np.exp(-N_LAG * DECADE / tau))
        h = np.diff(np.concatenate([[0.0], step])) / DECADE   # sums to 1/10
        fields = []
        for mode in MODES:
            ds = xr.open_dataset(str(MEAN_DYE).format(mode=mode))
            ds = ds.rename({"lat": "latitude", "lon": "longitude"})
            eq = xr.concat([ds[d] for d in DYES], dim="dye").assign_coords(dye=DYES)
            eq = eq.sel(latitude=slice(lat_min, None))
            fields.append(eq)
        eq = xr.concat(fields, dim="mode").assign_coords(mode=MODES)
        lag = xr.DataArray(h, coords={"lag": np.arange(N_LAG) * DECADE}, dims="lag")
        da = (eq * lag).transpose("mode", "dye", "lag", "latitude", "longitude").astype("float32")
        info = {"tau_years": tau, "note": "equilibrium field x first-order response; timing invented"}
        return cls(da.load(), source=f"placeholder (tau={tau:g} yr)", placeholder=True, info=info)

    @classmethod
    def load(cls, which="auto"):
        """'bundle', 'placeholder', or 'auto' (bundle if present)."""
        if which == "bundle" or (which == "auto" and Path(KERNELS).exists()):
            return cls.from_bundle()
        return cls.placeholder_from_equilibrium()

    # ---- site extraction ---------------------------------------------------
    def _uptake_weights(self):
        if self._uptake is None:
            self._uptake = xr.open_dataarray(UPTAKE_WEIGHTS).load()
        return self._uptake

    def point(self, site):
        """Point kernel (mode, dye, lag) at a site of myconfig.SITES."""
        s = SITES[site]
        if s["domain"] == "ocean":
            return st.extract_ocean(self.field, s["lat"], s["lon"], s.get("box", 2)).reset_coords(drop=True)
        up = self._uptake_weights()
        per_mode = []
        for mode in self.field.mode.values:
            w = st.land_weights(up.sel(site=s["uptake"], run=UPTAKE_RUN[str(mode)]))
            per_mode.append(st.extract_land(self.field.sel(mode=mode), w).reset_coords(drop=True))
        return xr.concat(per_mode, dim="mode").assign_coords(mode=self.field.mode.values)

    def points(self, site_names=None):
        names = list(site_names or SITES)
        return xr.concat([self.point(n) for n in names], dim="site").assign_coords(site=names)
