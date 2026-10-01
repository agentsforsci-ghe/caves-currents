#!/usr/bin/env python3
"""
Outlines of the nine HadCM3 dye input regions, as GeoJSON for the app.

Every ocean grid cell where a region's normalised dye input
(meltmodel dye_regions_norm.nc) is above zero becomes a 1.25 deg box; the
boxes of each region are merged. The result is the exact model input area,
not a hand-drawn outline. Properties match app/regions_approx.geojson
(`code`, `name_2026b`, ...), with `approximate: false`. Rings are wound
clockwise, as d3-geo expects, so the app needs no rewinding.

Run from the repository root:
    python3 scripts/make_region_geojson.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from shapely.geometry import MultiPolygon, Polygon, box, mapping
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from myconfig.DYES import DYE_TABLE
from myconfig.PATHS import DYE_REGIONS

OUT = ROOT / "app/regions_hadcm3.geojson"
APPROX = ROOT / "app/regions_approx.geojson"
D = 1.25


def d3_winding(geom):
    """Clockwise exterior rings, the convention d3-geo expects on the sphere.

    RFC 7946 asks for the opposite; the only consumer here is the d3 globe.
    """
    if isinstance(geom, Polygon):
        return orient(geom, sign=-1.0)
    return MultiPolygon([orient(p, sign=-1.0) for p in geom.geoms])


def main():
    ds = xr.open_dataset(DYE_REGIONS)
    lat = ds.latitude.values
    lon = ((ds.longitude.values + 180) % 360) - 180
    old = {f["properties"]["code"]: f["properties"] for f in json.load(open(APPROX))["features"]
           if "code" in f["properties"]}
    feats = []
    for _, row in DYE_TABLE.iterrows():
        m = ds[row["code"]].squeeze().values > 0
        iy, ix = np.nonzero(m)
        cells = [box(lon[j] - D / 2, lat[i] - D / 2, lon[j] + D / 2, lat[i] + D / 2) for i, j in zip(iy, ix)]
        geom = d3_winding(unary_union(cells).simplify(0.01))
        props = {k: v for k, v in old.get(row["code"], {}).items() if k != "drainage"}
        props.update({"code": row["code"], "name_2026b": row["region"], "n_cells": int(m.sum()),
                      "approximate": False,
                      "source": "HadCM3 dye input cells (dye_regions_norm.nc), Endres et al. 2026b"})
        feats.append({"type": "Feature", "properties": props, "geometry": mapping(geom)})
        print(f"  {row['code']:9s} {int(m.sum()):4d} cells  {geom.geom_type}")
    OUT.write_text(json.dumps({"type": "FeatureCollection", "features": feats}, separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
