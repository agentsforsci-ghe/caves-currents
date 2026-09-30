#!/usr/bin/env python3
"""
Convolve meltwater forcings with the HadCM3 impulse responses.

For every forcing (myconfig/FORCINGS.py), AMOC pathway (myconfig/PATHWAYS.py)
and proxy site (myconfig/SITES.py) this writes the regional contributions to
the surface d18O anomaly, and for every forcing and pathway the total surface
anomaly field at decadal resolution.

Outputs (in outputs/convolution/, or outputs/convolution/placeholder/ when
the placeholder kernels are used):
    site_anomaly.csv.gz   forcing, pathway, site, time_bp, one column per
                          region code, total  (per mil)
    scaling_check.csv     500-yr constant-input response vs the equilibrium
                          weights in data/meltmodel_site_weights.csv
    fields/field_<forcing>_<pathway>.nc
                          total anomaly (time_bp, latitude, longitude) over
                          FIELD_DOMAIN; large, not in git
    run_info.json         kernel source, settings, git commit

Run from the repository root:
    python3 scripts/run_convolution.py                   # everything
    python3 scripts/run_convolution.py --kernels placeholder --no-fields
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from myconfig.DYES import CODE_TO_REGION, CODES
from myconfig.FORCINGS import FORCINGS
from myconfig.PATHS import FIELDS, OUTPUTS, SITE_WEIGHTS
from myconfig.PATHWAYS import FIELD_DOMAIN, PATHWAYS
from myconfig.SITES import SITES
from mymodules.convolution import SCALE, convolve_field, convolve_series
from mymodules.forcing import load_forcing
from mymodules.kernels import Kernels


def crop_domain(field, dom):
    """Crop a 0-360 longitude field to FIELD_DOMAIN, with longitudes in -180..180."""
    lon = ((field.longitude + 180) % 360) - 180
    field = field.assign_coords(longitude=lon).sortby("longitude")
    return field.sel(latitude=slice(dom["lat_min"], dom["lat_max"]),
                     longitude=slice(dom["lon_min"], dom["lon_max"]))


def scaling_check(points):
    """Response at 500 yr to a constant unit input, per site/mode/region."""
    step = (SCALE * points.sum("lag")).to_dataframe(name="response_500yr").reset_index()
    step["region"] = step["dye"].map(dict(zip([f"dye{i:02d}" for i in range(9)],
                                              [CODE_TO_REGION[c] for c in CODES])))
    ref = pd.read_csv(SITE_WEIGHTS).rename(columns={"amoc_mode": "mode", "weight": "equilibrium_weight"})
    out = step.merge(ref[["site", "mode", "region", "equilibrium_weight"]],
                     on=["site", "mode", "region"], how="left")
    out["ratio"] = out["response_500yr"] / out["equilibrium_weight"]
    return out[["site", "mode", "region", "response_500yr", "equilibrium_weight", "ratio"]]


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernels", choices=["auto", "bundle", "placeholder"], default="auto")
    ap.add_argument("--forcing", nargs="+", default=list(FORCINGS), help="keys of FORCINGS")
    ap.add_argument("--pathway", nargs="+", default=list(PATHWAYS))
    ap.add_argument("--no-fields", action="store_true", help="skip the surface fields")
    ap.add_argument("--field-window", nargs=2, type=int, default=[23000, 9000],
                    metavar=("OLDEST_BP", "YOUNGEST_BP"), help="years BP kept in the field files")
    args = ap.parse_args()

    t0 = time.time()
    kern = Kernels.load(args.kernels)
    outdir = OUTPUTS / "placeholder" if kern.placeholder else OUTPUTS
    fielddir = outdir / "fields" if kern.placeholder else FIELDS
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Kernels: {kern.source}")
    if kern.placeholder:
        print("  !! PLACEHOLDER kernels: timing is invented. Outputs go to", outdir)

    points = kern.points()                                  # (site, mode, dye, lag)
    chk = scaling_check(points)
    chk.to_csv(outdir / "scaling_check.csv", index=False)
    r = chk["ratio"].dropna()
    print(f"Scaling check: response/equilibrium ratio median {r.median():.4f}, "
          f"range {r.min():.4f}-{r.max():.4f} over {len(r)} site-mode-region cells")

    rows = []
    for fk in args.forcing:
        F = load_forcing(fk)
        for pw in args.pathway:
            for site in SITES:
                A = convolve_series(F, points.sel(site=site), pw)
                df = pd.DataFrame(A, columns=CODES)
                df.insert(0, "time_bp", -F.index.to_numpy())
                df.insert(0, "site", site)
                df.insert(0, "pathway", pw)
                df.insert(0, "forcing", fk)
                df["total"] = df[CODES].sum(axis=1)
                rows.append(df)
            print(f"  series {fk:7s} {pw:6s} done")
    series = pd.concat(rows, ignore_index=True)
    series = series[series["time_bp"] >= 0]           # ICE-6G runs to +2 kyr model time
    num = CODES + ["total"]
    series[num] = series[num].astype("float32")
    series.to_csv(outdir / "site_anomaly.csv.gz", index=False, float_format="%.5g",
                  compression={"method": "gzip", "mtime": 0})

    field_files = []
    if not args.no_fields:
        fielddir.mkdir(parents=True, exist_ok=True)
        kf = crop_domain(kern.field, FIELD_DOMAIN)
        old, young = args.field_window
        for fk in args.forcing:
            F = load_forcing(fk)
            bp = -F.index.to_numpy()
            keep = (bp <= old) & (bp >= young)
            for pw in args.pathway:
                tot = convolve_field(F, kf, pw, keep=keep).astype("float32")
                da = xr.DataArray(tot, dims=("time_bp", "latitude", "longitude"),
                                  coords={"time_bp": bp[keep], "latitude": kf.latitude.values,
                                          "longitude": kf.longitude.values},
                                  name="anomaly",
                                  attrs={"units": "permil", "forcing": FORCINGS[fk]["label"],
                                         "pathway": pw, "kernels": kern.source,
                                         "placeholder": int(kern.placeholder)})
                path = fielddir / f"field_{fk}_{pw}.nc"
                da.to_netcdf(path, encoding={"anomaly": {"zlib": True, "complevel": 4}})
                field_files.append(str(path.relative_to(ROOT)))
                print(f"  field  {fk:7s} {pw:6s} {tot.shape} -> {path.name}")

    info = {
        "kernels": kern.source,
        "placeholder": kern.placeholder,
        "kernel_info": {k: (v if isinstance(v, (int, float, str)) else str(v)) for k, v in kern.info.items()},
        "forcings": {k: FORCINGS[k]["label"] for k in args.forcing},
        "pathways": {k: PATHWAYS[k] for k in args.pathway},
        "sites": list(SITES),
        "scale": SCALE,
        "field_window_bp": args.field_window,
        "field_domain": FIELD_DOMAIN,
        "field_files": field_files,
        "scaling_ratio_median": float(r.median()),
        "git_commit": git_commit(),
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "runtime_s": round(time.time() - t0, 1),
    }
    (outdir / "run_info.json").write_text(json.dumps(info, indent=2))
    print(f"Done in {info['runtime_s']} s. Outputs in {outdir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
