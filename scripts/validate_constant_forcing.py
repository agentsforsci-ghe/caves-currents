#!/usr/bin/env python3
"""
Validate the convolution against the meltmodel paper (Endres et al. 2026b).

Each region's source-region d18O anomaly of a paper scenario (-35 per mil
end-member, regionalmeltdischarge_withd18O.pkl) is held constant for 500 yr
and run through the same convolution code as the pipeline
(mymodules.convolution.convolve_series), with the real pulse kernels. After
500 yr the site anomaly should equal the paper's equilibrium site anomaly
(proxymag.pkl, the values behind its Fig. 7).

Writes outputs/convolution/validation_constant_forcing.csv with, per scenario,
mode and site: the model and paper totals and their ratio for factor 1 and
for factor 10.

Result on 2026-10-01: factor 1 reproduces the paper (median ratio 0.90-1.02
per scenario; 18.2 ka cold within 0.96-1.01 at all nine sites); factor 10
gives ten times the paper. Laura chose factor 1.

Run from the repository root:
    python3 scripts/validate_constant_forcing.py
"""

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore", category=UserWarning)

from myconfig.DYES import CODES
from myconfig.PATHS import MELTMODEL_DATA, OUTPUTS
from myconfig.SITES import SITES
from mymodules.convolution import convolve_series
from mymodules.kernels import Kernels

SCENARIOS = {"17.8k": "17.8ka", "18.2k": "18.2 ka", "19.4k": "19.4 ka", "20.7k": "20.7 ka"}
MODES = ["cold", "zonal"]          # the paper has no meridional site results
YEARS = 500


def main():
    pm = pickle.load(open(MELTMODEL_DATA / "intermediates/dyestuff_modelpaper/proxymag.pkl", "rb"))
    gdg = pd.read_pickle(MELTMODEL_DATA / "intermediates/regionalmeltdischarge_withd18O.pkl")
    kern = Kernels.load("bundle")
    points = kern.points()
    n = YEARS // 10
    rows = []
    for scen, col in SCENARIOS.items():
        # One value per region, in dye order (the paper's table is in dye order).
        f0 = gdg[f"mean (-35.0) region d18O anomaly {col}"].to_numpy(float)
        forcing = pd.DataFrame(np.tile(f0, (n, 1)), columns=CODES,
                               index=pd.Index(np.arange(-YEARS, 0, 10), name="t"))
        for mode in MODES:
            for site in SITES:
                if site not in pm[scen][mode]:
                    continue
                a = convolve_series(forcing, points.sel(site=site), mode, scale=1.0)[-1].sum()
                paper = pm[scen][mode][site]["total"]
                rows.append({"scenario": scen, "mode": mode, "site": site, "model_x1": a, "paper": paper,
                             "ratio_x1": a / paper, "ratio_x10": 10 * a / paper})
    d = pd.DataFrame(rows)
    out = OUTPUTS / "validation_constant_forcing.csv"
    d.to_csv(out, index=False, float_format="%.5g")
    s = d.groupby(["scenario", "mode"]).agg(x1_median=("ratio_x1", "median"), x1_min=("ratio_x1", "min"),
                                            x1_max=("ratio_x1", "max"), x10_median=("ratio_x10", "median"))
    print(s.round(2).to_string())
    print(f"\nAll: factor 1 median ratio {d.ratio_x1.median():.2f}, factor 10 median ratio {d.ratio_x10.median():.2f}")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
