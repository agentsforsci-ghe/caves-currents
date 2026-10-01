"""Where the pipeline reads and writes.

The meltmodel figure-scripts repository (Endres et al. 2026b; Zenodo
10.5281/zenodo.21828703) supplies the HadCM3 grid, the dye input regions,
the equilibrium dye fields and the moisture-uptake trajectories. Set the
environment variable MELTMODEL_REPO to use another copy.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

MELTMODEL_REPO = Path(os.environ.get(
    "MELTMODEL_REPO",
    Path.home() / "Documents/gitrepos/gh-lrndrs/gh-lrndrs-dyetracer_palaeo_figurescripts",
))
MELTMODEL_DATA = MELTMODEL_REPO / "data"
DYE_REGIONS = MELTMODEL_DATA / "intermediates/dyestuff_modelpaper/dye_regions_norm.nc"
MEAN_DYE = MELTMODEL_DATA / "intermediates/dyestuff_modelpaper/mean_dye_{mode}.nc"
TRAJECTORIES = MELTMODEL_DATA / "trajectories/{loc}_{exp}_th00_UTOT_weighted.nc"
LAND_SEA_MASK = MELTMODEL_DATA / "inputs/temev.qrparm.omask.nc"

KERNELS = DATA / "kernels/pulse_kernels_surface.nc"
UPTAKE_WEIGHTS = DATA / "meltmodel/uptake_weights.nc"
SITE_WEIGHTS = DATA / "meltmodel_site_weights.csv"

OUTPUTS = ROOT / "outputs/convolution"
FIELDS = OUTPUTS / "fields"          # large, gitignored
