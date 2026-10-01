"""The nine dye-tracer input regions, in dye order (dye00 ... dye08).

`code` is the column prefix in the discharge files, `region` the name used in
Endres et al. (2026b). Colours match `reg_col` in nisa_meltwater_sources.qmd
and the --r-* tokens in app/template.html (keep all three in sync).
"""

import pandas as pd

DYE_TABLE = pd.DataFrame([
    {"dye": "dye00", "code": "Med",      "region": "EIS_MedSea",       "colour": "#1f77b4"},
    {"dye": "dye01", "code": "Bri",      "region": "EIS_BayOfBiscay",  "colour": "#f1780b"},
    {"dye": "dye02", "code": "Fen",      "region": "EIS_NorwegianSea", "colour": "#2ca02c"},
    {"dye": "dye03", "code": "EurArc",   "region": "EIS_Arctic",       "colour": "#d62728"},
    {"dye": "dye04", "code": "AmeArc",   "region": "LAU_Arctic",       "colour": "#9467bd"},
    {"dye": "dye05", "code": "GIS",      "region": "GIS_GreenlandSea", "colour": "#9b4c3c"},
    {"dye": "dye06", "code": "NLau",     "region": "LAU_LabradorSea",  "colour": "#de72bd"},
    {"dye": "dye07", "code": "SLau",     "region": "LAU_StLawrence",   "colour": "#a5a608"},
    {"dye": "dye08", "code": "GulofMex", "region": "LAU_GulfOfMexico", "colour": "#00b3c3"},
])

DYES = list(DYE_TABLE["dye"])
CODES = list(DYE_TABLE["code"])
CODE_TO_REGION = dict(zip(DYE_TABLE["code"], DYE_TABLE["region"]))
