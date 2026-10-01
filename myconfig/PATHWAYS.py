"""AMOC pathways: which mode's impulse response applies when.

`cold` and `zonal` use one mode throughout. `mixed` follows the schedule of
the 2025.05 convolution notebook (cell 73), in years BP. Segment edges are
half-open [older, younger). Before the first and after the last segment the
first and last modes continue, so the series is spun up from the start of
the forcing instead of starting from zero at 21.5 ka.
"""

MIXED_SCHEDULE = [
    (21500, 18000, "zonal"),
    (18000, 14700, "cold"),
    (14700, 13700, "merid"),
    (13700, 12600, "zonal"),
    (12600, 11300, "cold"),
    (11300, 10000, "zonal"),
]

PATHWAYS = {
    "cold": [(None, None, "cold")],
    "zonal": [(None, None, "zonal")],
    "mixed": MIXED_SCHEDULE,
}

PATHWAY_LABELS = {
    "cold": "Cold (weak AMOC)",
    "zonal": "Zonal (strong AMOC)",
    "mixed": "Mixed (deglacial AMOC sequence)",
}

# Map domain kept for the surface-field animation (degrees).
FIELD_DOMAIN = {"lat_min": 0.0, "lat_max": 90.0, "lon_min": -110.0, "lon_max": 60.0}
