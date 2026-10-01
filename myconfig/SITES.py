"""Proxy sites (the nine of Endres et al. 2026b, Fig. 7).

Coordinates are copied from the meltmodel repo's myconfig/PROXYSITES.py.
Ocean sites use a latitude-weighted mean over a +/- `box` degree box around
the site. Land sites weight the surface field by the moisture-uptake field
of the trajectory run `uptake` (file prefix in the meltmodel data).
Llarga Cave is left out: it has no trajectory run.
"""

SITES = {
    "MD95-2010":           {"lat": 66.6833, "lon": -4.5667,  "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "PS2644-5":            {"lat": 67.8667, "lon": -21.7667, "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "MD03-2664":           {"lat": 57.4333, "lon": -48.6000, "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "ODP980":              {"lat": 55.4833, "lon": -14.7000, "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "MD95-2042":           {"lat": 37.8000, "lon": -10.1667, "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "MD02-2552":           {"lat": 26.9500, "lon": -91.3500, "type": "Marine sediment core", "domain": "ocean", "box": 2},
    "La Vallina (NISA)":   {"lat": 43.4167, "lon": -4.8000,  "type": "Speleothem", "domain": "land", "uptake": "NISA_LaVallina"},
    "Cave Without a Name": {"lat": 29.8667, "lon": -98.6333, "type": "Speleothem", "domain": "land", "uptake": "NonameCave"},
    "NGRIP":               {"lat": 75.1000, "lon": -42.3167, "type": "Ice core",   "domain": "land", "uptake": "NGRIP"},
}

# High-resolution trajectory runs per AMOC mode (meltmodel EXPERIMENTS.py).
# There is no merid trajectory run; the merid segment of the mixed pathway
# borrows the zonal (strong-AMOC) uptake field. This is an assumption.
UPTAKE_RUN = {"cold": "xqeic", "zonal": "xqeie", "merid": "xqeie"}
