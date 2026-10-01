"""Meltwater forcing series.

Each entry says how to turn a file into a source-region d18O anomaly per
region on a model-time axis (years, negative = BP). To add a forcing, add an
entry here. Any CSV with a `time_bp` column and one column per region code
(Med, Bri, ..., GulofMex) also works with the defaults of `generic_csv`.

    file        CSV path, relative to the repository root
    time        name of the time column
    to_bp       factor that turns the time column into years BP
    column      template for a region's anomaly column ({code} = region code)
    discharge   template for the region's discharge column in Sv (optional)
    label       name shown in figures and the app
"""

FORCINGS = {
    "glac1d": {
        "file": "data/Discharge_Glac1d_regional-withd18O.csv",
        "time": "time",
        "to_bp": -1,
        "column": "{code} d18O (-35.0)",
        "discharge": "{code}",
        "label": "GLAC-1D",
    },
    "ice6g": {
        "file": "data/Discharge_ice6g_regional-withd18O.csv",
        "time": "t_adj",
        "to_bp": -1,
        "column": "{code} d18O (-35.0)",
        "discharge": "discharge_{code}",
        "label": "ICE-6G",
    },
}


def generic_csv(path, label=None):
    """Forcing entry for a CSV with `time_bp` and one column per region code."""
    return {"file": str(path), "time": "time_bp", "to_bp": 1,
            "column": "{code}", "discharge": None, "label": label or str(path)}
