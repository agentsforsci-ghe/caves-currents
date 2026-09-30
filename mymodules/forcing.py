"""Load a meltwater forcing and put it on the 10-yr convolution grid."""

import numpy as np
import pandas as pd

from myconfig.DYES import CODES
from myconfig.FORCINGS import FORCINGS
from myconfig.PATHS import ROOT

STEP = 10  # years; the impulse responses are decadal


def _entry(forcing):
    return FORCINGS[forcing] if isinstance(forcing, str) else forcing


def load_forcing(forcing, kind="column"):
    """Source-region series on a 10-yr grid.

    Parameters
    ----------
    forcing : str or dict
        A key of FORCINGS, or an entry like those in FORCINGS.
    kind : "column" or "discharge"
        The d18O anomaly (the convolution input) or the discharge in Sv.

    Returns
    -------
    pd.DataFrame
        Index `t` = model time in years (negative = BP, ascending), one column
        per region code. Values are interpolated linearly between the native
        time steps, as in the 2025.05 convolution notebook.
    """
    e = _entry(forcing)
    df = pd.read_csv(ROOT / e["file"])
    t = np.asarray(df[e["time"]], dtype=float) * e["to_bp"] * -1  # model time
    cols = {code: e[kind].format(code=code) for code in CODES}
    missing = [c for c in cols.values() if c not in df.columns]
    if missing:
        raise KeyError(f"{e['label']}: missing columns {missing}")
    native = pd.DataFrame({code: df[col].to_numpy(float) for code, col in cols.items()}, index=t)
    native = native.sort_index()
    grid = np.arange(np.ceil(native.index.min() / STEP) * STEP, native.index.max() + 1, STEP)
    out = pd.DataFrame(
        {c: np.interp(grid, native.index, native[c]) for c in CODES},
        index=pd.Index(grid.astype(int), name="t"),
    )
    out.attrs["label"] = e["label"]
    return out
