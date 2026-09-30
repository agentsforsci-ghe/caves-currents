#!/usr/bin/env python3
"""
Turn the decadal surface fields into compact frames for the explorer app.

Reads outputs/convolution[/placeholder]/fields/field_<forcing>_<pathway>.nc
(written by run_convolution.py) and writes, into app/frames/:

    <forcing>_<pathway>.txt   uint8 frames (frame-major: n_frames x n_cells),
                              gzip-compressed and base64-encoded, because
                              claude.ai artifacts serve text but not raw binary
    meta.json                 cell centres, frame ages, colour-scale limits,
                              kernel source and the placeholder flag

Only ocean cells with a value at every frame are kept. Each byte encodes the
freshening magnitude |min(A, 0)| on one logarithmic scale shared by all
files, so the six combinations compare directly:

    q = 0                          |A| < VMIN (and the rare, tiny A > 0)
    q = 1 + round(254 * log(|A|/VMIN) / log(VMAX/VMIN)),  capped at 255

Run from the repository root after run_convolution.py:
    python3 scripts/export_app_frames.py [--window 21500 10000]
"""

import argparse
import base64
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from myconfig.PATHS import OUTPUTS

OUT = ROOT / "app/frames"
VMIN = 0.005   # per mil; below this a cell shows as open ocean


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", nargs=2, type=int, default=[21500, 10000], metavar=("OLDEST_BP", "YOUNGEST_BP"))
    args = ap.parse_args()

    real, ph = OUTPUTS / "run_info.json", OUTPUTS / "placeholder/run_info.json"
    src = OUTPUTS if real.exists() else OUTPUTS / "placeholder"
    info = json.loads((src / "run_info.json").read_text())
    files = sorted((src / "fields").glob("field_*_*.nc")) if info["placeholder"] else \
        sorted((OUTPUTS / "fields").glob("field_*_*.nc"))
    if not files:
        sys.exit(f"No field files next to {src}; run scripts/run_convolution.py first.")

    old, young = args.window
    fields = {}
    for f in files:
        _, forcing, pathway = f.stem.split("_", 2)
        da = xr.open_dataarray(f).sel(time_bp=slice(old, young))
        fields[(forcing, pathway)] = da
    first = next(iter(fields.values()))
    ok = np.ones(first.shape[1:], bool)
    vmax = 0.0
    for da in fields.values():
        v = da.values
        ok &= np.isfinite(v).all(axis=0)
        vmax = max(vmax, float(np.nanmax(-v)))
    iy, ix = np.nonzero(ok)
    lat = first.latitude.values[iy]
    lon = first.longitude.values[ix]

    OUT.mkdir(parents=True, exist_ok=True)
    for old in list(OUT.glob("*.bin")) + list(OUT.glob("*.txt")):
        old.unlink()
    span = np.log(vmax / VMIN)
    sizes = {}
    for (forcing, pathway), da in fields.items():
        mag = np.clip(-da.values[:, iy, ix], 0, None)
        q = np.zeros(mag.shape, np.uint8)
        on = mag >= VMIN
        q[on] = np.minimum(255, 1 + np.rint(254 * np.log(mag[on] / VMIN) / span)).astype(np.uint8)
        path = OUT / f"{forcing}_{pathway}.txt"
        path.write_bytes(base64.b64encode(gzip.compress(q.tobytes(), compresslevel=9, mtime=0)))
        sizes[path.name] = path.stat().st_size

    meta = {
        "cells": {"lat": [round(float(a), 3) for a in lat], "lon": [round(float(a), 3) for a in lon], "d": 1.25},
        "ages": [int(a) for a in first.time_bp.values],
        "files": {f"{fo}_{pw}": f"frames/{fo}_{pw}.txt" for fo, pw in fields},
        "encoding": "gzip+base64",
        "scale": {"vmin": VMIN, "vmax": round(vmax, 4), "levels": 255, "type": "log"},
        "placeholder": bool(info["placeholder"]),
        "kernels": info["kernels"],
        "source_run": info.get("created"),
        "git_commit": info.get("git_commit"),
    }
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    n_frames, n_cells = len(meta["ages"]), len(lat)
    print(f"{len(fields)} files, {n_frames} frames x {n_cells} cells, "
          f"{max(sizes.values()) / 1e6:.1f} MB each (gzip+base64), scale {VMIN}-{vmax:.2f} per mil"
          + ("  [PLACEHOLDER kernels]" if meta["placeholder"] else ""))


if __name__ == "__main__":
    main()
