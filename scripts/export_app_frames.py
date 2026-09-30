#!/usr/bin/env python3
"""
Turn the decadal surface fields into compact frames for the explorer app.

Reads outputs/convolution[/placeholder]/fields/field_<forcing>_<pathway>.nc
(written by run_convolution.py) and writes, into app/frames/:

    <forcing>_<pathway>.txt   uint8 frames (frame-major: n_frames x n_cells),
                              gzip-compressed and base64-encoded, because
                              claude.ai artifacts serve text but not raw binary
    diff_<pathway>.txt        same, for the difference between the first two
                              forcings (GLAC-1D minus ICE-6G), per pathway
    meta.json                 cell centres, frame ages, colour-scale limits,
                              kernel source and the placeholder flag

Only ocean cells with a value at every frame are kept. Each byte encodes the
freshening magnitude |min(A, 0)| on one logarithmic scale shared by all
files, so the six combinations compare directly:

    q = 0                          |A| < VMIN (and the rare, tiny A > 0)
    q = 1 + round(254 * log(|A|/VMIN) / log(VMAX/VMIN)),  capped at 255

The difference D = A(GLAC-1D) - A(ICE-6G) uses a symmetric log scale around
q = 128: |D| < DVMIN is 128, negative D (GLAC-1D fresher) is 128 - k and
positive D is 128 + k, with k = 1 ... 127 logarithmic from DVMIN to the
largest |D|.

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
DVMIN = 0.005  # per mil; differences smaller than this show as open ocean


def pack(q):
    """uint8 array -> gzip + base64 bytes (artifacts serve text, not binary)."""
    return base64.b64encode(gzip.compress(q.astype(np.uint8).tobytes(), compresslevel=9, mtime=0))


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
        path.write_bytes(pack(q))
        sizes[path.name] = path.stat().st_size

    # Differences between the first two forcings (in FORCINGS order), per pathway.
    from myconfig.FORCINGS import FORCINGS
    fa, fb = [k for k in FORCINGS if any(k == fo for fo, _ in fields)][:2]
    pathways = sorted({pw for _, pw in fields})
    diffs = {pw: fields[(fa, pw)].values[:, iy, ix] - fields[(fb, pw)].values[:, iy, ix] for pw in pathways}
    dvmax = max(float(np.abs(d).max()) for d in diffs.values())
    dspan = np.log(dvmax / DVMIN)
    for pw, d in diffs.items():
        mag = np.abs(d)
        k = np.zeros(d.shape)
        on = mag >= DVMIN
        k[on] = np.minimum(127, 1 + np.rint(126 * np.log(mag[on] / DVMIN) / dspan))
        q = 128 + np.sign(d) * k
        path = OUT / f"diff_{pw}.txt"
        path.write_bytes(pack(q))
        sizes[path.name] = path.stat().st_size

    meta = {
        "cells": {"lat": [round(float(a), 3) for a in lat], "lon": [round(float(a), 3) for a in lon], "d": 1.25},
        "ages": [int(a) for a in first.time_bp.values],
        "files": {f"{fo}_{pw}": f"frames/{fo}_{pw}.txt" for fo, pw in fields},
        "diff": {"a": fa, "b": fb, "files": {pw: f"frames/diff_{pw}.txt" for pw in pathways},
                 "scale": {"vmin": DVMIN, "vmax": round(dvmax, 4), "levels": 127, "type": "symlog"}},
        "encoding": "gzip+base64",
        "scale": {"vmin": VMIN, "vmax": round(vmax, 4), "levels": 255, "type": "log"},
        "placeholder": bool(info["placeholder"]),
        "kernels": info["kernels"],
        "source_run": info.get("created"),
        "git_commit": info.get("git_commit"),
    }
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    n_frames, n_cells = len(meta["ages"]), len(lat)
    print(f"{len(sizes)} files ({len(fields)} anomaly + {len(diffs)} difference), {n_frames} frames x {n_cells} cells, "
          f"{max(sizes.values()) / 1e6:.1f} MB each (gzip+base64), scale {VMIN}-{vmax:.2f}, difference {DVMIN}-{dvmax:.2f} per mil"
          + ("  [PLACEHOLDER kernels]" if meta["placeholder"] else ""))


if __name__ == "__main__":
    main()
