#!/usr/bin/env python3
"""Build the Meltwater Discharge Explorer (v2).

Left: the surface d18O anomaly spreading over the North Atlantic and Arctic,
decade by decade, from the convolution of a meltwater forcing with the
HadCM3 impulse responses, over the exact HadCM3 input regions.
Right: the anomaly at one proxy site, split by source region, or the
regional discharge itself.

Inputs (made by the pipeline; see CLAUDE.md):
  outputs/convolution[/placeholder]/site_anomaly.csv.gz, run_info.json
  app/frames/meta.json and app/frames/*.txt   (scripts/export_app_frames.py)
  app/regions_hadcm3.geojson                  (scripts/make_region_geojson.py)
  data/Discharge_*.csv, myconfig/SITES.py, myconfig/PATHWAYS.py
  scripts/make_forward_model_schematic.py     ("How the model works" panel)

Outputs
  app/discharge_explorer.html   full HTML document. Serve app/ over http to
                                see the map (the frames are fetched):
                                python3 -m http.server 8765 --bind 127.0.0.1
  --publish-out PATH            second copy without the document wrapper, as
                                the claude.ai Artifact format expects. Refused
                                while the placeholder kernels are in use.

Standard library only. Run from the repository root:
  python3 app/build_app.py [--publish-out PATH]
"""
import argparse
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
APP = ROOT / "app"
OUTPUTS = ROOT / "outputs/convolution"
sys.path.insert(0, str(ROOT))

from myconfig.FORCINGS import FORCINGS          # noqa: E402  (plain dicts, no dependencies)
from myconfig.PATHWAYS import PATHWAYS, PATHWAY_LABELS  # noqa: E402
from myconfig.SITES import SITES, UPTAKE_RUN    # noqa: E402

CODES = ["Med", "Bri", "Fen", "EurArc", "AmeArc", "GIS", "NLau", "SLau", "GulofMex"]
SHORT = {"GulofMex": "GoM"}
SERIES_WINDOW = (24000, 8000)   # yr BP shown in the chart
SERIES_STEP = 20                # yr; the pipeline is decadal, every 2nd step is kept
SERIES_UNIT = 1e-4              # per mil; series are stored as integers of this unit


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def rnd(x, sig=5):
    return float(f"{x:.{sig}g}")


def discharge(key):
    e = FORCINGS[key]
    rows = sorted(read_rows(ROOT / e["file"]), key=lambda r: float(r[e["time"]]))
    ages, Q = [], {c: [] for c in CODES}
    for r in rows:
        age = float(r[e["time"]]) * e["to_bp"]
        if age < 0:
            continue
        ages.append(round(age))
        for c in CODES:
            Q[c].append(rnd(float(r[e["discharge"].format(code=c)])))
    step = ages[1] - ages[0] if ages[0] < ages[1] else ages[0] - ages[1]
    return {"label": e["label"], "ages": ages, "Q": Q, "step": abs(step)}


def site_series(src):
    """series[forcing][pathway][site][code or 'total'] -> ints (SERIES_UNIT), oldest first."""
    old, young = SERIES_WINDOW
    out = {}
    with gzip.open(src / "site_anomaly.csv.gz", "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            t = int(float(r["time_bp"]))
            if t > old or t < young or (old - t) % SERIES_STEP:
                continue
            node = out.setdefault(r["forcing"], {}).setdefault(r["pathway"], {}).setdefault(r["site"], {})
            for k in CODES + ["total"]:
                node.setdefault(k, []).append(round(float(r[k]) / SERIES_UNIT))
    return out


def howto():
    """The forward-model schematic as a scoped fragment (scripts/make_forward_model_schematic.py)."""
    spec = importlib.util.spec_from_file_location("fm_schematic", ROOT / "scripts/make_forward_model_schematic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build(fragment=True)


def build():
    real = OUTPUTS / "run_info.json"
    src = OUTPUTS if real.exists() else OUTPUTS / "placeholder"
    info = json.loads((src / "run_info.json").read_text(encoding="utf-8"))
    meta = json.loads((APP / "frames/meta.json").read_text(encoding="utf-8"))
    if bool(meta["placeholder"]) != bool(info["placeholder"]):
        sys.exit("app/frames and the site series come from different kernels; rerun export_app_frames.py")

    approx = {f["properties"]["code"]: f["properties"]
              for f in json.loads((APP / "regions_approx.geojson").read_text(encoding="utf-8"))["features"]
              if "code" in f["properties"]}
    geo = json.loads((APP / "regions_hadcm3.geojson").read_text(encoding="utf-8"))
    regions = [{
        "code": c, "short": SHORT.get(c, c),
        "name_draft": approx[c]["name_draft"], "name_2026b": approx[c]["name_2026b"],
        "ice_sheet": approx[c]["ice_sheet"], "drainage": approx[c]["drainage"],
    } for c in CODES]

    schedule = {k: [[a, b, m] for a, b, m in v] for k, v in PATHWAYS.items()}
    payload = {
        "placeholder": bool(info["placeholder"]),
        "kernels": info["kernels"],
        "run": {"created": info.get("created"), "git_commit": info.get("git_commit")},
        "regions": regions,
        "geo": {"type": "FeatureCollection", "features": geo["features"]},
        "sites": [{"name": n, "lat": s["lat"], "lon": s["lon"], "type": s["type"], "domain": s["domain"]}
                  for n, s in SITES.items()],
        "uptake_run": UPTAKE_RUN,
        "forcings": {k: discharge(k) for k in FORCINGS},
        "pathways": {k: {"label": PATHWAY_LABELS[k], "schedule": schedule[k]} for k in PATHWAYS},
        "series": {"window": SERIES_WINDOW, "step": SERIES_STEP, "unit": SERIES_UNIT, "data": site_series(src)},
        "frames": meta,
        "land": json.loads((APP / "vendor" / "land-110m.json").read_text(encoding="utf-8")),
    }
    return payload


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--publish-out", type=Path, help="also write the Artifact copy (no document wrapper) here")
    args = ap.parse_args()

    payload = build()
    if args.publish_out and payload["placeholder"]:
        sys.exit("Refusing --publish-out: the outputs use PLACEHOLDER kernels. "
                 "Put the pulse-kernel bundle in data/kernels/ and rerun the pipeline first.")

    blob = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    template = (APP / "template.html").read_text(encoding="utf-8")
    assert template.count("/*__DATA__*/null") == 1, "data placeholder missing"
    page = template.replace("/*__DATA__*/null", blob)
    fm_css, fm_html = howto()
    assert page.count("/*__HOWTO_CSS__*/") == 1 and page.count("<!--__HOWTO__-->") == 1, "howto placeholders missing"
    page = page.replace("/*__HOWTO_CSS__*/", fm_css).replace("<!--__HOWTO__-->", fm_html)
    head, body = page.split("<!--HEAD-END-->", 1)

    full = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
            "<!-- Generated by app/build_app.py from app/template.html; do not edit by hand. -->\n"
            + head.strip() + "\n</head>\n<body>\n" + body.strip() + "\n</body>\n</html>\n")
    out = APP / "discharge_explorer.html"
    out.write_text(full, encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size / 1024:.0f} KB)"
          + ("  [PLACEHOLDER kernels: local use only]" if payload["placeholder"] else ""))

    if args.publish_out:
        args.publish_out.write_text(head.strip() + "\n" + body.strip() + "\n", encoding="utf-8")
        print(f"wrote {args.publish_out} ({args.publish_out.stat().st_size / 1024:.0f} KB)")
        print("publish with files: " + ", ".join(f"{p} <- app/{p}" for p in payload["frames"]["files"].values()))

    # Spot checks
    for k, r in payload["forcings"].items():
        tot = [sum(r["Q"][c][i] for c in CODES) for i in range(len(r["ages"]))]
        i = max(range(len(tot)), key=tot.__getitem__)
        print(f"{r['label']}: peak total discharge {tot[i]:.3f} Sv at {r['ages'][i] / 1000:.1f} ka")
    s = payload["series"]["data"]["glac1d"]["mixed"]["La Vallina (NISA)"]["total"]
    print(f"NISA mixed GLAC-1D: {len(s)} steps, min {min(s) * SERIES_UNIT:.3f} per mil")


if __name__ == "__main__":
    main()
