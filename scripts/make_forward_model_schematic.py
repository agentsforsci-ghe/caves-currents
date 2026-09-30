#!/usr/bin/env python3
"""
One-page schematic of the meltwater forward model.

Five numbered stages (forcing, pulse response, convolution, AMOC pathway,
reading at a site) and the equation that chains them. The forcing panel
draws the real GLAC-1D regional discharge; the pulse-response and pathway
curves are schematic shapes, not model output.

Used two ways:
  * app/build_app.py imports build(fragment=True) and inlines the result as
    the explorer's "How the model works" panel;
  * run as a script, it writes a standalone page for the claude.ai Artifact
    https://claude.ai/artifact/CqfLtTKp2MCwDi9L8YFa6x (no document wrapper):

        python3 scripts/make_forward_model_schematic.py --out /tmp/forward_model.html

All classes and SVG ids carry an `fm-` prefix, so the fragment cannot clash
with the explorer's own styles. Colours are the explorer's tokens. Standard
library only.
"""
import argparse
import csv
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

STACK = ["Bri", "NLau", "SLau", "GulofMex", "EurArc", "Med", "Fen", "AmeArc", "GIS"]
SHORT = {"GulofMex": "GoM"}


def f(v):
    return f"{v:.1f}"


# ---------------------------------------------------------------- panel 1: real forcing
def panel_forcing():
    W, H, L, R, T, B = 300, 190, 34, 10, 16, 30
    rows = list(csv.DictReader(open(REPO / "data/Discharge_Glac1d_regional-withd18O.csv")))
    pts = []
    for r in rows:
        age = -float(r["time"])
        if 8000 <= age <= 24000:
            pts.append((age, {c: float(r[c]) for c in STACK}))
    pts.sort(key=lambda p: -p[0])
    x = lambda a: L + (24000 - a) / 16000 * (W - L - R)
    ymax = 0.28
    y = lambda v: H - B - v / ymax * (H - B - T)
    parts, base = [], [0.0] * len(pts)
    for c in STACK:
        top = [base[i] + pts[i][1][c] for i in range(len(pts))]
        up = " ".join(f"{f(x(pts[i][0]))},{f(y(top[i]))}" for i in range(len(pts)))
        down = " ".join(f"{f(x(pts[i][0]))},{f(y(base[i]))}" for i in reversed(range(len(pts))))
        parts.append(f'<polygon class="r-{c}" points="{up} {down}"/>')
        base = top
    peak = max(range(len(pts)), key=lambda i: base[i])
    pa, pv = pts[peak][0], base[peak]
    ticks = "".join(
        f'<line class="tick" x1="{f(x(a))}" x2="{f(x(a))}" y1="{H-B}" y2="{H-B+4}"/>'
        f'<text class="ax" x="{f(x(a))}" y="{H-B+15}" text-anchor="middle">{a//1000}</text>'
        for a in (24000, 20000, 16000, 12000, 8000))
    yt = "".join(
        f'<line class="grid" x1="{L}" x2="{W-R}" y1="{f(y(v))}" y2="{f(y(v))}"/>'
        f'<text class="ax" x="{L-5}" y="{f(y(v)+3.5)}" text-anchor="end">{int(v*1000)}</text>'
        for v in (0.1, 0.2))
    return f'''<svg viewBox="0 0 {W} {H}" role="img" aria-label="GLAC-1D meltwater discharge by region, 24 to 8 ka, peaking at 260 mSv at 14.4 ka">
  {yt}
  <g class="stack">{''.join(parts)}</g>
  <line class="axis" x1="{L}" x2="{W-R}" y1="{H-B}" y2="{H-B}"/>
  {ticks}
  <text class="ax" x="{W-R}" y="{H-3}" text-anchor="end">ka BP</text>
  <text class="ax" x="4" y="{T-4}">mSv</text>
  <line class="lead" x1="{f(x(pa)-2)}" y1="{f(y(pv))}" x2="{f(x(pa)-18)}" y2="{f(y(pv))}"/>
  <text class="note" x="{f(x(pa)-21)}" y="{f(y(pv)+3.5)}" text-anchor="end">MWP-1A, {pv*1000:.0f} mSv, {pa/1000:.1f} ka</text>
</svg>'''


# ------------------------------------------------------------- kernels (schematic)
def kern(t, a, b):
    return (t ** a) * math.exp(-t / b) if t > 0 else 0.0


MODES = {"cold": (1.6, 95.0), "zonal": (1.3, 45.0), "merid": (1.4, 65.0)}


def panel_kernel():
    W, H, L, R, T, B = 300, 190, 16, 12, 22, 30
    x0 = 96
    x = lambda t: x0 + t / 500 * (W - R - x0)
    lines, labels = [], []
    ts = range(0, 501, 5)
    for mode, (a, b) in MODES.items():
        vals = [kern(t, a, b) for t in ts]
        m = max(vals)
        y = lambda v: H - B - v / m * (H - B - T - 8)
        lines.append(f'<polyline class="k k-{mode}" points="' + " ".join(f"{f(x(t))},{f(y(v))}" for t, v in zip(ts, vals)) + '"/>')
        # Line legend in the empty top-right corner (the peaks are too close for direct labels).
        ly = T + 4 + 13 * len(labels)
        labels.append(f'<line class="k k-{mode}" x1="{W-R-70}" x2="{W-R-48}" y1="{ly-3.5}" y2="{ly-3.5}"/>'
                      f'<text class="note" x="{W-R-42}" y="{ly}">{mode}</text>')
    ticks = "".join(
        f'<line class="tick" x1="{f(x(t))}" x2="{f(x(t))}" y1="{H-B}" y2="{H-B+4}"/>'
        f'<text class="ax" x="{f(x(t))}" y="{H-B+15}" text-anchor="middle">{t}</text>' for t in (0, 250, 500))
    return f'''<svg viewBox="0 0 {W} {H}" role="img" aria-label="A 10-year dye pulse in one region produces a delayed, spread-out response at a site, different for each AMOC mode">
  <defs><marker id="ah2" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" class="head"/></marker></defs>
  <text class="note" x="{L}" y="{H-B-84}">10-yr</text>
  <text class="note" x="{L}" y="{H-B-72}">pulse</text>
  <rect class="pulse" x="{L+4}" y="{H-B-62}" width="14" height="62"/>
  <text class="note" x="{(L+22+x0-8)/2}" y="{H-B-38}" text-anchor="middle">HadCM3</text>
  <line class="flow" x1="{L+24}" y1="{H-B-30}" x2="{x0-8}" y2="{H-B-30}" marker-end="url(#ah2)"/>
  {''.join(lines)}
  {''.join(labels)}
  <line class="axis" x1="{x0}" x2="{W-R}" y1="{H-B}" y2="{H-B}"/>
  {ticks}
  <text class="ax" x="{W-R}" y="{H-3}" text-anchor="end">years after pulse</text>
  <text class="ax" x="{x0}" y="{T-9}">h(k): dye at one site</text>
</svg>'''


# --------------------------------------------------------- panel 3: superposition
def panel_convolution():
    W, H, L, R, T, B = 300, 190, 14, 12, 34, 30
    F = [0.15, 0.2, 0.2, 0.9, 1.0, 0.45, 0.25, 0.15, 0.1, 0.1]
    n = 34
    a, b = 1.3, 2.2          # kernel in decades, zonal-like
    h = [kern(k, a, b) for k in range(n)]
    hs = sum(h)
    h = [v / hs for v in h]
    resp = [sum(F[j] * h[i - j] for j in range(len(F)) if 0 <= i - j < n) for i in range(n)]
    x = lambda i: L + i / (n - 1) * (W - L - R)
    base = H - B
    fb = 38        # forcing bar height scale
    zero = base - 50
    rmax = max(resp)
    yr = lambda v: zero - v / rmax * (zero - T)
    bars = "".join(
        f'<rect class="fbar" x="{f(x(i)-3)}" y="{f(base - v*fb)}" width="6" height="{f(v*fb)}"/>'
        for i, v in enumerate(F))
    copies, cmax, ic = [], 0.0, 0
    for j, v in enumerate(F):
        copies.append('<polyline class="copy" points="' + " ".join(f"{f(x(j+k))},{f(yr(v*h[k]))}" for k in range(n - j)) + '"/>')
    kp = max(range(n), key=lambda k: h[k])
    cmax = max(v * h[kp] for v in F)
    tot = " ".join(f"{f(x(i))},{f(yr(r))}" for i, r in enumerate(resp))
    ipk = max(range(n), key=lambda i: resp[i])
    il = ipk + 4
    ypk = yr(resp[ipk])
    return f'''<svg viewBox="0 0 {W} {H}" role="img" aria-label="Each decade of forcing launches a scaled copy of the pulse response; their sum is the site anomaly, delayed and smoothed">
  <defs><marker id="ah3" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" class="head"/></marker></defs>
  <line class="grid" x1="{L}" x2="{W-R}" y1="{zero}" y2="{zero}"/>
  {''.join(copies)}
  <polyline class="sum" points="{tot}"/>
  {bars}
  <line class="axis" x1="{L}" x2="{W-R}" y1="{base}" y2="{base}"/>
  <line class="lead" x1="{f(x(4))}" x2="{f(x(4))}" y1="{f(base - F[4]*fb - 3)}" y2="{f(ypk - 10)}" stroke-dasharray="2 3"/>
  <line class="flow" x1="{f(x(4))}" y1="{f(ypk - 10)}" x2="{f(x(ipk) - 2)}" y2="{f(ypk - 10)}" marker-end="url(#ah3)"/>
  <text class="note" x="{f(x(4) - 4)}" y="{f(ypk - 7)}" text-anchor="end">delay</text>
  <text class="note" x="{f(x(il) + 6)}" y="{f(yr(resp[il]) - 2)}">sum = A(t)</text>
  <text class="note" x="{f(x(len(F) + 3))}" y="{f(yr(cmax) - 2)}">copies, one per decade</text>
  <text class="note" x="{f(x(len(F)) + 4)}" y="{base-12}">F: forcing, per decade</text>
  <text class="ax" x="{W-R}" y="{H-3}" text-anchor="end">time (decades) →</text>
</svg>'''


# ---------------------------------------------------------------- panel 4: pathway
SCHEDULE = [(21500, 18000, "zonal"), (18000, 14700, "cold"), (14700, 13700, "merid"),
            (13700, 12600, "zonal"), (12600, 11300, "cold"), (11300, 10000, "zonal")]


def panel_pathway():
    W, H, L, R = 460, 214, 16, 16
    # Top: the whole mixed schedule. Bottom: a zoom around the switch to merid.
    x = lambda a: L + (22000 - a) / 12500 * (W - L - R)
    sy, sh = 18, 20
    segs = [(22000, 21500, "zonal")] + SCHEDULE + [(10000, 9500, "zonal")]
    rects = "".join(
        f'<rect class="m-{m}" x="{f(x(o))}" y="{sy}" width="{f(x(yng)-x(o))}" height="{sh}"/>'
        + (f'<text class="modelab" x="{f((x(o)+x(yng))/2)}" y="{sy+14}" text-anchor="middle">{m}</text>' if x(yng) - x(o) > 34 else "")
        for o, yng, m in segs)
    ticks = "".join(
        f'<text class="ax" x="{f(x(a))}" y="{sy+sh+12}" text-anchor="middle">{a/1000:g}</text>'
        for a in (21500, 18000, 14700, 12600, 10000))
    z0, z1 = 15300, 13900
    zx = lambda a: L + (z0 - a) / (z0 - z1) * (W - L - R)
    zt, zb = 92, 176          # zoom plot top / baseline
    zoom_band = (f'<rect class="zoomwin" x="{f(x(z0))}" y="{sy-3}" width="{f(x(z1)-x(z0))}" height="{sh+6}"/>'
                 f'<line class="lead" x1="{f(x(z0))}" y1="{sy+sh+3}" x2="{L}" y2="{zt-18}"/>'
                 f'<line class="lead" x1="{f(x(z1))}" y1="{sy+sh+3}" x2="{W-R}" y2="{zt-18}"/>')
    strip = (f'<rect class="m-cold" x="{f(zx(z0))}" y="{zb}" width="{f(zx(14700)-zx(z0))}" height="10"/>'
             f'<rect class="m-merid" x="{f(zx(14700))}" y="{zb}" width="{f(zx(z1)-zx(14700))}" height="10"/>')
    def release(rel, mode, amp):
        a_, b_ = MODES[mode]
        ts = list(range(0, 1301, 10))
        vals = [kern(t, a_, b_) for t in ts]
        mv = max(vals)
        pts = [(rel - t, zb - 4 - v / mv * amp) for t, v in zip(ts, vals) if rel - t >= z1]
        return pts, " ".join(f"{f(zx(a))},{f(yv)}" for a, yv in pts)
    p1, c1 = release(15100, "cold", (zb - zt) * 0.72)
    p2, c2 = release(14550, "merid", (zb - zt) * 0.62)
    zticks = "".join(
        f'<text class="ax" x="{f(zx(a))}" y="{zb+22}" text-anchor="middle">{a/1000:g}</text>' for a in (15200, 14700, 14200))
    return f'''<svg viewBox="0 0 {W} {H}" role="img" aria-label="The mixed pathway switches AMOC mode on a fixed schedule; water released before a switch keeps spreading with the kernel of its release mode">
  {rects}
  {ticks}
  {zoom_band}
  <text class="ax" x="{f((x(z0)+x(z1))/2)}" y="{zt-24}" text-anchor="middle">zoom: 15.3–13.9 ka</text>
  <line class="switch" x1="{f(zx(14700))}" x2="{f(zx(14700))}" y1="{zt-6}" y2="{zb+10}"/>
  <text class="note" x="{f(zx(14700)+4)}" y="{zt}">14.7 ka: switch to merid</text>
  <polyline class="k k-cold" points="{c1}"/>
  <circle class="rel" cx="{f(zx(15100))}" cy="{zb-4}" r="3.5"/>
  <text class="note" x="{L}" y="{zt+2}">released under cold:</text>
  <text class="note" x="{L}" y="{zt+14}">keeps the cold kernel</text>
  <polyline class="k k-merid" points="{c2}"/>
  <circle class="rel" cx="{f(zx(14550))}" cy="{zb-4}" r="3.5"/>
  <text class="note" x="{f(zx(14300))}" y="{f(zb-(zb-zt)*0.62-8)}">released under merid</text>
  {strip}
  {zticks}
  <text class="ax" x="{W-R}" y="{H-2}" text-anchor="end">ka BP</text>
</svg>'''


# ------------------------------------------------------------- panel 5: to the site
def panel_site():
    W, H = 460, 190
    return f'''<svg viewBox="0 0 {W} {H}" role="img" aria-label="The surface anomaly field is read at an ocean core as a box mean and at a land site as a moisture-uptake weighted sum">
  <defs>
    <radialGradient id="plume" cx="0.35" cy="0.45" r="0.65">
      <stop offset="0" class="pl0"/><stop offset="0.55" class="pl1"/><stop offset="1" class="pl2"/>
    </radialGradient>
    <marker id="ah5" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" class="head"/></marker>
  </defs>
  <rect class="sea" x="10" y="16" width="238" height="150" rx="3"/>
  <ellipse cx="100" cy="86" rx="128" ry="72" fill="url(#plume)"/>
  <path class="landmass" d="M196,16 L248,16 L248,166 L214,166 C206,140 222,120 208,100 C196,84 212,60 200,44 Z"/>
  <rect class="core" x="70" y="52" width="26" height="26"/>
  <circle class="dot" cx="83" cy="65" r="2.6"/>
  <text class="note" x="66" y="47" text-anchor="end">ocean core</text>
  <text class="note" x="66" y="59" text-anchor="end">±2° box</text>
  <ellipse class="uptake u1" cx="170" cy="118" rx="60" ry="32"/>
  <ellipse class="uptake u2" cx="182" cy="116" rx="36" ry="19"/>
  <ellipse class="uptake u3" cx="194" cy="114" rx="15" ry="9"/>
  <path class="dot" d="M218,108 l6,6 l-6,6 l-6,-6 z"/>
  <text class="note" x="229" y="101">land site</text>
  <text class="note" x="120" y="160">moisture uptake (trajectories)</text>
  <text class="ax" x="16" y="30">surface field A(x, t)</text>
  <line class="flow" x1="100" y1="56" x2="290" y2="56" marker-end="url(#ah5)"/>
  <line class="flow" x1="226" y1="114" x2="290" y2="124" marker-end="url(#ah5)"/>
  <text class="note" x="296" y="52">box mean</text>
  <text class="lbl" x="296" y="66">A_core(t)</text>
  <text class="note" x="296" y="120">uptake-weighted sum</text>
  <text class="lbl" x="296" y="134">A_site(t)</text>
  <text class="note" x="296" y="162">… for each of the 9 regions,</text>
  <text class="note" x="296" y="175">then summed</text>
</svg>'''



# ------------------------------------------------------------------ page parts
TOKENS_CSS = """
:root {
  --bg: #f3f6f8; --panel: #fbfcfd; --fg: #16212a; --fg-2: #4a5a67; --fg-3: #74838f; --rule: #d7dfe5; --grat: #d9e2ea;
  --accent: #1c5cab; --focus: #eb6834; --land: #dfe4e2; --land-edge: #b9c3c1; --ocean: #eef3f7;
  --m-cold: #c9d3db; --m-zonal: #e7ecf0; --m-merid: #aebbc6;
  --plume-0: #1c5cab; --plume-1: #9ec5f4;
  --r-Med: #1f77b4; --r-Bri: #f1780b; --r-Fen: #2ca02c; --r-EurArc: #d62728; --r-AmeArc: #9467bd; --r-GIS: #9b4c3c; --r-NLau: #de72bd; --r-SLau: #a5a608; --r-GulofMex: #00b3c3;
  --font-display: "Barlow Condensed", "Arial Narrow", "Roboto Condensed", sans-serif;
  --font-body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --font-data: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0f151a; --panel: #131a20; --fg: #e8eef2; --fg-2: #aab8c2; --fg-3: #7d8c97; --rule: #26313a; --grat: #1f2a33;
    --accent: #86b6ef; --focus: #f08a5d; --land: #28313a; --land-edge: #3c4853; --ocean: #151e25;
    --m-cold: #2b3640; --m-zonal: #1c252d; --m-merid: #3a4855;
    --plume-0: #6da7ec; --plume-1: #184f95;
    --r-Bri: #db6b00; --r-GIS: #a15142; --r-NLau: #cd63ad; --r-SLau: #959600; --r-GulofMex: #03a2b1;
    color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #0f151a; --panel: #131a20; --fg: #e8eef2; --fg-2: #aab8c2; --fg-3: #7d8c97; --rule: #26313a; --grat: #1f2a33;
  --accent: #86b6ef; --focus: #f08a5d; --land: #28313a; --land-edge: #3c4853; --ocean: #151e25;
  --m-cold: #2b3640; --m-zonal: #1c252d; --m-merid: #3a4855;
  --plume-0: #6da7ec; --plume-1: #184f95;
  --r-Bri: #db6b00; --r-GIS: #a15142; --r-NLau: #cd63ad; --r-SLau: #959600; --r-GulofMex: #03a2b1;
  color-scheme: dark;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg); font: 14px/1.5 var(--font-body); }
"""

# Component styles: every selector is anchored on a class, so prefix() scopes them all.
COMPONENT_CSS = """
.wrap { max-width: 1280px; margin: 0 auto; padding-inline: 16px; padding-block: 22px 32px; display: grid; gap: 16px; }
.body { display: grid; gap: 16px; min-width: 0; }
.top { display: grid; gap: 6px; }
.eyebrow { font: 500 11px/1 var(--font-body); letter-spacing: 0.08em; text-transform: uppercase; color: var(--fg-3); }
.top h1 { font: 600 32px/1.05 var(--font-display); letter-spacing: 0.01em; margin: 0; text-wrap: balance; }
.lede { margin: 0; color: var(--fg-2); max-width: 78ch; }
.stages { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 12px; }
.stage { grid-column: span 2; margin: 0; background: var(--panel); border: 1px solid var(--rule); border-radius: 6px; padding: 12px; display: grid; gap: 6px; align-content: start; min-width: 0; }
.wide { grid-column: span 3; }
@media (max-width: 1000px) { .stages { grid-template-columns: repeat(2, minmax(0, 1fr)); } .stage, .wide { grid-column: span 1; } .last { grid-column: span 2; } }
@media (max-width: 620px) { .stages { grid-template-columns: minmax(0, 1fr); } .stage, .wide, .last { grid-column: span 1; } }
.stage h2 { margin: 0; font: 600 18px/1.1 var(--font-display); letter-spacing: 0.02em; display: flex; gap: 8px; align-items: baseline; }
.n { font: 500 12px var(--font-data); color: var(--accent); border: 1px solid var(--accent); border-radius: 999px; padding: 0 7px; }
.stage svg { display: block; width: 100%; height: auto; max-width: 100%; }
.stage figcaption { font-size: 13px; color: var(--fg-2); display: grid; gap: 4px; }
.eq { font: 12.5px var(--font-data); color: var(--fg); }
.chain { background: var(--panel); border: 1px solid var(--rule); border-radius: 6px; padding: 14px 16px; display: grid; gap: 8px; min-width: 0; }
.chain h2 { margin: 0; font: 600 18px var(--font-display); letter-spacing: 0.02em; }
.big { font: 500 clamp(14px, 2.1vw, 21px)/1.5 var(--font-data); overflow-x: auto; white-space: nowrap; padding-block: 4px; }
.big sub { font-size: 0.7em; }
.terms { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 4px 20px; margin: 0; font-size: 13px; color: var(--fg-2); }
.terms div { display: grid; grid-template-columns: 5.5em 1fr; gap: 8px; }
.terms dt { font: 500 12.5px var(--font-data); color: var(--fg); }
.terms dd { margin: 0; }
.foot { font-size: 12px; color: var(--fg-3); max-width: 120ch; display: grid; gap: 4px; }
.ax { font: 10.5px var(--font-data); fill: var(--fg-3); }
.note { font: 11px var(--font-body); fill: var(--fg-2); }
.lbl { font: 600 13px var(--font-data); fill: var(--fg); }
.modelab { font: 500 11px var(--font-body); fill: var(--fg-2); }
.axis, .tick { stroke: var(--fg-3); stroke-width: 1; }
.grid { stroke: var(--grat); stroke-width: 1; }
.lead { stroke: var(--fg-3); stroke-width: 1; }
.flow { stroke: var(--fg-2); stroke-width: 1.4; fill: none; }
.head { fill: var(--fg-2); }
.stack polygon { stroke: none; }
""" + "".join(f".r-{c} {{ fill: var(--r-{c}); }}\n" for c in STACK) + """
.pulse, .fbar { fill: var(--accent); }
.k { fill: none; stroke-width: 2; stroke: var(--fg); }
.k-zonal { stroke-dasharray: 1.5 3; stroke: var(--fg-2); }
.k-merid { stroke-dasharray: 6 3; stroke: var(--fg-2); }
.k-cold { stroke: var(--fg); }
.copy { fill: none; stroke: var(--accent); stroke-width: 1; opacity: 0.45; }
.sum { fill: none; stroke: var(--fg); stroke-width: 2.2; }
.m-cold { fill: var(--m-cold); stroke: var(--panel); }
.m-zonal { fill: var(--m-zonal); stroke: var(--panel); }
.m-merid { fill: var(--m-merid); stroke: var(--panel); }
.zoomwin { fill: none; stroke: var(--fg-2); stroke-width: 1.2; }
.switch { stroke: var(--focus); stroke-width: 1.4; stroke-dasharray: 4 3; }
.rel { fill: var(--focus); }
.sea { fill: var(--ocean); stroke: var(--rule); }
.landmass { fill: var(--land); stroke: var(--land-edge); }
.pl0 { stop-color: var(--plume-0); stop-opacity: 0.85; }
.pl1 { stop-color: var(--plume-1); stop-opacity: 0.6; }
.pl2 { stop-color: var(--ocean); stop-opacity: 0; }
.core { fill: none; stroke: var(--fg); stroke-width: 1.6; }
.dot { fill: var(--focus); stroke: var(--panel); stroke-width: 1; }
.uptake { fill: none; stroke: var(--focus); stroke-dasharray: 3 3; }
.u2 { stroke-dasharray: none; opacity: 0.8; }
.u3 { stroke-dasharray: none; stroke-width: 1.6; }
"""

TOP = """<header class="top">
    <div class="eyebrow">Endres et al., in prep. · impulse-response forward model</div>
    <h1>From ice-sheet meltwater to a δ¹⁸O signal at a proxy site</h1>
    <p class="lede">LEDE</p>
  </header>"""
LEDE = ("The forward model has five steps. A reconstruction gives the meltwater per drainage region. "
        "The ocean model says how a short pulse from each region spreads. Adding up those pulses, decade "
        "by decade and in the AMOC mode of the time, gives the surface anomaly everywhere, which is then "
        "read at the core or cave.")


def body():
    return f"""<div class="stages">
    <figure class="stage">
      <h2><span class="n">1</span>Forcing</h2>
      {panel_forcing()}
      <figcaption>
        <span>Discharge from the ice-sheet reconstruction, per drainage region (GLAC-1D shown, 9 regions, 100-yr steps). It is turned into a source-region δ¹⁸O anomaly, interpolated to 10-yr steps.</span>
        <span class="eq">F<sub>r</sub>(t) = Q<sub>r</sub>(t) · Δt · δ¹⁸O<sub>ice</sub> / V<sub>r</sub>,  δ¹⁸O<sub>ice</sub> = −35 ‰</span>
      </figcaption>
    </figure>
    <figure class="stage">
      <h2><span class="n">2</span>Pulse response</h2>
      {panel_kernel()}
      <figcaption>
        <span>In HadCM3, 10 years of dye enter one region, then stop. The surface dye over the next 500 years is the region's impulse response. There is one per region and AMOC mode, and it is a map, shown here at one site. Curve shapes are schematic.</span>
        <span class="eq">h<sub>r,m</sub>(x, k),  k = 0 … 490 yr, decadal means</span>
      </figcaption>
    </figure>
    <figure class="stage">
      <h2><span class="n">3</span>Convolution</h2>
      {panel_convolution()}
      <figcaption>
        <span>Each decade of forcing launches a copy of the pulse response, scaled by its size. The copies overlap and add up. A short pulse arrives late and spread out, and a long input builds up to its equilibrium.</span>
        <span class="eq">A<sub>r</sub>(x, t) = Σ<sub>k</sub> 10 · F<sub>r</sub>(t − k) · h<sub>r</sub>(x, k)</span>
      </figcaption>
    </figure>
    <figure class="stage wide">
      <h2><span class="n">4</span>AMOC pathway</h2>
      {panel_pathway()}
      <figcaption>
        <span>Cold and zonal use one mode throughout. The mixed pathway follows the deglacial sequence on the strip; the zoom shows the switch to merid at 14.7 ka. Each decade's meltwater takes the kernel of the mode at its release and keeps it, so water released just before a switch still spreads the old way.</span>
      </figcaption>
    </figure>
    <figure class="stage wide last">
      <h2><span class="n">5</span>Reading it at a site</h2>
      {panel_site()}
      <figcaption>
        <span>Ocean cores take the mean of the surface field in a ±2° box. Caves and ice cores take the field weighted by where their rain picks up its moisture, from back-trajectories. Doing this for each region gives the contributions that are stacked in the explorer.</span>
      </figcaption>
    </figure>
  </div>
  <section class="chain" aria-label="The whole chain in one equation">
    <h2>The whole chain</h2>
    <div class="big">A<sub>site</sub>(t) = Σ<sub>r</sub> Σ<sub>x</sub> w<sub>site</sub>(x) · Σ<sub>k</sub> 10 · F<sub>r</sub>(t − k) · h<sub>r, m(t−k)</sub>(x, k)</div>
    <dl class="terms">
      <div><dt>F<sub>r</sub></dt><dd>source-region δ¹⁸O anomaly of region r (step 1)</dd></div>
      <div><dt>h<sub>r,m</sub></dt><dd>surface response to a 10-yr pulse, AMOC mode m (step 2)</dd></div>
      <div><dt>m(t − k)</dt><dd>mode at the time the water was released (step 4)</dd></div>
      <div><dt>w<sub>site</sub></dt><dd>±2° box for cores, moisture uptake for land sites (step 5)</dd></div>
      <div><dt>Σ<sub>k</sub>, Σ<sub>r</sub></dt><dd>over the 50 decades of the response, and over the 9 regions</dd></div>
      <div><dt>A<sub>site</sub></dt><dd>δ¹⁸O anomaly at the site, to compare with the proxy record</dd></div>
    </dl>
  </section>
  <footer class="foot">
    <div>Assumptions: the response is linear in discharge (no saturation), the kernels come from 17.8 ka boundary conditions, and only the surface layer is used. Pulse runs: xpran (cold), xprao (zonal) and xpujc (merid), each after 10 years of its constant-input parent (xpraj, xprak, xpram). The merid parent is not the paper's merid run, xpral.</div>
    <div>Forcing panel: GLAC-1D regional discharge, 24–8 ka, with regions in the explorer's colours. The pulse-response and pathway curves are schematic shapes, not model output.</div>
  </footer>"""


def prefix(text, css=False):
    """Give every class and SVG id an fm- prefix."""
    if css:
        return re.sub(r"\.([A-Za-z][\w-]*)", r".fm-\1", text)
    text = re.sub(r'class="([^"]*)"', lambda m: 'class="' + " ".join("fm-" + c for c in m.group(1).split()) + '"', text)
    text = re.sub(r'id="([^"]*)"', r'id="fm-\1"', text)
    return re.sub(r"url\(#([^)]*)\)", r"url(#fm-\1)", text)


def build(fragment=False):
    """Return (css, html). The fragment omits the page tokens and title block."""
    css = prefix(COMPONENT_CSS, css=True)
    if fragment:
        html = prefix(f'<div class="body"><p class="lede">{LEDE}</p>{body()}</div>')
        return css, html
    html = prefix(f'<div class="wrap">{TOP.replace("LEDE", LEDE)}{body()}</div>')
    return TOKENS_CSS + css, html


FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">\n'
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Barlow+Condensed:wght@500;600'
         '&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True, help="standalone page (Artifact format: no document wrapper)")
    args = ap.parse_args()
    css, html = build()
    page = (f"<title>Meltwater Forward Model</title>\n{FONTS}\n<style>\n"
            "/* Layout: five numbered stages in reading order (3 across, then 2 wide), closed by the one "
            "equation that chains them. Same tokens as the explorer app. */\n"
            f"{css}</style>\n{html}\n")
    args.out.write_text(page, encoding="utf-8")
    print(f"wrote {args.out} ({args.out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
