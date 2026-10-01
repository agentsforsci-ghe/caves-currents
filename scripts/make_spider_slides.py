#!/usr/bin/env python3
"""
Three explanatory slides for the spider figures (claude.ai Slides deck).

1. explain    How to read a spider: the five spokes as plain questions,
              drawn with a real example (LAU_LabradorSea at NISA).
2. nisa       The nine NISA spiders (mixed pathway), GLAC-1D against ICE-6G.
3. takeaways  Five site-level results across the nine proxy sites.

All numbers come from the convolution outputs through R/site_metrics.R (the
same metrics as site_anomalies_v2.qmd), so the slides follow the data. The
spiders are inline SVG, drawn from the metrics, not pasted images.

The deck lives at https://claude.ai/artifact/TeDmKyFH5oj1oW9dpwhAEW (private
until shared; it downloads as PowerPoint or PDF). To update it:

    python3 scripts/make_spider_slides.py --out <folder>

then republish the three files `project/slides/{explain,nisa,takeaways}.html`
to that URL with `root` = <folder>. Pass `--index` only for a new deck: it also
writes `project/deck.json`, which the existing deck keeps.

Needs Rscript with the packages of R/site_metrics.R, and pandas.
"""
import argparse
import json
import math
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def load_metrics():
    """The five metrics per forcing x pathway x site x region, from R/site_metrics.R."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "metrics.csv"
        code = (f'source("R/site_metrics.R"); '
                f'write.csv(site_metrics(to_grid(read_site_anomaly())), "{out}", row.names = FALSE)')
        subprocess.run(["Rscript", "-e", code], cwd=REPO, check=True)
        return pd.read_csv(out)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out", type=Path, required=True, help="folder to write project/... into (the publish root)")
ap.add_argument("--index", action="store_true", help="also write project/deck.json (new deck only)")
args = ap.parse_args()
ROOT = args.out
(ROOT / "project/slides").mkdir(parents=True, exist_ok=True)

m = load_metrics()
KEYS = ["share", "dom", "sd", "vshare", "sd_diff"]
LABELS = {"share": "Share", "dom": "Dominant", "sd": "SD", "vshare": "Var. share", "sd_diff": "SD Δ"}
REGIONS = ["EIS_MedSea", "EIS_BayOfBiscay", "EIS_NorwegianSea", "EIS_Arctic", "LAU_Arctic",
           "GIS_GreenlandSea", "LAU_LabradorSea", "LAU_StLawrence", "LAU_GulfOfMexico"]
nisa = m[(m.site == "La Vallina (NISA)") & (m.pathway == "mixed")].copy()
for k in KEYS:
    nisa[k] = nisa[k].clip(lower=0)
mx = {k: nisa[k].max() for k in KEYS}          # each spoke scaled to the site maximum
val = lambda f, r, k: float(nisa[(nisa.forcing == f) & (nisa.region == r)][k].iloc[0])

# Palette and type (same family as the explorer and the report)
BG, PANEL, INK, BODY, MUTED, RULE, ACC = "#F3F6F8", "#FBFCFD", "#16212A", "#4A5A67", "#74838F", "#D7DFE5", "#1C5CAB"
GL, IC = "#184F95", "#C26A00"                  # GLAC-1D, ICE-6G (as in the report)
DK_BG, DK_INK, DK_BODY, DK_MUTED, DK_GL, DK_IC = "#16212A", "#E8EEF2", "#BFD0DB", "#8FA0AC", "#86B6EF", "#F0A04B"
HEAD = "font-family:'Barlow Condensed', 'Arial Narrow', Arial, sans-serif"
TEXT = "'IBM Plex Sans', Arial, sans-serif"


def f(x):
    return f"{x:.1f}"


def ang(k):
    return math.pi / 2 - 2 * math.pi * k / 5


def pt(cx, cy, r, k):
    return cx + r * math.cos(ang(k)), cy - r * math.sin(ang(k))


def spider_svg_parts(cx, cy, R, polys, spoke_w=1.5):
    """Rings, spokes and polygons for one spider. polys = [(radii[5], colour)]."""
    out = []
    for rr in (0.5, 1.0):
        pts = " ".join(f"{f(x)},{f(y)}" for x, y in (pt(cx, cy, R * rr, k) for k in range(5)))
        out.append(f'<polygon points="{pts}" fill="none" stroke="{RULE}" stroke-width="{spoke_w}"/>')
    for k in range(5):
        x, y = pt(cx, cy, R, k)
        out.append(f'<line x1="{f(cx)}" y1="{f(cy)}" x2="{f(x)}" y2="{f(y)}" stroke="{RULE}" stroke-width="{spoke_w}"/>')
    for radii, col in polys:
        pts = " ".join(f"{f(x)},{f(y)}" for x, y in (pt(cx, cy, R * max(r, 0.015), k) for k, r in enumerate(radii)))
        out.append(f'<polygon points="{pts}" fill="{col}" fill-opacity="0.13" stroke="{col}" stroke-width="{spoke_w + 1}" stroke-linejoin="round"/>')
    return out


def radii(forcing, region):
    return [val(forcing, region, k) / mx[k] for k in KEYS]


# ------------------------------------------------------------------ slide 1: how to read it
W1, H1, C1 = 700, 540, (350, 285)
R1 = 210
ex = "LAU_LabradorSea"
svg1 = (f'<svg aria-label="Example spider with five spokes: Share at the top, then clockwise Dominant, SD, Var. share and SD Δ. '
        f'Two polygons for {ex} at NISA: GLAC-1D reaches far out on every spoke, ICE-6G stays near the centre." '
        f'width="{W1}" height="{H1}" viewBox="0 0 {W1} {H1}" xmlns="http://www.w3.org/2000/svg">'
        + "".join(spider_svg_parts(*C1, R1, [(radii("ICE-6G", ex), IC), (radii("GLAC-1D", ex), GL)], spoke_w=2))
        + "</svg>")
tips = {k: pt(*C1, R1, i) for i, k in enumerate(KEYS)}
lab = lambda text, left, top, width, align: (
    f'<p style="position:absolute;left:{left}px;top:{top}px;width:{width}px;font-size:32px;{HEAD};'
    f'font-weight:600;color:{INK};text-align:{align};line-height:1">{text}</p>')
spoke_labels = "".join([
    lab("Share", int(tips["share"][0] - 100), int(tips["share"][1] - 48), 200, "center"),
    lab("Dominant", int(tips["dom"][0] + 14), int(tips["dom"][1] - 16), 140, "left"),
    lab("SD", int(tips["sd"][0] - 4), int(tips["sd"][1] + 12), 90, "left"),
    lab("Var. share", int(tips["vshare"][0] - 150), int(tips["vshare"][1] + 12), 150, "right"),
    lab("SD Δ", int(tips["sd_diff"][0] - 154), int(tips["sd_diff"][1] - 16), 140, "right"),
])
rows = [
    ("Share", "How much of the signal comes from here?", "Mean contribution ÷ sum over all 9 regions"),
    ("Dominant", "How often is it the biggest contributor?", "Share of 500-yr steps in which it is the largest"),
    ("SD", "How much does it swing over time?", "Standard deviation over time, in ‰"),
    ("Var. share", "Which part of the total's swings is it?", "Cov(region, total) ÷ Var(total); sums to 100 %"),
    ("SD Δ", "How abrupt are its changes?", "Standard deviation of the 500-yr changes, in ‰"),
]
row_html = "".join(
    f'<div style="display:flex;gap:24px;align-items:start">'
    f'<h3 style="width:190px;font-size:44px;{HEAD};font-weight:600;color:{ACC};line-height:1">{n}</h3>'
    f'<div style="display:flex;flex-direction:column;gap:4px;flex:1">'
    f'<p style="font-size:28px;color:{INK};line-height:1.3">{q}</p>'
    f'<p style="font-size:24px;color:{BODY};line-height:1.3">{h}</p></div></div>'
    for n, q, h in rows)
slide1 = f"""<section id="explain" data-transition="fade" style="background:{BG};color:{INK};font-family:{TEXT};padding:128px;display:flex;flex-direction:column;gap:40px">
<div style="display:flex;flex-direction:column;gap:8px">
<p style="font-size:24px;color:{ACC};text-transform:uppercase;letter-spacing:2px;font-weight:600">Reading the spider figures</p>
<h2 style="font-size:72px;{HEAD};font-weight:600;line-height:1.1;color:{INK}">One spider, five questions about a source region</h2>
</div>
<div style="display:flex;gap:64px;align-items:start;flex:1">
<div style="position:relative;width:700px;height:640px">
{svg1.replace('<svg ', '<svg style="position:absolute;left:0px;top:0px;width:700px;height:540px" ', 1)}
{spoke_labels}
<p style="position:absolute;left:0px;top:560px;width:700px;font-size:24px;color:{BODY};line-height:1.35;text-align:center">Example: {ex} at NISA · <span style="color:{GL}">GLAC-1D</span> · <span style="color:{IC}">ICE-6G</span></p>
</div>
<div style="display:flex;flex-direction:column;gap:26px;flex:1">
{row_html}
<p style="font-size:24px;color:{MUTED};line-height:1.35">Computed on a 500-yr grid, 21.5–10 ka. Each spoke is scaled to the site's largest value, so the outer ring is the site maximum and the inner ring half of it.</p>
</div>
</div>
<aside>How to read the spider figures. Each spider describes one source region at one site. The five spokes each ask one question. Share: how much of the site's meltwater signal comes from this region, its mean contribution divided by the sum over all nine regions. Dominant: in how many 500-year steps it is the single largest contributor. SD: how much its contribution swings over time. Variance share: how much of the ups and downs of the total it causes, the covariance with the total divided by the variance of the total; these add up to 100 percent. SD Delta: how abrupt it is, the standard deviation of the step-to-step changes. Every spoke is scaled to the largest value at the site, so a polygon that reaches the outer ring is the site's leader on that question. The example is the Labrador Sea region at NISA: under GLAC-1D it is large on every spoke, under ICE-6G it stays small.</aside>
</section>"""

# ------------------------------------------------------------------ slide 2: nine NISA spiders
W2, H2 = 1000, 650
cw, ch, title_h = W2 / 3, H2 / 3, 36
R2 = min(cw, ch - title_h) / 2 - 12
parts, titles = [], []
for i, reg in enumerate(REGIONS):
    col, row = i % 3, i // 3
    x0, y0 = col * cw, row * ch
    cx, cy = x0 + cw / 2, y0 + title_h + (ch - title_h) / 2
    parts += spider_svg_parts(cx, cy, R2, [(radii("ICE-6G", reg), IC), (radii("GLAC-1D", reg), GL)])
    titles.append(f'<p style="position:absolute;left:{int(x0)}px;top:{int(y0)}px;width:{int(cw)}px;font-size:24px;'
                  f'font-weight:600;color:{INK};text-align:center;line-height:1.2">{reg}</p>')
svg2 = (f'<svg aria-label="Nine spiders for La Vallina (NISA), one per source region, GLAC-1D in blue and ICE-6G in orange. '
        f'Under ICE-6G, LAU_GulfOfMexico fills its spider on every spoke; under GLAC-1D the large polygons are '
        f'LAU_LabradorSea, EIS_NorwegianSea and EIS_BayOfBiscay." style="position:absolute;left:0px;top:0px;width:{W2}px;height:{H2}px" '
        f'width="{W2}" height="{H2}" viewBox="0 0 {W2} {H2}" xmlns="http://www.w3.org/2000/svg">' + "".join(parts) + "</svg>")
KW, KH, KC, KR = 360, 250, (180, 135), 88
key_svg = (f'<svg aria-label="Key: Share at the top, then clockwise Dominant, SD, Var. share and SD Δ." '
           f'style="position:absolute;left:0px;top:0px;width:{KW}px;height:{KH}px" width="{KW}" height="{KH}" viewBox="0 0 {KW} {KH}" '
           f'xmlns="http://www.w3.org/2000/svg">' + "".join(spider_svg_parts(*KC, KR, [])) + "</svg>")
kt = {k: pt(*KC, KR, i) for i, k in enumerate(KEYS)}
klab = lambda text, left, top, width, align: (
    f'<p style="position:absolute;left:{left}px;top:{top}px;width:{width}px;font-size:24px;color:{BODY};'
    f'text-align:{align};line-height:1">{text}</p>')
key_labels = "".join([
    klab("Share", int(kt["share"][0] - 60), int(kt["share"][1] - 32), 120, "center"),
    klab("Dominant", int(kt["dom"][0] + 8), int(kt["dom"][1] - 12), 110, "left"),
    klab("SD", int(kt["sd"][0]), int(kt["sd"][1] + 8), 60, "left"),
    klab("Var. share", int(kt["vshare"][0] - 124), int(kt["vshare"][1] + 8), 124, "right"),
    klab("SD Δ", int(kt["sd_diff"][0] - 88), int(kt["sd_diff"][1] - 12), 80, "right"),
])
pc = lambda f_, r, k: f"{100 * val(f_, r, k):.0f} %"
slide2 = f"""<section id="nisa" data-transition="fade" style="background:{PANEL};color:{INK};font-family:{TEXT};padding:128px;display:flex;flex-direction:column;gap:32px">
<div style="display:flex;flex-direction:column;gap:8px">
<p style="font-size:24px;color:{ACC};text-transform:uppercase;letter-spacing:2px;font-weight:600">La Vallina (NISA) · mixed AMOC pathway · 21.5–10 ka</p>
<h2 style="font-size:72px;{HEAD};font-weight:600;line-height:1.1;color:{INK}">Two reconstructions, two main sources</h2>
</div>
<div style="display:flex;gap:64px;align-items:start;flex:1">
<div style="position:relative;width:{W2}px;height:{H2}px">
{svg2}
{''.join(titles)}
</div>
<div style="display:flex;flex-direction:column;gap:24px;flex:1">
<div style="position:relative;width:{KW}px;height:{KH}px">
{key_svg}
{key_labels}
</div>
<p style="font-size:28px;line-height:1.3;font-weight:600"><span style="color:{GL}">GLAC-1D</span>  ·  <span style="color:{IC}">ICE-6G</span></p>
<p style="font-size:24px;color:{BODY};line-height:1.4"><b>GLAC-1D</b> spreads the signal: LAU_LabradorSea leads ({pc('GLAC-1D','LAU_LabradorSea','share')}), with EIS_NorwegianSea ({pc('GLAC-1D','EIS_NorwegianSea','share')}) and EIS_BayOfBiscay ({pc('GLAC-1D','EIS_BayOfBiscay','share')}) close behind.</p>
<p style="font-size:24px;color:{BODY};line-height:1.4"><b>ICE-6G</b> concentrates it: LAU_GulfOfMexico gives {pc('ICE-6G','LAU_GulfOfMexico','share')} of the signal, leads in {pc('ICE-6G','LAU_GulfOfMexico','dom')} of the steps and causes {pc('ICE-6G','LAU_GulfOfMexico','vshare')} of the variance.</p>
</div>
</div>
<aside>The nine spiders for NISA, one per source region, both reconstructions, mixed AMOC pathway. Read them against the key on the right: a polygon that reaches the outer ring leads the site on that question. Under GLAC-1D the signal is spread: the Labrador Sea region has the largest share, about a quarter, and the Norwegian Sea and Bay of Biscay regions follow closely; the Norwegian Sea is the most abrupt. Under ICE-6G one region dominates: the Gulf of Mexico fills its spider on every spoke, with about a third of the signal, more than half of the time steps and more than half of the variance. So the answer to which ice sheet NISA records depends on the reconstruction.</aside>
</section>"""

# ------------------------------------------------------------------ slide 3: takeaways
# Site summary, as site_summary() in R/site_metrics.R.
def leader(g, col):
    return g.loc[g[col].idxmax()]
summ = (m.groupby(["forcing", "pathway", "site"])
         .apply(lambda g: pd.Series({"important": leader(g, "share").region, "share": g.share.max(),
                                     "abrupt": leader(g, "sd_diff").region}), include_groups=False)
         .reset_index())
mix = summ[summ.pathway == "mixed"]
piv = mix.pivot(index="site", columns="forcing", values="important")
n_sites = len(piv)
n_agree = int((piv["GLAC-1D"] == piv["ICE-6G"]).sum())
n_cases = len(mix)
n_fen = int((mix.abrupt == "EIS_NorwegianSea").sum())
n_robust = int(summ.groupby(["forcing", "site"]).important.nunique().eq(1).sum())
top = lambda f_, s: mix[(mix.forcing == f_) & (mix.site == s)].iloc[0]
g_n, i_n = top("GLAC-1D", "La Vallina (NISA)"), top("ICE-6G", "La Vallina (NISA)")
rng = lambda s: sorted(mix[mix.site == s].share)
md, cw_ = rng("MD02-2552"), rng("Cave Without a Name")
P = lambda x: f"{100 * x:.0f}"
take = [
    (f"{n_agree} / {n_sites}", "sites where GLAC-1D and ICE-6G agree on the most important source region."),
    (f'<span style="color:{DK_GL}">{P(g_n.share)} %</span> vs <span style="color:{DK_IC}">{P(i_n.share)} %</span>',
     f"NISA: {g_n.important} leads under GLAC-1D, {i_n.important} under ICE-6G."),
    (f"{P(md[0])}–{P(md[1])} %",
     f"LAU_GulfOfMexico's share at the Gulf core MD02-2552; {P(cw_[0])}–{P(cw_[1])} % at Cave Without a Name (Texas)."),
    (f"{n_fen} / {n_cases}", "site and reconstruction cases in which EIS_NorwegianSea is the most abrupt source."),
    (f"{n_robust} / {n_cases}", "cases in which the main source is the same in all three AMOC pathways, so the pathway matters."),
]
take_html = "".join(
    f'<div style="display:flex;gap:48px;align-items:center">'
    f'<p style="width:420px;font-size:72px;{HEAD};font-weight:600;color:{DK_GL};line-height:1.05;white-space:nowrap">{n}</p>'
    f'<p style="flex:1;font-size:32px;color:{DK_BODY};line-height:1.35">{t}</p></div>'
    for n, t in take)
slide3 = f"""<section id="takeaways" data-transition="fade" style="background:{DK_BG};color:{DK_INK};font-family:{TEXT};padding:128px 128px 160px;display:flex;flex-direction:column;gap:40px">
<div style="display:flex;flex-direction:column;gap:8px">
<p style="font-size:24px;color:{DK_GL};text-transform:uppercase;letter-spacing:2px;font-weight:600">Takeaways · nine proxy sites</p>
<h2 style="font-size:72px;{HEAD};font-weight:600;line-height:1.1;color:{DK_INK}">Where the meltwater arrives, site by site</h2>
</div>
<div style="display:flex;flex-direction:column;gap:22px">
{take_html}
</div>
<p style="position:absolute;left:128px;bottom:64px;width:1664px;font-size:24px;color:{DK_MUTED}">Forward model v2 · HadCM3 pulse responses · mixed AMOC pathway unless stated · 21.5–10 ka, 500-yr grid</p>
<aside>Five takeaways from the nine sites. First, the two ice-sheet reconstructions agree on the most important source region at only five of the nine sites. Second, at NISA they disagree: GLAC-1D points to the Labrador Sea region with about a quarter of the signal, ICE-6G to the Gulf of Mexico with about a third. Third, the two sites near the Gulf of Mexico are dominated by the Gulf of Mexico region in both reconstructions. Fourth, the Norwegian Sea region is the most abrupt contributor in eleven of eighteen site and reconstruction cases. Fifth, the main source stays the same across the cold, zonal and mixed AMOC pathways in only nine of eighteen cases, so the choice of pathway matters as much as the choice of reconstruction.</aside>
</section>"""

for sid, html in (("explain", slide1), ("nisa", slide2), ("takeaways", slide3)):
    (ROOT / f"project/slides/{sid}.html").write_text(html, encoding="utf-8")

deck = {
    "v": 4,
    "createdOnFiles": {"v": 1, "at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
    "lists": "css",
    "title": "NISA Spider Slides",
    "order": ["explain", "nisa", "takeaways"],
    "cover": "explain",
    "sections": {"s1": {"description": "How to read the spider, NISA's nine spiders, and the takeaways across the nine sites",
                        "start": "explain"}},
    "faces": {
        "barlow-condensed": {"family": "Barlow Condensed",
                             "href": "https://fonts.googleapis.com/css2?family=Barlow+Condensed:wght@500;600&display=swap"},
        "ibm-plex-sans": {"family": "IBM Plex Sans",
                          "href": "https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&display=swap"},
    },
    "designSystems": [],
}
if args.index:
    (ROOT / "project/deck.json").write_text(json.dumps(deck, ensure_ascii=False, indent=1), encoding="utf-8")
for p in sorted(ROOT.rglob("project/**/*.*")):
    print(p.relative_to(ROOT), f"{p.stat().st_size / 1024:.1f} KB")
print("takeaways:", [t[0] for t in take])
