# CLAUDE.md

Research repo that answers the questions in `questions.md` about the NISA
stalagmite **Glas** (La Vallina cave, NW Iberia) and deglacial North Atlantic
meltwater. Each answer is a Quarto manuscript that renders to DOCX. The owner
is Laura Endres.

## Status of the questions

| # | Question | Output |
|---|---|---|
| 1 | Age range of Glas | `glas_age_range.qmd` / `.docx` |
| 2 | Most important and most variable meltwater source at NISA; alignment of abrupt Glas δ¹⁸O shifts with discharge | `nisa_meltwater_sources.qmd` / `.docx`, plus the explorer in `app/` |
| 3 | Age range of the Lake Gerzensee record | Open. No Gerzensee data is in `data/` yet |

## Layout

- `data/`: SISAL-format database excerpt (schema in `data/schema.dbml`), the two regional discharge series, and `meltmodel_site_weights.csv` (with provenance in `meltmodel_site_weights.md`).
- `*.qmd` and `*.docx` at the root: one manuscript per question. The DOCX is committed alongside its source.
- `app/`: Meltwater Discharge Explorer, an interactive map and time series of the discharge. See `app/README.md`. `regions_hadcm3.geojson` holds the exact HadCM3 input regions (made by `scripts/make_region_geojson.py`).
- **Convolution pipeline (v2)**, laid out like the meltmodel repo:
  - `myconfig/`: dyes and colours, sites, forcings, pathways and paths.
  - `mymodules/`: the forcing loader, site extraction, kernels and the convolution.
  - `scripts/`: `export_pulse_kernels.py` (runs on the Leeds server), `make_uptake_weights.py`, `make_region_geojson.py` and `run_convolution.py`.
  - Outputs go to `outputs/convolution/`.
  - New forcings are one entry in `myconfig/FORCINGS.py`, and new sites one entry in `myconfig/SITES.py`.
- `data/kernels/`: the HadCM3 pulse-kernel bundle. It is **untracked** (see its README for how to make it). `data/meltmodel/uptake_weights.nc` holds the land-site moisture-uptake fields regridded from the trajectories.
- `plans/`: plans written before large tasks, dated.
- `prompts/`: verbatim prompt archive written by the commit skill (`YYYY-MM-DD-NNN-slug.md`).
- `Impulse_Response_Revisions_v1.pdf`: unpublished forward-model draft (Endres et al., in prep.). It is deliberately **untracked**, so never commit it. The same goes for Word lock files (`~$*.docx`).

## Commands

Quarto is not on PATH. Use the copy bundled with RStudio:

```sh
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render nisa_meltwater_sources.qmd   # about 3.5 min (Bchron and surrogates)
python3 app/build_app.py                        # regenerate app/discharge_explorer.html
python3 app/build_app.py --publish-out <path>   # also write the artifact copy (no <html>/<head> wrapper)
```

After a render, read the text back to check every inline value and claim:

```sh
pandoc nisa_meltwater_sources.docx -t plain --wrap=none
```

Convolution pipeline (needs the conda `base` Python with xarray and scipy; about 30 s):

```sh
python3 scripts/make_uptake_weights.py      # once: trajectories -> data/meltmodel/uptake_weights.nc
python3 scripts/make_region_geojson.py      # once: HadCM3 input cells -> app/regions_hadcm3.geojson
python3 scripts/run_convolution.py          # uses data/kernels/ if present, else the placeholder
```

Without the kernel bundle, `run_convolution.py` falls back to a **placeholder**: the equilibrium field spread over time with an invented 100-yr response time. Its outputs go to `outputs/convolution/placeholder/`, which is gitignored. Never report, commit or publish them as results.

R 4.6 is installed with readr, dplyr, tidyr, purrr, ggplot2, patchwork, strucchange, Bchron, zoo, boot and knitr. **Not** installed: changepoint, flextable and kableExtra. Tables use `knitr::kable()`. The app build uses only the Python standard library.

## Data conventions

- **Glas** is `entity_id` 903. It has 532 samples joined through `sample.csv`, `original_chronology.csv` (published Bchron ages, 11,939–23,819 yr BP) and `d18O.csv`, plus 25 U–Th dates in `dating.csv`. Depths are in mm from the top. `sisal_chronology.csv` has no Glas rows.
- **Discharge files.** Time runs forward in model years: negate `time` (GLAC-1D, 100-yr steps) or `t_adj` (ICE-6G, 500-yr steps) to get yr BP. As a sanity check, total discharge peaks at Meltwater Pulse 1A: 0.260 Sv at 14.4 ka and 0.243 Sv at 14.0 ka.
- **Region codes** in file order are `Med, Bri, Fen, EurArc, AmeArc, GIS, NLau, SLau, GulofMex`. They map one to one onto Endres et al. (2026b) names: `EIS_MedSea, EIS_BayOfBiscay, EIS_NorwegianSea, EIS_Arctic, LAU_Arctic, GIS_GreenlandSea, LAU_LabradorSea, LAU_StLawrence, LAU_GulfOfMexico`. `LAU_Arctic` drains to the Beaufort Sea.
- **Anomaly columns** (`<code> d18O (-30/-35/-40)`) hold the *source-region* δ¹⁸O anomaly, which is linear in discharge. They are the forcing that the forward-model draft convolves, not the anomaly at NISA. Use the −35 ‰ end-member; the others only rescale everything by 6/7 or 8/7. Endres et al. (2026b) use a saturating mixing fraction instead, so their scenario values are smaller for large pulses.
- **Site anomalies (v1).** Multiply source anomalies by `meltmodel_site_weights.csv` (site × region × AMOC mode, cold or zonal). This is an *equilibrium* approximation: it has no transit delay and overstates short pulses, so always say so.
- **Site anomalies (v2).** Convolve the decadal forcing with the pulse kernels: `A(t) = Σ_k 10·F(t−k)·h(k)`, as in the 2025.05 notebook (`old script/`).
  - The pulse runs are xpran (cold), xprao (zonal) and xpujc (merid). Each starts after its 10-yr pulse, so the export prepends the parent run's first 10 yr (xpraj, xprak, xpram). xpram was confirmed from the time axes: it starts in 4511 and xpujc in 4521.
  - The mixed pathway switches modes on the notebook's schedule. Its first and last modes extend to the ends of the forcing, which gives a spin-up.
  - Land sites use the trajectory uptake fields. The merid segment borrows the zonal uptake, because there is no merid trajectory run.
  - `scaling_check.csv` compares the 500-yr step response with the equilibrium weights, and should be about 1.
- **Published abrupt Glas events** (Endres et al. 2026a, Table 1) are hard-coded in both `nisa_meltwater_sources.qmd` and `app/build_app.py`. Define events by sample **depth**, so that they are re-dated in every age-model ensemble member.

## Analysis standards (question 2)

The design was reviewed from Kira Rehfeld's perspective: irregular, age-uncertain series. Keep these standards when extending it:

- Propagate Glas age uncertainty with a Bchron ensemble rerun on the U–Th dates. Report the median and the 5–95 % range, not point estimates.
- Use Gaussian-kernel correlation (Rehfeld et al. 2011) on the irregular samples rather than interpolating. Separate the levels, which carry the millennial trend, from the high-passed centennial part.
- Test against red noise: AR(1) surrogates of Glas, phase-randomised surrogates of the discharge, and a family-wise maximum over lags and regions.
- Frame results as "consistent with" or "detectable at the ~500-yr scale". Centennial alignment could not be shown.

## Manuscript conventions

- The YAML follows `glas_age_range.qmd`: `format: docx`, `toc: false`, `fig-dpi: 300`, `echo: false`, `warning: false`, `message: false`, `date: today`, author Laura Endres.
- Every number in the prose is inline R, so re-render and re-check the prose whenever an analysis changes. Qualitative claims ("negligible", "strongly", "in both modes") need checking too, not just the numbers.
- Chunks that print progress bars, such as `Bchronology`, need `#| results: hide`.
- In figure labels write δ¹⁸O as plotmath, `expression("Glas" ~ delta^18*O ~ "(‰)")`. The Unicode superscript "⁸" renders as a box in the PNG device.
- Never put prose or TODO lines inside an R chunk, because the render breaks.
- Write in short, plain sentences, and give a plain-language explanation next to each statistical method.
- Verify references before citing them: use Zotero for Laura's papers and web search otherwise. The QSR 2019 paper "Correlating paleoclimate time series" is by Franke and Donner, not Rehfeld's group.

## Colours

The region colours match the paper figures (matplotlib tab10), with lightness tuned to pass colour-blindness checks. They are defined in **two places that must stay in sync**: `reg_col` in `nisa_meltwater_sources.qmd` and the `--r-*` tokens in `app/template.html`.

| Region | Light | Dark (app only) |
|---|---|---|
| EIS_MedSea | `#1f77b4` | `#1f77b4` |
| EIS_BayOfBiscay | `#f1780b` | `#db6b00` |
| EIS_NorwegianSea | `#2ca02c` | `#2ca02c` |
| EIS_Arctic | `#d62728` | `#d62728` |
| LAU_Arctic | `#9467bd` | `#9467bd` |
| GIS_GreenlandSea | `#9b4c3c` | `#a15142` |
| LAU_LabradorSea | `#de72bd` | `#cd63ad` |
| LAU_StLawrence | `#a5a608` | `#959600` |
| LAU_GulfOfMexico | `#00b3c3` | `#03a2b1` |

Stacked charts use this bottom-to-top order, which keeps neighbouring colours distinguishable: Bri, NLau, SLau, GoM, EurArc, Med, Fen, AmeArc, GIS. The untuned paper colours fail the check, because Bri and Fen look alike to colour-blind readers. If the palette changes, re-run the dataviz skill's `validate_palette.js` in both themes, and show Laura the palette before implementing it.

## Explorer app

- `discharge_explorer.html` is generated, so edit `template.html` or `build_app.py` and rebuild.
- `regions_approx.geojson` holds hand-digitised, approximate outlines from Endres et al. (2026b, Fig. 2a). Exact HadCM3 masks with the same `code` property can replace it without code changes.
- The shared page is https://claude.ai/artifact/HMQMC2ACsbpvvDQNyRgcPG, which is private until shared. Republish the `--publish-out` copy to the same URL to keep the link.
- The browser tools cannot open `file://`. To test locally, serve the folder with `python3 -m http.server 8765 --bind 127.0.0.1` from `app/` and stop the server afterwards.

## Git workflow

- Work on `dev`, commit with the `ghe-skills:commit` skill, and open PRs from `dev` into `main` with `ghe-skills:open-pr`. PRs are merged with a merge commit, and `dev` is kept.
- Commits follow Conventional Commits, archive prompts under `prompts/`, and end with `Assisted-by: Claude <model-id>` rather than Co-Authored-By.
- The commit and PR skills do not push. Push only when Laura asks.
- PR bodies carry no commit trailers and no emoji. Put `Closes #N` in the body so issues close on merge.

## Related sources

- **Endres et al. (2026a)**, Clim. Past 22, 797–824: the Glas record and Table 1 events. Zotero key `4E5CEBT8`.
- **Endres et al. (2026b)**, "Tracing meltwater…", the "meltmodel" paper: input regions, Tables 1–2, proxy sites. Zotero key `KWBMYG69`. Figure scripts and cached site data are in `~/Documents/gitrepos/gh-lrndrs/gh-lrndrs-dyetracer_palaeo_figurescripts`, archived at Zenodo 10.5281/zenodo.21828703.
- **Paper registry.** Laura's paper configs (Overleaf paths, figures) are in `~/Documents/Obsidian_Matriarchy/Research_Notes/Agentic/Pacco/papers.md`.
