# Plan: fresh answer to question 2 as a Quarto → DOCX manuscript (age-uncertainty aware)

## Context

`questions.md`, question 2: *At the NISA location, which meltwater source region is the most important and which the most variable contributor to the δ¹⁸O anomaly over time? Is there alignment between abrupt δ¹⁸O changes in Glas and the discharge, and how can that be quantified?*

A first answer exists (`nisa_meltwater_sources.qmd`, merged in PR #2). The user chose **a fresh, independent analysis that replaces it** (same filenames; old version stays in git history), and asked for two additions after the first plan draft:

1. Revise the design from the perspective of **Kira Rehfeld as reviewer** (Tübingen; palaeoclimate variability, irregular and age-uncertain time series, proxy–model comparison, SISAL).
2. **Propagate Glas age uncertainty**: Glas δ¹⁸O can shift within its chronological uncertainty, which may change every correlation and alignment statistic.

Nothing has been run or written yet. Toolchain verified read-only: R 4.6.1; readr/dplyr/tidyr/purrr/ggplot2/patchwork/strucchange/zoo/boot/matrixStats/knitr present; **Bchron 4.7.8 and rbacon 4.0.0 installed**; **Quarto 1.9.38 only at `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`** (not on PATH). Missing (do not use): changepoint, nest, PaleoSpec, flextable, kableExtra.

## Reviewer perspective (Rehfeld) → what changes in the design

Themes of her programme, each mapped to a concrete change. Attribution kept to work she authored or co-authored: Rehfeld et al. 2011 (NPG, correlation of irregularly sampled series; kernel estimators beat interpolation), Rehfeld & Kurths 2014 (CP, similarity estimators for irregular *and age-uncertain* series via age ensembles), [wrongly attributed in the first version of this plan: the QSR 2019 paper "Correlating paleoclimate time series: sources of uncertainty and potential pitfalls" is by Franke & Donner, not by Rehfeld's group; it is still used for the multiple-testing point], Rehfeld et al. 2018 (Nature, timescale-dependent variability), Bühler et al. 2021 (CP, speleothem δ¹⁸O vs iHadCM3: low proxy signal-to-noise), Comas-Bru et al. 2020 (ESSD, SISAL v2 age-model ensembles).

| Likely reviewer objection | Change |
|---|---|
| Point age model treated as truth; "alignment within ±300 yr" is meaningless if Glas ages carry ±150–600 yr uncertainty | Age-model **ensemble** for Glas (Sect. B0). Every alignment statistic is a distribution over ensemble members, reported as median and 5–95 % range. Events defined by **depth**, not age, so they are re-dated in each member |
| Discharge chronologies also uncertain (GLAC-1D/ICE-6G tied to sea-level and ¹⁴C-dated margins) but not in the data | Systematic **chronology-offset scan**: shift the discharge series by −1000…+1000 yr (100-yr steps) and show how each statistic responds; state which offsets the Glas ensemble can and cannot absorb. Acknowledge that reconstruction uncertainty is not quantified in the repo |
| Interpolation/binning of an irregular series biases correlation and persistence estimates | Report Glas sampling resolution (median step, range) against the 100-yr grid; use **Gaussian-kernel cross-correlation** (Rehfeld et al. 2011; implemented inline, ~15 lines) on the raw irregular ages as the primary estimator, binned Pearson r only as a cross-check |
| Shared deglacial trend (and the ~−1 ‰ global-ocean δ¹⁸O drift) inflates level correlations | Analyse **first differences and high-pass (< 2 kyr) components**; report level correlations only with the caveat. Timescale-dependent results (centennial vs millennial) rather than one r |
| Significance from AR(1) surrogates of one series only, and multiple testing across 9 regions × 41 lags | Surrogates for **both** series (AR(1) fitted to Glas; phase-randomised surrogates of the discharge preserving its spectrum); **max-statistic null** over regions and lags for the family-wise test; report effective sample size |
| Comparing SD across GLAC-1D (100 yr) and ICE-6G (500 yr) confounds variability with resolution | Compute all variability metrics on a **common 500-yr grid** for the cross-reconstruction comparison; give 100-yr GLAC-1D values separately |
| Proxy–model comparison in absolute units; speleothem δ¹⁸O has low signal-to-noise and carries temperature and source effects | Standardise (z-score within the Glas interval) before comparison; add a **detectability check**: given the GLAC-1D forcing plus white/red noise at the SNR implied by the residuals, what correlation would be expected? Frame the result as "consistent with / not distinguishable from" rather than "aligned" |
| Ad hoc event lists | Events from an objective detector (mean-shift breakpoints, `strucchange`) **and** the published Table 1 events, both re-dated per ensemble member; sensitivity table over detector and tolerance |
| Reproducibility | Fixed seed, all code in the qmd, ensemble sizes stated, runtime kept ≤ ~10 min |

## Data facts established (read-only)

| Item | Fact |
|---|---|
| `data/Discharge_Glac1d_regional-withd18O.csv` | 261 rows, `time` = −26000…0 (100-yr steps). Discharge cols `Med, Bri, Fen, EurArc, AmeArc, GIS, NLau, SLau, GulofMex` (Sv); anomaly cols `<code> d18O (-35.0/-40.0/-30.0)` |
| `data/Discharge_ice6g_regional-withd18O.csv` | 61 rows, `t_adj` = −28000…+2000 (500-yr steps). Discharge cols prefixed `discharge_`; same anomaly cols; ignore `t, depth, discharge_lsm, discharge_Sv_residual` |
| Time convention | negate `time` / `t_adj` → years BP (check: MWP-1A peak ≈ 14.6 ka in both) |
| Region order / names | Med, Bri, Fen, EurArc, AmeArc, GIS, NLau, SLau, GoM = input regions of Endres et al. (2026b) in the same order (EIS_MedSea, EIS_BayOfBiscay, EIS_NorwegianSea, EIS_Arctic, LAU_Arctic, GIS_GreenlandSea, LAU_LabradorSea, LAU_StLawrence, LAU_GulfOfMexico) |
| Anomaly meaning | Draft `Impulse_Response_Revisions_v1.pdf` Sect. 5.2.3/5.3.5: anomaly = discharge scaled to ice end-member ÷ region volume = **source-region forcing** of the Green's-function model. The convolved NISA series (Fig. 5.13 C/E) is not in the repo → "contribution at NISA" is evaluated on the forcing side and stated as such |
| Glas | entity_id 903; 532 samples, depths 3.9–39.6 (units as in database); `original_chronology.csv` Bchron ages 11,939–23,819 yr BP with asymmetric uncertainties (median ≈ 148 / 130 yr, up to ≈ 1,200 yr at the young end); `d18O.csv` 532 values |
| Glas dates | `dating.csv`: 25 U-Th dates, depths 3.3–38.3 (same units as samples, bracket them), `corr_age` 10,956–23,524 yr BP, symmetric 2σ-style uncertainties 47–97 yr, all `date_used = yes`; one minor stratigraphic reversal (6.0 vs 7.7) |
| Reference events | Endres et al. (2026a, CP 22, 797–824; `reference.csv` ref_id 563) Table 1 abrupt Glas events (23.4kS, 17.8kS, 16.4kS, 16.1kS, 15.3kS, 14.7kS). Ages/CIs to be **verified from the paper via Zotero (read-only)** before hard-coding; then converted to **depths** using the original chronology |

## Deliverables

1. `nisa_meltwater_sources.qmd` — overwritten (YAML as `glas_age_range.qmd`: `format: docx`, `toc: false`, `fig-dpi: 300`, `echo: false`, `warning/message: false`, `date: today`, author Laura Endres).
2. `nisa_meltwater_sources.docx` — rendered.
3. No other files, no commit unless asked (then via `ghe-skills:commit`, which handles `prompts/`). If Bchron runtime forces caching, the cache dir would need a `.gitignore` line; flag this to the user rather than adding silently.

## Analysis design (R chunks in the qmd; `set.seed(20260929)`)

### 0. Setup and sanity checks
- Tidy both discharge files to long `(recon, age, code, region, ice_sheet, Q, anom)` for −35 ‰; note −30/−40 are exact rescalings.
- Checks printed in text: `anom/Q` constant per region; MWP-1A peak age; Glas interval and n; Glas median sampling step.
- Common grids: 100 yr (GLAC-1D) and 500 yr (both).

### A. Importance and variability of source regions (forcing side)
Over the Glas interval and the full series; on the 500-yr common grid for cross-reconstruction comparison, 100-yr for GLAC-1D detail.
- **Importance:** time-integrated anomaly and share of the all-region total; **dominance fraction** (share of time steps in which the region is the largest contributor).
- **Variability:** SD; **variance-share decomposition** Cov(anomᵢ, total)/Var(total); SD of first differences (abruptness); all on the common grid.
- Ice-sheet aggregates. Stacked-area figure; ranking tables (kable).

### B0. Glas age-model ensemble (new; feeds everything in B)
- Primary: rerun **Bchron** (`Bchronology`, `calCurves = "normal"`, default outlier handling) on the 25 `dating.csv` dates (depth, `corr_age`, uncertainties), `predictPositions` = the 532 sample depths; keep N = 1,000 realisations. Validate: ensemble median vs `interp_age` (report max deviation), ensemble 2.5–97.5 % vs `interp_age_uncert_neg/pos`. If the match is poor, report it and use the fallback as primary.
- Fallback / cross-check: coherent perturbation of the published chronology, age = interp_age + z·(uncert_pos if z>0 else uncert_neg) with one z per realisation plus small independent jitter, enforced monotonic.
- Also a **deterministic offset scan** Δ ∈ [−1000, +1000] yr applied to the discharge series, with the ensemble telling which Δ are inside Glas's uncertainty.
- Figure: age–depth ensemble envelope with dates; ensemble spread through time.

### B. Alignment of abrupt Glas δ¹⁸O changes with discharge
Glas events (defined by **depth**, dated per ensemble member):
1. Published Table 1 events (verified), located to depth.
2. Objective mean-shift breakpoints (`strucchange::breakpoints`, BIC, min segment ≈ 300 yr) on Glas δ¹⁸O vs depth-ordered samples.

Discharge events: local extrema of the 100-yr rate of change of the GLAC-1D total (and per region) above k·SD (k = 1, 1.5, 2), split by sign.

Quantifications, each computed for every ensemble member and for the offset scan, reported as median [5–95 %] with the fraction of members significant:
- **B1 Kernel cross-correlation** (Rehfeld et al. 2011 Gaussian kernel, bandwidth ≈ 0.25× mean step) of standardised, high-passed Glas δ¹⁸O vs total and per-region anomaly, lags ±2 kyr; null from Glas AR(1) surrogates and discharge phase-randomised surrogates; family-wise max-statistic p across regions × lags. Binned Pearson r and first-difference r as cross-checks.
- **B2 Event matching**: nearest same-sign pulse and offset per event; match count at ±200/300/500 yr; nulls (a) random placement, (b) age ensemble → match probability per event; sensitivity over k and tolerance.
- **B3 Superposed epoch analysis**: composite Glas Δδ¹⁸O in −500…+1000 yr around GLAC-1D freshening pulses (and recoveries), block-bootstrap null, computed over the age ensemble.
- **B4 Detectability**: correlation expected if Glas = standardised GLAC-1D total + AR(1) noise at the observed residual SNR (100 synthetic proxies), to calibrate what "alignment" could look like at this SNR.
- **B5 ICE-6G**: 500-yr change correlations only; state resolution limits.
- Figures: three-panel time series (Glas + ensemble spread; GLAC-1D total with pulses; ICE-6G total) with event lines; kernel-CCF with envelope; SEA composite; offset-scan curve; tables for events/offsets and sensitivity.

### C. Answer and caveats
- Direct numeric answers to each sub-question; which quantification is appropriate (event-level + SEA on an age ensemble, with kernel CCF of high-passed series; not level correlation).
- Caveats: forcing-side evaluation (no transport/AMOC-mode weighting; draft shows Bay of Biscay/Norwegian Sea reach NISA efficiently, Arctic sectors weakly, MWP-1A cancelled at NISA by AMOC recovery); reconstruction chronology uncertainty unquantified; Glas δ¹⁸O also records temperature and global-ocean drift; low proxy SNR.
- References: Endres et al. 2026a, 2026b, in prep. (draft); Rehfeld et al. 2011; Rehfeld & Kurths 2014; Franke & Donner 2019; Bühler et al. 2021; Comas-Bru et al. 2020; Haslett & Parnell 2008 (Bchron); GLAC-1D/ICE-6G/Ivanovic 2018/Romé 2022 as cited in the draft.

## Execution steps
1. Verify Table 1 event ages via Zotero read tools (search "Endres 2026 Interplay North Atlantic freshening"; fulltext/PDF pages). If unavailable, fall back to the CP DOI and say so in the manuscript.
2. Prototype Bchron on the Glas dates in the scratchpad to confirm runtime and agreement with `interp_age`; choose ensemble size accordingly (target 1,000; minimum 500).
3. Write `nisa_meltwater_sources.qmd` (overwrite): Question → Data → Methods (incl. age-uncertainty treatment) → Results A, B → Answer → Caveats → References.
4. Render: `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render nisa_meltwater_sources.qmd`; fix chunk errors; keep total runtime ≤ ~10 min (surrogates ≤ 1,000 per member where nested; use vectorised matrix ops).
5. Read back the DOCX text (`pandoc -t plain` or python-docx) to confirm every inline value resolved and prose matches numbers.

## Verification
- Render exits 0; DOCX non-trivial; all figures/tables present.
- Sanity checks hold: constant anom/Q; MWP-1A ≈ 14.6 ka; Glas 11.9–23.8 ka, 532 samples; Bchron ensemble median within reported tolerance of `interp_age`.
- Sign conventions consistent (negative anomaly = freshening; negative Glas Δδ¹⁸O = freshening).
- Re-render reproduces identical numbers (seed).
- Final report to user: headline answers, how much the age ensemble widened/changed them versus point-chronology results, and one sentence on agreement with the previous manuscript.
