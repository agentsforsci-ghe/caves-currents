# caves-currents

Where did the deglacial North Atlantic meltwater go, and what did the NISA
stalagmite **Glas** (La Vallina cave, northern Iberia) record of it? This
repository answers the questions in `questions.md` with Quarto reports, an
interactive explorer and a forward model. The model convolves ice-sheet
meltwater with HadCM3 dye-tracer pulse responses.

Laura Endres. Repository conventions for contributors and coding agents are in
[`CLAUDE.md`](CLAUDE.md).

## Results

| Question | Where |
|---|---|
| 1. Age range of Glas | `glas_age_range.qmd` / `.docx` |
| 2. Most important and most variable meltwater source at NISA; alignment with abrupt Glas δ¹⁸O shifts (equilibrium approximation, v1) | `nisa_meltwater_sources.qmd` / `.docx` |
| 2 (v2). Source-region contributions at nine proxy sites from the forward model, with spider figures | `site_anomalies_v2.qmd` / `.docx` |
| 3. Age range of the Lake Gerzensee record | open (no data yet) |

These pages are private links that the owner can share:

- [Meltwater Discharge Explorer](https://claude.ai/artifact/HMQMC2ACsbpvvDQNyRgcPG):
  - the surface anomaly spreading decade by decade, with a GLAC-1D − ICE-6G difference view;
  - the anomaly at nine proxy sites;
  - a "How the model works" panel. Source in `app/`.
- [Meltwater Forward Model](https://claude.ai/artifact/CqfLtTKp2MCwDi9L8YFa6x): a one-page schematic of the model.
- [NISA Spider Slides](https://claude.ai/artifact/TeDmKyFH5oj1oW9dpwhAEW): how to read a spider, the nine NISA spiders, and the takeaways.

## The forward model in brief

1. **Forcing.** Regional meltwater discharge from GLAC-1D or ICE-6G
   (`data/Discharge_*.csv`) becomes a source-region δ¹⁸O anomaly per year,
   interpolated to 10-yr steps.
2. **Pulse responses.** HadCM3 surface dye after a 10-yr pulse from each of
   the nine input regions, for three AMOC modes (cold, zonal, meridional).
   They are exported on the Leeds server with `scripts/export_pulse_kernels.py`
   (see `data/kernels/README.md`).
3. **Convolution.** Each decade of forcing launches a scaled copy of the pulse
   response. The mixed pathway switches modes on a deglacial schedule.
4. **Sites.** Ocean cores take a ±2° box mean. Caves and the ice core weight the
   field by the moisture uptake of their rain, from back-trajectories.

**Open question:** should the forcing be multiplied by 10, as in the 2025.05
notebook? The pipeline and reports use ×1, and the explorer can switch to ×10.
See "Open questions" in `CLAUDE.md`.

## Running it

Python with xarray and scipy (conda `base` here), and R 4.6 with the
tidyverse packages. Quarto is the copy bundled with RStudio.

```sh
python3 scripts/make_uptake_weights.py      # once: trajectories -> data/meltmodel/uptake_weights.nc
python3 scripts/make_region_geojson.py      # once: HadCM3 input cells -> app/regions_hadcm3.geojson
python3 scripts/run_convolution.py          # site series + decadal fields -> outputs/convolution/ (about 30 s)
python3 scripts/export_app_frames.py        # fields -> app/frames/ for the explorer
python3 app/build_app.py                    # -> app/discharge_explorer.html (serve app/ over http to view)
quarto render site_anomalies_v2.qmd         # the v2 report (about 10 s)
python3 scripts/make_spider_slides.py --out <folder>   # the three spider slides
```

The pulse-kernel bundle in `data/kernels/` is not in git. Without it,
`run_convolution.py` uses a clearly flagged placeholder that must never be
reported.

## Layout

| Path | Contents |
|---|---|
| `data/` | SISAL-format Glas data, the two discharge series, the meltmodel site weights, the land-site uptake weights and the kernel README |
| `myconfig/` | Dyes and colours, sites, forcings, pathways, paths. A new forcing or site is one entry |
| `mymodules/` | Forcing loader, site extraction, kernels, convolution |
| `scripts/` | Pipeline steps, the Leeds export, the schematic and slide generators |
| `R/site_metrics.R` | The five site metrics and spider plots, reusable from any report |
| `app/` | The explorer (see `app/README.md`) |
| `outputs/convolution/` | Site series, scaling check and run info. The large field files are not in git |
| `plans/`, `prompts/` | Plans written before large tasks, and the archive of prompts behind Claude-assisted commits |

## Sources

- Endres et al. (2026a), Clim. Past 22, 797–824: the Glas record.
- Endres et al. (2026b), "Tracing meltwater…": dye regions, proxy sites and
  trajectories. Figure scripts and data: Zenodo 10.5281/zenodo.21828703.
- Endres et al. (in prep.): the impulse-response forward model.
