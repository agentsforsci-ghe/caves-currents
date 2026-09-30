# Meltwater Discharge Explorer (v2)

An interactive page for the deglacial meltwater forcing and where its
freshwater goes:

- **Left: the globe.** The surface δ¹⁸O anomaly spreads decade by decade over
  the North Atlantic and Arctic (21.5–10 ka). It comes from convolving a
  forcing (GLAC-1D or ICE-6G) with the HadCM3 pulse responses, along a cold,
  zonal or mixed AMOC pathway. The outlines are the exact HadCM3 cells where
  each region's meltwater enters the ocean.
- **Right: one site.** Pick one of the nine proxy sites of Endres et al.
  (2026b) to see its anomaly split by source region, with the totals of all
  three pathways for comparison. The tenth option shows the regional
  discharge itself.

The shared page is https://claude.ai/artifact/HMQMC2ACsbpvvDQNyRgcPG, which is
private until shared. It still shows v1 until the real kernels are in.

## Build

The page is generated from the convolution pipeline's outputs. From the
repository root:

```sh
python3 scripts/run_convolution.py      # site series + decadal fields (about 30 s)
python3 scripts/export_app_frames.py    # fields -> app/frames/*.bin + meta.json
python3 app/build_app.py                # -> app/discharge_explorer.html
```

The map frames are separate files, so the page must be served over http
(`file://` shows a message instead of the map):

```sh
cd app && python3 -m http.server 8765 --bind 127.0.0.1   # open http://127.0.0.1:8765/discharge_explorer.html
```

Stop the server afterwards.

**Placeholder kernels.** Until `data/kernels/pulse_kernels_surface.nc`
exists, the pipeline uses placeholder kernels. The page then shows a warning
banner, and `build_app.py --publish-out` refuses to write the Artifact copy.
Do not commit a page built from placeholder outputs.

**Publishing.** Run `python3 app/build_app.py --publish-out <path>`, then
republish that file to the same artifact URL to keep the link. Pass the six
`frames/<forcing>_<pathway>.bin` files as the artifact's `files`, each about
5 MB. The build prints the mapping.

## Files

| File | Role |
|---|---|
| `template.html` | Page markup, styles and script, with a data placeholder |
| `build_app.py` | Reads the pipeline outputs, the discharge files and the geometry, and writes the page. Standard library only |
| `discharge_explorer.html` | Generated page; do not edit by hand |
| `frames/` | Generated 8-bit frames and `meta.json` (`scripts/export_app_frames.py`); not in git |
| `regions_hadcm3.geojson` | Exact HadCM3 input cells per region (`scripts/make_region_geojson.py`) |
| `regions_approx.geojson` | v1 hand-digitised outlines; now only a source of region names and drainage notes |
| `vendor/land-110m.json` | Natural Earth 110 m land, from the world-atlas package (ISC licence) |

## Notes

- **Map colours.** They show the size of the negative (freshening) anomaly,
  on one logarithmic scale for all forcings and pathways, from 0.005 ‰ to the
  largest value. Positive values reach at most about 0.01 ‰ and show as zero.
- **Chart sampling.** Site series are decadal in the pipeline. The chart keeps
  every second step (20 yr) between 24 and 8 ka, to keep the page small.
- **Colours.** Region colours follow the paper figures, adjusted in lightness
  to pass colour-vision checks in light and dark themes. Keep them in sync
  with `reg_col` in `nisa_meltwater_sources.qmd`. The chart stacks the regions
  in an order that keeps neighbouring colours distinguishable.
- **Merid kernel.** The merid pulse kernel comes from a different parent run
  (xpram) than the paper's merid run. It only affects the mixed pathway at
  14.7–13.7 ka. See `data/kernels/README.md`.
