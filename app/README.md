# Meltwater Discharge Explorer

An interactive page for the GLAC-1D and ICE-6G regional meltwater discharge
series in `data/`, with a map of the nine dye-tracer source regions and the
Glas δ¹⁸O record from NISA.

- **Open locally:** `app/discharge_explorer.html` in any browser. The data is
  inlined. An internet connection is needed only for the two script libraries
  (d3, topojson-client) and the web fonts.
- **Shared page (private, owner-only until shared):**
  https://claude.ai/artifact/HMQMC2ACsbpvvDQNyRgcPG

## Rebuild

The page is generated. Edit `template.html` or the inputs, then run from the
repository root:

```sh
python3 app/build_app.py
```

To also write the copy used for the shared page, which has no `<html>`/`<head>`
wrapper, pass `--publish-out <path>`. Republish that file to the same artifact
URL to keep the link.

## Files

| File | Role |
|---|---|
| `template.html` | Page markup, styles and script, with a data placeholder |
| `build_app.py` | Reads the CSVs, the Glas samples and the geometry, and writes the page (standard library only) |
| `discharge_explorer.html` | Generated page; do not edit by hand |
| `regions_approx.geojson` | Approximate outlines of the nine input regions and the NISA and PS2644-5 sites |
| `vendor/land-110m.json` | Natural Earth 110 m land, from the world-atlas package (ISC licence) |

## Notes

- Region outlines are hand-digitised from Endres et al. (2026b, Fig. 2a) and
  the forward-model draft (Fig. 5.1). They are not the HadCM3 masks. Replacing
  `regions_approx.geojson` with the real masks, same `code` property, needs no
  code change.
- The δ¹⁸O anomaly is the source-region anomaly before transport to NISA, as
  given in the discharge files.
- Region colours follow the paper figures, adjusted in lightness so that the
  palette passes colour-vision checks in light and dark themes. The chart stacks
  the regions in an order that keeps neighbouring colours distinguishable.
