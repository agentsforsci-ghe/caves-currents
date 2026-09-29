# meltmodel_site_weights.csv

Transport weights from source region to proxy site, taken from Endres et al.
(2026b), "Tracing meltwater from northern ice sheets to palaeoclimate archives
during the early last deglaciation: a conservative tracer approach" (ESS Open
Archive, doi:10.22541/essoar.15007256/v1).

| Column | Meaning |
|---|---|
| `site` | Proxy site as named in the paper (9 of its 10 sites; Llarga Cave has no moisture-uptake mask and is also absent from the paper's Fig. 7) |
| `type` | Marine sediment core, speleothem or ice core |
| `lat`, `lon` | Site coordinates (°N, °E) |
| `amoc_mode` | HadCM3 AMOC mode of the dye simulations: `cold` (weak) or `zonal` (strong) |
| `region` | Dye input region (names as in the paper) |
| `weight` | Dimensionless. Site δ¹⁸O anomaly per unit source-region δ¹⁸O anomaly, after 500 yr of constant dye input in that region |

**How the weights were obtained.** The paper's figure-scripts repository
(`dyetracer_palaeo_figurescripts`, commit d58001f, Zenodo
10.5281/zenodo.21828703) caches the per-site, per-dye δ¹⁸O contribution as
`proxymag[scenario][mode][site][dye] = weight × source anomaly`. Each weight is
that contribution divided by the scenario's source-region anomaly (−35 ‰
end-member). The 18.2k and 20.7k scenarios use the same dye fields and give
identical weights (maximum relative difference 2 × 10⁻¹⁶). Marine sites are
latitude-weighted means over a ±2° box. Land sites are weighted by the
moisture-uptake mask of each AMOC mode.
