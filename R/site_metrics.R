# Site metrics and spider plots for the convolved meltwater anomalies.
#
# Reusable by any Quarto report: source("R/site_metrics.R").
# Input: outputs/convolution/site_anomaly.csv.gz from scripts/run_convolution.py
# (one row per forcing, pathway, site and decade; one column per region code).
#
# The five metrics follow nisa_meltwater_sources.qmd (v1), computed per
# forcing x pathway x site on a common 500-yr grid:
#   share   mean anomaly of the region / sum of the regional means
#   dom     fraction of time steps in which the region is the largest contributor
#   sd      temporal standard deviation of the region's anomaly (per mil)
#   vshare  Cov(region, total) / Var(total); sums to 1 over the regions
#   sd_diff standard deviation of the 500-yr changes (per mil), a measure of abruptness

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

# Region code (discharge files) -> name in Endres et al. (2026b), in dye order.
region_names <- c(Med = "EIS_MedSea", Bri = "EIS_BayOfBiscay", Fen = "EIS_NorwegianSea",
                  EurArc = "EIS_Arctic", AmeArc = "LAU_Arctic", GIS = "GIS_GreenlandSea",
                  NLau = "LAU_LabradorSea", SLau = "LAU_StLawrence", GulofMex = "LAU_GulfOfMexico")

# Region colours: the paper figures' colours, lightness-tuned to pass colour-vision
# checks; identical to reg_col in nisa_meltwater_sources.qmd and app/template.html.
reg_col <- c(EIS_MedSea = "#1f77b4", EIS_BayOfBiscay = "#f1780b", EIS_NorwegianSea = "#2ca02c",
             EIS_Arctic = "#d62728", LAU_Arctic = "#9467bd", GIS_GreenlandSea = "#9b4c3c",
             LAU_LabradorSea = "#de72bd", LAU_StLawrence = "#a5a608", LAU_GulfOfMexico = "#00b3c3")

forcing_labels <- c(glac1d = "GLAC-1D", ice6g = "ICE-6G")
rec_col <- c("GLAC-1D" = "#184f95", "ICE-6G" = "#c26a00")

metric_labels <- c(share = "Share", dom = "Dominant", sd = "SD", vshare = "Var. share", sd_diff = "SD Δ")

#' Read the convolved site series as a long table (one row per region and decade).
read_site_anomaly <- function(path = "outputs/convolution/site_anomaly.csv.gz") {
  read_csv(path, show_col_types = FALSE) |>
    select(-total) |>
    pivot_longer(all_of(names(region_names)), names_to = "code", values_to = "anom") |>
    mutate(region = factor(region_names[code], levels = unname(region_names)),
           forcing = recode(forcing, !!!forcing_labels))
}

#' Average onto a common grid of `step` years (bins centred on multiples of `step`),
#' keeping ages inside `window` (oldest, youngest; yr BP).
to_grid <- function(d, step = 500, window = c(21500, 10000)) {
  d |> filter(time_bp <= window[1], time_bp >= window[2]) |>
    mutate(age = round(time_bp / step) * step) |>
    group_by(forcing, pathway, site, region, age) |>
    summarise(anom = mean(anom), .groups = "drop")
}

#' The five metrics per forcing x pathway x site x region.
site_metrics <- function(g) {
  g |>
    group_by(forcing, pathway, site, age) |>
    mutate(tot = sum(anom), dom = anom == min(anom)) |>          # largest contributor = most negative
    ungroup() |>
    arrange(forcing, pathway, site, region, desc(age)) |>
    group_by(forcing, pathway, site, region) |>
    mutate(dA = anom - lag(anom)) |>
    summarise(mean = mean(anom), sd = sd(anom), dom = mean(dom),
              vshare = cov(anom, tot) / var(tot), sd_diff = sd(dA, na.rm = TRUE), .groups = "drop") |>
    group_by(forcing, pathway, site) |>
    mutate(share = mean / sum(mean)) |>
    ungroup()
}

#' Spider plot for one site: one panel per region, one polygon per forcing.
#' Each spoke is scaled to the largest value of that metric at this site
#' (over all regions and both forcings), so the outer ring is the site maximum.
#' Negative shares and variance shares are drawn at zero.
spider_plot <- function(m, site_name, pathway_name = "mixed", ncol = 3) {
  keys <- names(metric_labels)
  rad <- m |> filter(site == site_name, pathway == pathway_name) |>
    select(forcing, region, all_of(keys)) |>
    pivot_longer(all_of(keys), names_to = "metric", values_to = "value") |>
    mutate(value = pmax(value, 0)) |>
    group_by(metric) |> mutate(r = if (max(value) > 0) value / max(value) else 0) |> ungroup() |>
    mutate(k = match(metric, keys), ang = pi / 2 - 2 * pi * (k - 1) / length(keys),
           x = r * cos(ang), y = r * sin(ang)) |>
    arrange(forcing, region, k)
  closed <- rad |> group_by(forcing, region) |> group_modify(~ bind_rows(.x, .x[1, ])) |> ungroup()
  spokes <- tibble(k = seq_along(keys), ang = pi / 2 - 2 * pi * (k - 1) / length(keys),
                   x = cos(ang), y = sin(ang), lab = unname(metric_labels),
                   lx = 1.28 * cos(ang), ly = 1.2 * sin(ang))
  rings <- expand_grid(r = c(0.5, 1), t = seq(0, 2 * pi, length.out = 73)) |>
    mutate(x = r * cos(t), y = r * sin(t))
  ggplot() +
    geom_path(data = rings, aes(x, y, group = r), colour = "grey85", linewidth = 0.3) +
    geom_segment(data = spokes, aes(0, 0, xend = x, yend = y), colour = "grey85", linewidth = 0.3) +
    geom_text(data = spokes, aes(lx, ly, label = lab), size = 2.3, colour = "grey35") +
    geom_polygon(data = closed, aes(x, y, fill = forcing, group = forcing), alpha = 0.12, colour = NA) +
    geom_path(data = closed, aes(x, y, colour = forcing, group = forcing), linewidth = 0.6) +
    facet_wrap(~region, ncol = ncol, drop = FALSE) +
    scale_colour_manual(values = rec_col, name = NULL, aesthetics = c("colour", "fill")) +
    coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.2, 1.3)) +
    theme_void(base_size = 9) +
    theme(strip.text = element_text(face = "bold", size = 8.5, margin = margin(b = 2)),
          legend.position = "bottom", panel.spacing = unit(4, "pt"))
}

#' Most important and most variable region per forcing x pathway x site.
site_summary <- function(m) {
  m |> group_by(forcing, pathway, site) |>
    summarise(important = as.character(region[which.max(share)]), share = max(share),
              variable = as.character(region[which.max(sd)]), sd = max(sd),
              abrupt = as.character(region[which.max(sd_diff)]), .groups = "drop")
}
