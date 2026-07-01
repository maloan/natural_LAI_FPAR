# ==============================================================================
# 11_zonal_diagnostics_combined_figure.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(here)
})

source(here("R", "helpers", "plotting.R"))

lat_label_fn <- get("lat_labels", mode = "function")
theme_pub_fn <- get("theme_pub", mode = "function")

tau <- "tau_0.1"
var <- "LAI"

f_rel <- here(
  "analysis",
  "results",
  "figures",
  "summaries",
  sprintf("zonal_relative_trends_all_masks_%s.csv", tau)
)
f_amp <- here(
  "analysis",
  "results",
  "figures",
  "summaries",
  sprintf(
    "zonal_yearamp_timeMean_%s_all_masks_%s.csv",
    toupper(var),
    tau
  )
)

if (!file.exists(f_rel)) {
  stop("Missing relative-trends table: ", f_rel)
}
if (!file.exists(f_amp)) {
  stop("Missing seasonal-amplitude table: ", f_amp)
}

df_rel <- read_csv(f_rel, show_col_types = FALSE)
df_amp <- read_csv(f_amp, show_col_types = FALSE)

scenario_levels <- c("Unmasked", "CCI tau=0.05", "CCI tau=0.1", "CCI tau=0.2", "GLC")
df_rel <- df_rel |>
  mutate(scenario = factor(.data$scenario, levels = scenario_levels))
df_amp <- df_amp |>
  mutate(scenario = factor(.data$scenario, levels = scenario_levels))

df_rel_mean <- df_rel |> filter(metric == "Annual mean")
df_rel_max <- df_rel |> filter(metric == "Annual maximum")
rel_ylim <- range(
  c(
    df_rel_mean$reltrend_pct_per_year,
    df_rel_max$reltrend_pct_per_year
  ),
  na.rm = TRUE
)
rel_pad <- 0.05 * diff(rel_ylim)
rel_ylim <- rel_ylim + c(-rel_pad, rel_pad)

p1 <- plot_zonal_diagnostics(
  df_rel_mean,
  "reltrend_pct_per_year",
  "Annual Mean Relative Trend",
  expression("Relative trend (% yr"^
    {
      -1
    } * ")"),
  y_limits = rel_ylim
)

p2 <- plot_zonal_diagnostics(
  df_rel_max,
  "reltrend_pct_per_year",
  "Annual Maximum Relative Trend",
  expression("Relative trend (% yr"^
    {
      -1
    } * ")"),
  y_limits = rel_ylim
)
p3 <- plot_zonal_diagnostics(
  df_amp,
  "mean_yearamp",
  "Seasonal Amplitude",
  expression("Mean seasonal amplitude LAI (" * m^2 ~ m^
    {
      -2
    } * ")")
)

# Combine panels
fig <- (p1 / p2 / p3) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.caption = element_text(
    size = 9,
    hjust = 0,
    margin = margin(t = 10)
  ))
outdir <- here("analysis", "results", "figures", "summaries")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_stub <- sprintf(
  "zonal_diagnostics_mean_max_amplitude_%s_all_masks_%s",
  toupper(var),
  tau
)
ggsave(
  file.path(outdir, paste0(out_stub, ".png")),
  fig,
  width = 10.5,
  height = 11.5,
  dpi = 400
)
ggsave(
  file.path(outdir, paste0(out_stub, ".pdf")),
  fig,
  width = 10.5,
  height = 11.5,
  device = cairo_pdf
)
