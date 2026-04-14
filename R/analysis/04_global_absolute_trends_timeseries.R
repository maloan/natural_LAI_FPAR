# ==============================================================================
# 04_global_absolute_trends_timeseries.R
# Figure 3: Global annual absolute LAI time series and linear trends
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(here)
})

# ---- config ------------------------------------------------------------------
var <- "LAI"
metrics <- c("yearmean", "yearmax")
year0 <- 1982L
show_raw <- TRUE

taus_cci <- c("tau_0.05", "tau_0.1", "tau_0.2")
tau_glc <- "tau_0.1"

outdir <- here("analysis", "results", "paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

area_path <- here("src", "area_0p25_validdomain_km2.nc")
if (!file.exists(area_path)) {
  stop("Missing valid-domain area raster: ", area_path)
}
area <- rast(area_path)[[1]]

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    )
}

path_unmasked <- function(metric) {
  here("analysis", "unmasked", "0p25", sprintf("%s_georef_%s_0p25.nc", var, metric))
}
path_cci <- function(metric, tau) {
  here("output", tau, "eval", sprintf("trend_%s_CCI", var), sprintf("%s_%s_0p25.nc", var, metric))
}
path_glc <- function(metric, tau) {
  here("output", tau, "eval", sprintf("trend_%s_GLC", var), sprintf("%s_%s_0p25.nc", var, metric))
}

series_wmean <- function(r, area, year0 = 1982L) {
  compareGeom(r, area, stopOnError = TRUE)
  yrs <- year0:(year0 + nlyr(r) - 1L)
  vals <- vector("numeric", length = nlyr(r))

  for (i in seq_len(nlyr(r))) {
    x <- r[[i]]
    ok <- is.finite(x) & is.finite(area) & area > 0
    num <- as.numeric(global(ifel(ok, x * area, NA), "sum", na.rm = TRUE)[1, 1])
    den <- as.numeric(global(ifel(ok, area, NA), "sum", na.rm = TRUE)[1, 1])
    vals[i] <- ifelse(is.finite(den) && den > 0, num / den, NA_real_)
  }

  tibble(year = yrs, value = vals)
}

rows <- list()
for (metric in metrics) {
  fU <- path_unmasked(metric)
  if (!file.exists(fU)) stop("Missing unmasked file: ", fU)
  rows[[length(rows) + 1]] <-
    series_wmean(rast(fU), area, year0) |>
    mutate(metric = metric, scenario = "Unmasked")

  for (tau in taus_cci) {
    fC <- path_cci(metric, tau)
    if (!file.exists(fC)) stop("Missing CCI file: ", fC)
    rows[[length(rows) + 1]] <-
      series_wmean(rast(fC), area, year0) |>
      mutate(metric = metric, scenario = sprintf("CCI %s", gsub("tau_", "tau=", tau)))
  }

  fG <- path_glc(metric, tau_glc)
  if (!file.exists(fG)) stop("Missing GLC file: ", fG)
  rows[[length(rows) + 1]] <-
    series_wmean(rast(fG), area, year0) |>
    mutate(metric = metric, scenario = "GLC")
}

scenario_levels <- c("Unmasked", sprintf("CCI tau=%s", sub("^tau_", "", taus_cci)), "GLC")

df <- bind_rows(rows) |>
  mutate(
    metric = factor(metric, levels = c("yearmean", "yearmax"), labels = c("Annual mean", "Annual maximum")),
    scenario = factor(scenario, levels = scenario_levels)
  )

x_breaks <- seq(floor(min(df$year, na.rm = TRUE) / 7) * 7, max(df$year, na.rm = TRUE), by = 7)

# Fit simple OLS trend per metric/scenario
trend_df <- df |>
  group_by(metric, scenario) |>
  group_modify(~ {
    d <- .x |> filter(is.finite(value), is.finite(year))
    if (nrow(d) < 3) {
      return(tibble(year = d$year, fit = NA_real_, slope = NA_real_))
    }
    mod <- lm(value ~ year, data = d)
    tibble(
      year = d$year,
      fit = predict(mod, newdata = d),
      slope = unname(coef(mod)[["year"]])
    )
  }) |>
  ungroup()

col_map <- c(
  "Unmasked" = "grey20",
  "CCI tau=0.05" = "#1b9e77",
  "CCI tau=0.1" = "#d95f02",
  "CCI tau=0.2" = "#7570b3",
  "GLC" = "#386cb0"
)

p <- ggplot() +
  {
    if (isTRUE(show_raw)) {
      geom_line(
        data = df,
        aes(year, value, colour = scenario),
        linewidth = 0.4,
        alpha = 0.4,
        na.rm = TRUE
      )
    }
  } +
  geom_line(
    data = trend_df,
    aes(year, fit, colour = scenario),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = col_map, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Year",
    y = "Global mean absolute LAI",
    colour = NULL,
    title = "Global absolute LAI trajectories under masking scenarios",
    subtitle = "Thin lines: annual values; thick lines: OLS trends"
  ) +
  theme_pub()

write_csv(df, file.path(outdir, "global_timeseries_absolute_trends_yearmean_yearmax.csv"))
write_csv(trend_df, file.path(outdir, "global_timeseries_absolute_trends_yearmean_yearmax_trendlines.csv"))

ggsave(
  file.path(outdir, "global_timeseries_absolute_trends_yearmean_yearmax.png"),
  p,
  width = 11, height = 5.8, dpi = 350
)
ggsave(
  file.path(outdir, "global_timeseries_absolute_trends_yearmean_yearmax.pdf"),
  p,
  width = 11, height = 5.8
)
