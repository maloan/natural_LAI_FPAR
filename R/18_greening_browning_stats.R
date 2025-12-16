
## =============================================================================
## 18_greening_browning_stats.R — Fractions of greening vs browning
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(scales)
})

ROOT    <- here::here()
IN_GEO  <- file.path(ROOT, "analysis", "unmasked", "0p25")
OUTDIR  <- file.path(ROOT, "analysis", "greening_browning")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

VARS    <- c("LAI", "FPAR")
METRICS <- c("yearmean", "yearmax")
MASKS   <- c("CCI")
TAU     <- c("tau_0.05", "tau_0.1", "tau_0.2")
# Classification threshold: sign-only greening/browning
EPS <- 0

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
load_slope_df <- function(filename) {
  if (!file.exists(filename)) {
    warning("Missing slope file: ", filename)
    return(NULL)
  }
  r <- rast(filename)
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  stopifnot(ncol(df) == 3)
  colnames(df) <- c("lon", "lat", "slope")
  df
}

add_area_weight <- function(df) {
  df |>
    mutate(
      weight = cos(lat * pi / 180),
      lat_band = round(lat)
    )
}

summarise_greening_browning <- function(df, eps, group_vars = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  df <- df |>
    filter(is.finite(slope), is.finite(weight))

  if (!is.null(group_vars)) {
    grp <- rlang::syms(group_vars)
  } else {
    grp <- NULL
  }

  df |>
    mutate(class = case_when(
      slope >  eps  ~ "Greening",
      slope < -eps  ~ "Browning",
      TRUE          ~ "Neutral"
    )) |>
    group_by(!!!grp) |>
    summarise(
      w_total   = sum(weight, na.rm = TRUE),
      w_green   = sum(weight[class == "Greening"],  na.rm = TRUE),
      w_brown   = sum(weight[class == "Browning"],  na.rm = TRUE),
      w_neutral = sum(weight[class == "Neutral"],   na.rm = TRUE),
      frac_green   = w_green   / w_total,
      frac_brown   = w_brown   / w_total,
      frac_neutral = w_neutral / w_total,
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# Main loops
# -----------------------------------------------------------------------------

res_global <- list()
res_zonal  <- list()

for (var in VARS) {
  for (met in METRICS) {
    message("=== ", var, " / ", met, " ===")

    ## ---- unmasked -----------------------------------------------------------
    file_geo <- file.path(IN_GEO,
                          sprintf("%s_georef_%s_trend_slope_peryear_0p25.nc", var, met))
    df_geo <- load_slope_df(file_geo)
    if (!is.null(df_geo)) {
      df_geo <- add_area_weight(df_geo)

      g_geo <- summarise_greening_browning(df_geo, EPS, group_vars = NULL) |>
        mutate(
          VAR    = var,
          METRIC = met,
          TAU    = NA_character_,
          MASK   = NA_character_,
          domain = "unmasked"
        )

      z_geo <- summarise_greening_browning(df_geo, EPS, group_vars = "lat_band") |>
        mutate(
          VAR    = var,
          METRIC = met,
          TAU    = NA_character_,
          MASK   = NA_character_,
          domain = "unmasked"
        )

      res_global[[length(res_global) + 1]] <- g_geo
      res_zonal[[length(res_zonal) + 1]]   <- z_geo
    }

    ## ---- masked (per tau & mask) -------------------------------------------
    for (tau in TAU) {
      for (mask in MASKS) {
        file_mask <- file.path(
          ROOT,
          "output",
          tau,
          "eval",
          sprintf("trend_%s_%s", var, mask),
          sprintf("%s_%s_trend_slope_peryear_0p25.nc", var, met)
        )

        df_mask <- load_slope_df(file_mask)
        if (is.null(df_mask))
          next

        df_mask <- add_area_weight(df_mask)

        g_mask <- summarise_greening_browning(df_mask, EPS, group_vars = NULL) |>
          mutate(
            VAR    = var,
            METRIC = met,
            TAU    = tau,
            MASK   = mask,
            domain = "masked"
          )

        z_mask <- summarise_greening_browning(df_mask, EPS, group_vars = "lat_band") |>
          mutate(
            VAR    = var,
            METRIC = met,
            TAU    = tau,
            MASK   = mask,
            domain = "masked"
          )

        res_global[[length(res_global) + 1]] <- g_mask
        res_zonal[[length(res_zonal) + 1]]   <- z_mask
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Save summaries
# -----------------------------------------------------------------------------
df_global <- bind_rows(res_global)
df_zonal  <- bind_rows(res_zonal)

write_csv(df_global, file.path(OUTDIR, "greening_browning_global.csv"))
write_csv(df_zonal,
          file.path(OUTDIR, "greening_browning_zonal_1deg.csv"))

message("Saved greening/browning summaries to: ", OUTDIR)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

ROOT   <- here::here()
INDIR  <- file.path(ROOT, "analysis", "greening_browning")
OUTDIR <- file.path(INDIR, "plots")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  file.path(INDIR, "greening_browning_zonal_1deg.csv"),
  show_col_types = FALSE
)

# --- reshape to long format ---------------------------------------------------
df_long <- df |>
  tidyr::pivot_longer(
    cols = c(frac_green, frac_brown, frac_neutral),
    names_to = "class",
    values_to = "fraction"
  ) |>
  mutate(
    class = recode(
      class,
      frac_green   = "Greening",
      frac_brown   = "Browning",
      frac_neutral = "Neutral"
    )
  )

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text  = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold")
    )
}

lab_deg <- scales::label_number(suffix = "°")

# Split once; no filtering inside the loop (avoids TAU/global confusion)
df_list <- df_long |>
  group_by(VAR, METRIC, domain, MASK, TAU) |>
  group_split()

for (d in df_list) {
  meta <- d |> slice(1)

  p <- ggplot(d, aes(fraction, lat_band, colour = class)) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(labels = lab_deg) +
    labs(
      x = "Area fraction",
      y = "Latitude (°)",
      title = sprintf(
        "Zonal fractions: %s %s (%s, MASK=%s, τ=%s)",
        meta$VAR, meta$METRIC, meta$domain, meta$MASK, meta$TAU
      ),
      colour = NULL
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

  outf <- sprintf(
    "%s_%s_%s_%s_%s_zonal_fractions.png",
    meta$VAR, meta$METRIC, meta$domain, meta$MASK, meta$TAU
  )

  ggsave(file.path(OUTDIR, outf), p, width = 6.5, height = 5.0, dpi = 330)
}

cat("Saved greening/browning zonal fraction plots to:\n  ", OUTDIR, "\n")
