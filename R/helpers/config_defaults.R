## =============================================================================
## config_defaults.R — Global configuration for LAI/FPAR trend analysis
##
## Single source of truth for all hardcoded parameters used across analysis scripts.
## Load this file at the top of any script to access standardized configuration.
## =============================================================================

# ==============================================================================
# TEMPORAL CONFIGURATION
# ==============================================================================

# Trend calculation period
temporal_config <- list(
  year_start = 1982L,
  year_end = 2024L,
  n_years = 2024L - 1982L + 1L, # 43 years
  lc_year_start = 1992L, # ESA-CCI land-cover availability
  lc_year_end = 2020L,
  lc_years = 1992L:2020L
)

# ==============================================================================
# MASKING CONFIGURATION
# ==============================================================================

# Cloud cover thresholds (tau parameters)
cci_tau_options <- c("tau_0.05", "tau_0.1", "tau_0.2")
default_cci_tau <- "tau_0.1"
default_glc_tau <- "tau_0.1"

mask_sources <- c("CCI", "GLC")

# Relative trend thresholds (avoid inflated trends where mean ≈ 0)
eps_thresholds <- list(
  LAI = 0.05, # LAI must be >= 0.05 for relative trend
  FPAR = 0.02 # FPAR must be >= 0.02 for relative trend
)

# ==============================================================================
# VARIABLES AND METRICS
# ==============================================================================

variables <- c("LAI", "FPAR")
metrics <- c("yearmean", "yearmax")
all_metrics <- c("yearmean", "yearmax", "yearmin", "yearamp")

# Unit labels for plots and tables
unit_labels <- list(
  LAI = list(
    absolute = "LAI units yr⁻¹",
    relative = "% yr⁻¹"
  ),
  FPAR = list(
    absolute = "fAPAR units yr⁻¹",
    relative = "% yr⁻¹"
  )
)

# ==============================================================================
# SPATIAL CONFIGURATION
# ==============================================================================

# Grid specification
grid_config <- list(
  resolution_deg = 0.25,
  n_rows = 720L,
  n_cols = 1440L,
  lat_min = -90,
  lat_max = 90,
  lon_min = -180,
  lon_max = 180,
  crs = "EPSG:4326" # WGS84 geographic
)

# Zonal aggregation
zonal_config <- list(
  latitude_band_width = 1L, # degrees
  n_latitude_bands = 180L # 180 1-degree bands
)

# ==============================================================================
# BOOTSTRAP CONFIGURATION
# ==============================================================================

bootstrap_config <- list(
  n_boot = 1000L, # Number of bootstrap samples
  confidence_level = 0.95, # 95% CI
  seed = 42L, # Fixed seed for reproducibility
  alpha = 0.05 # Significance level
)

# ==============================================================================
# MANN-KENDALL CONFIGURATION
# ==============================================================================

mk_config <- list(
  alpha = 0.05, # Significance level (p < 0.05)
  method = "Sen's slope", # Robust median slope estimator
  two_tailed = TRUE # Two-tailed test
)

# ==============================================================================
# AGGREGATION SCENARIOS (Script 02, 03, 04 etc.)
# ==============================================================================

scenario_spec <- tibble::tibble(
  scenario = c(
    "Unmasked",
    "CCI tau=0.05",
    "CCI tau=0.1",
    "CCI tau=0.2",
    "GLC"
  ),
  source = c("unmasked", "CCI", "CCI", "CCI", "GLC"),
  run_tag = c(NA_character_, "tau_0.05", "tau_0.1", "tau_0.2", "tau_0.1"),
  mask_type = c(NA_character_, "CCI", "CCI", "CCI", "GLC"),
  tau_param = c(NA_character_, 0.05, 0.1, 0.2, 0.1)
)

# ==============================================================================
# LAND-COVER CONFIGURATION (Script 12)
# ==============================================================================

# IPCC-style land-cover class aggregation mapping
lc_merge_mapping <- list(
  # Cropland + mosaics -> 10
  `11` = 10, `12` = 10, `20` = 10, `30` = 10, `40` = 10,
  # Broadleaf: 61/62 -> 60, then 60 -> 50
  `61` = 50, `62` = 50, `60` = 50,
  # Needleleaf: 71/72 -> 70; 81/82 -> 80; then 80 -> 70
  `71` = 70, `72` = 70, `81` = 70, `82` = 70, `80` = 70,
  # Mosaic herbaceous/tree -> 100
  `110` = 100,
  # Shrubland -> 120
  `121` = 120, `122` = 120,
  # Sparse vegetation -> 150
  `151` = 150, `152` = 150, `153` = 150,
  # Flooded -> 180
  `160` = 180, `170` = 180,
  # Bare areas -> 200
  `201` = 200, `202` = 200
)

lc_legend <- tibble::tibble(
  lc_id = c(10L, 50L, 70L, 90L, 100L, 120L, 130L, 140L, 150L, 180L, 190L, 200L, 210L, 220L),
  lc_name = c(
    "Cropland (incl. mosaics)",
    "Tree cover, broadleaved",
    "Tree cover, needleleaved",
    "Tree cover, mixed",
    "Mosaic tree/shrub and herbaceous",
    "Shrubland",
    "Grassland",
    "Lichens and mosses",
    "Sparse vegetation (<15%)",
    "Flooded vegetation / wetlands",
    "Urban areas",
    "Bare areas",
    "Water bodies",
    "Permanent snow and ice"
  )
)

# ==============================================================================
# KÖPPEN-GEIGER CONFIGURATION (Script 13)
# ==============================================================================

kg_resolution <- "course" # "course" or "fine"
kg_lookup_chunk_size <- 50000L # Process pixels in chunks

kg_names_2letter <- list(
  "Af" = "Tropical rainforest",
  "Am" = "Tropical monsoon",
  "As" = "Tropical savanna (dry summer)",
  "Aw" = "Tropical savanna (dry winter)",
  "BW" = "Desert",
  "BS" = "Steppe",
  "Cf" = "Temperate, fully humid",
  "Cl" = "Unknown climate zone",
  "Cs" = "Temperate, dry summer",
  "Cw" = "Temperate, dry winter",
  "Df" = "Cold, fully humid",
  "Ds" = "Cold, dry summer",
  "Dw" = "Cold, dry winter",
  "ET" = "Tundra",
  "EF" = "Ice cap"
)

# ==============================================================================
# OUTPUT CONFIGURATION
# ==============================================================================

output_config <- list(
  base_dir = "analysis/results",
  figure_format = c("png", "pdf"),
  png_dpi = 320,
  png_width = 1100,
  png_height = 550,
  csv_decimal_places = list(
    absolute_trend = 5,
    relative_trend = 4,
    ci_bounds = 4
  )
)

# ==============================================================================
# HELPER FUNCTION: Get configuration value with validation
# ==============================================================================

get_config <- function(key, default = NULL) {
  # Safely retrieve configuration values with fallback to default
  # Usage: get_config("bootstrap_config$n_boot")

  tryCatch(
    {
      eval(parse(text = key))
    },
    error = function(e) {
      if (is.null(default)) {
        stop(sprintf("Configuration key not found: %s", key))
      }
      default
    }
  )
}

# ==============================================================================
# Validation function: ensure config matches data
# ==============================================================================

validate_config_with_data <- function(r, var, mask = NULL) {
  # Check that raster matches expected configuration
  # Args:
  #   r: SpatRaster to validate
  #   var: "LAI" or "FPAR"
  #   mask: mask type (optional, for informational output)
  # Returns: invisible(r) if valid

  if (nrow(r) != grid_config$n_rows || ncol(r) != grid_config$n_cols) {
    stop(
      sprintf(
        "Raster dimensions mismatch: got %d×%d, expected %d×%d",
        nrow(r), ncol(r), grid_config$n_rows, grid_config$n_cols
      )
    )
  }

  if (!var %in% variables) {
    stop(sprintf("Unknown variable: %s (must be one of: %s)", var, paste(variables, collapse = ", ")))
  }

  invisible(r)
}

# ==============================================================================
# Print configuration summary
# ==============================================================================

print_config_summary <- function() {
  cat(
    "\n=== LAI/FPAR TREND ANALYSIS CONFIGURATION ===\n\n",
    "TEMPORAL:\n",
    sprintf("  Trend period: %d-%d (%d years)\n", temporal_config$year_start, temporal_config$year_end, temporal_config$n_years),
    sprintf("  LC stability: %d-%d\n\n", temporal_config$lc_year_start, temporal_config$lc_year_end),
    "SPATIAL:\n",
    sprintf("  Grid: %s, 0.25° resolution (%d × %d pixels)\n", grid_config$crs, grid_config$n_rows, grid_config$n_cols),
    sprintf("  Zonal bands: %d × %d° latitude\n\n", zonal_config$n_latitude_bands, zonal_config$latitude_band_width),
    "STATISTICAL:\n",
    sprintf("  Bootstrap: n=%d, confidence=%.0f%%, seed=%d\n", bootstrap_config$n_boot, bootstrap_config$confidence_level * 100, bootstrap_config$seed),
    sprintf("  Mann-Kendall: α=%.2f, %s\n\n", mk_config$alpha, mk_config$method),
    "SCENARIOS:\n",
    sprintf("  CCI taus: %s\n", paste(cci_tau_options, collapse = ", ")),
    sprintf("  Masks: %s\n\n", paste(mask_sources, collapse = ", ")),
    "RELATIVE TREND THRESHOLDS (avoid 0/0):\n",
    sprintf("  LAI: ≥%.2f\n", eps_thresholds$LAI),
    sprintf("  FPAR: ≥%.2f\n", eps_thresholds$FPAR)
  )
  invisible(NULL)
}
