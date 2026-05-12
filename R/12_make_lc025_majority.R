# ==============================================================================
# make_lc025_majority.R
# Generate 0.25° majority landcover maps from ESACCI raw data (300m)
#
# Workflow:
#   1. Load annual ESACCI landcover GeoTIFFs (300m/~1km native resolution)
#   2. Rasterize to reference 0.25° grid
#   3. Compute majority class using focal statistics
#   4. Save as lc025_majority_YYYY.tif
#
# Outputs:
#   analysis/tmp/lc025_majority_yearly/lc025_majority_YYYY.tif (1992-2020)
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

terraOptions(progress = 1, memfrac = 0.6, todisk = TRUE)

# ==============================================================================
# Configuration
# ==============================================================================

year_start <- 1992L
year_end <- 2020L
years <- year_start:year_end

# Paths
esacci_dir <- here("data-raw", "ESACCI", "ESACCI_1992-2020")
out_dir <- here("analysis", "tmp", "lc025_majority_yearly")
ref_025 <- rast(here("src", "ref_0p25.nc"))

# Create output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(esacci_dir)) {
  stop("ESACCI data directory not found: ", esacci_dir, call. = FALSE)
}

if (!file.exists(here("src", "ref_0p25.nc"))) {
  stop("Reference grid not found: src/ref_0p25.nc", call. = FALSE)
}

# ==============================================================================
# Main workflow
# ==============================================================================

for (year in years) {
  cat(sprintf("Processing ESACCI landcover %d...", year))

  # Find ESACCI file for this year
  # Pattern: ESACCI-LC-L4-LCCS-Map-300m-P1Y-YYYY-v2.0.7.tif or
  #          C3S-LC-L4-LCCS-Map-300m-P1Y-YYYY-v2.1.1.tif
  esacci_files <- list.files(
    esacci_dir,
    pattern = sprintf(".*%d.*\\.tif$", year),
    full.names = TRUE
  )

  if (length(esacci_files) == 0) {
    warning(sprintf("No ESACCI file found for year %d, skipping", year))
    next
  }

  if (length(esacci_files) > 1) {
    # Use the one with more recent version if multiple found
    esacci_files <- sort(esacci_files, decreasing = TRUE)[1]
  }

  esacci_file <- esacci_files[1]
  cat(sprintf(" (file: %s)\n", basename(esacci_file)))

  # Load ESACCI landcover
  lc_native <- rast(esacci_file)

  if (is.null(lc_native)) {
    warning(sprintf("Failed to load ESACCI data for year %d", year))
    next
  }

  # ===========================================================================
  # Rasterize to 0.25° with majority class
  # ===========================================================================

  # Method: use focal statistics to compute majority class at 0.25° resolution
  # First, we need to establish the aggregation factor

  # Get the native resolution (should be ~300m, ~0.003 degrees)
  native_res <- res(lc_native)[1]
  target_res <- res(ref_025)[1]

  # Aggregation factor (e.g., 300m / 0.003deg ≈ 100 cells per side)
  agg_factor <- round(target_res / native_res)

  if (agg_factor < 2) {
    warning(
      sprintf(
        "Native resolution %f is already finer than or equal to target 0.25° (factor: %d)",
        native_res, agg_factor
      )
    )
    agg_factor <- 10 # Use reasonable default
  }

  # Compute majority class via aggregation
  lc_agg <- aggregate(
    lc_native,
    fact = agg_factor,
    fun = "modal",
    na.rm = TRUE
  )

  # Reproject to match reference grid exactly
  lc_025 <- project(
    lc_agg,
    ref_025,
    method = "near"
  )

  # Ensure same extent and resolution
  lc_025 <- resample(lc_025, ref_025, method = "near")

  # ===========================================================================
  # Save output
  # ===========================================================================

  out_file <- file.path(
    out_dir,
    sprintf("lc025_majority_%d.tif", year)
  )

  writeRaster(
    lc_025,
    filename = out_file,
    overwrite = TRUE,
    gdal = c(
      "COMPRESS=DEFLATE",
      "ZLEVEL=6"
    ),
    datatype = "INT2U",
    NAflag = 0
  )

  cat(sprintf("  -> Saved to %s\n", basename(out_file)))
}

cat(sprintf("\nComplete! Generated %d files in %s\n", length(years), out_dir))
