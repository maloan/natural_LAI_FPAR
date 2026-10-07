# =============================================================================
# 08_agg_nonveg_0p25.R — Aggregate the water/ice-only baseline to 0.25°
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()
terraOptions(progress = 1, memfrac = 0.25)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)
area005 <- align_to_template(area005, ref005, method = "bilinear")
source_dependencies <- c(
  cfg$grids$grid_005$area_raster,
  cfg$grids$grid_005$ref_raster,
  cfg$grids$grid_025$ref_raster,
  here("config", sprintf("config_%s.yml", Sys.getenv("RUN_TAG", "alpha_0.1"))),
  here("R", "08_agg_nonveg_0p25.R"),
  here("R", "helpers", "netcdf.R"),
  here("R", "helpers", "io.R")
)

months <- format(
  seq(
    as.Date(sprintf("%d-01-01", cfg$project$years$lai_start)),
    as.Date(sprintf("%d-12-01", cfg$project$years$lai_end)),
    by = "month"
  ),
  "%Y%m"
)

for (var in c("LAI", "FPAR")) {
  input_dir <- here("output", "nonvegetated_only_0p05", var)
  output_dir <- here("output", "nonvegetated_only_0p25", var)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  input_files <- file.path(
    input_dir,
    sprintf("%s_%s_0p05_masked_nonvegetated.tif", var, months)
  )
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files)) {
    stop("Missing expected monthly files:\n", paste(missing_files, collapse = "\n"))
  }

  for (i in seq_along(input_files)) {
    output_file <- file.path(
      output_dir,
      sprintf("%s_%s_0p25_masked_nonvegetated.tif", var, months[i])
    )
    dependencies <- c(input_files[i], source_dependencies)
    if (file.exists(output_file) &&
        file.mtime(output_file) >= max(file.mtime(dependencies))) {
      next
    }

    r005 <- rast(input_files[i])
    r005 <- align_to_template(r005, ref005, method = "bilinear")

    numerator <- aggregate(
      r005 * area005,
      fact = 5,
      fun = "sum",
      na.rm = TRUE
    )
    denominator <- aggregate(
      (!is.na(r005)) * area005,
      fact = 5,
      fun = "sum",
      na.rm = TRUE
    )
    r025 <- ifel(denominator == 0, NA, numerator / denominator)
    r025 <- align_to_template(r025, ref025, method = "near")

    writeRaster(
      r025,
      output_file,
      overwrite = TRUE,
      wopt = wopt_f32(FALSE)
    )
  }
}

gc()
