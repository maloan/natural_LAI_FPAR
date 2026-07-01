# ==============================================================================
# make_lc025_fractions.R — Generate 0.25° fractional land-cover composition from
# ESACCI 300 m data
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(terra)
})

source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()

year_start <- as.integer(cfg$project$years$cci_start)
year_end <- as.integer(cfg$project$years$cci_end)

lc_year <- Sys.getenv("LC_YEAR", unset = Sys.getenv("lc_year", unset = NA_character_))
if (!is.na(lc_year) && nzchar(lc_year)) {
  years <- as.integer(lc_year)
  if (anyNA(years)) {
    stop("Invalid lc_year: ", lc_year, call. = FALSE)
  }
} else {
  years <- year_start:year_end
}

scratch_dir <- Sys.getenv("SLURM_TMPDIR", unset = file.path(tempdir(), "terra_lc025_fractions"))

dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(
  progress = 1,
  todisk = TRUE,
  tempdir = scratch_dir
)

esacci_dir <- cfg$paths$cci_dir
ref_025 <- rast(cfg$grids$grid_025$ref_raster)

out_dir_frac <- here("analysis", "tmp", "lc025_fraction_yearly")
out_dir_majority <- here("analysis", "tmp", "lc025_majority_yearly")

dir.create(out_dir_frac, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_majority,
  recursive = TRUE,
  showWarnings = FALSE
)

merge_to_parent <- c(
  `0` = 0,
  # No Data
  `10` = 10,
  # Cropland, rainfed
  `11` = 10,
  # Herbaceous cropland
  `12` = 10,
  # Tree or shrub cropland
  `20` = 10,
  # Irrigated or post-flooded cropland
  `30` = 10,
  # Mosaic cropland
  `40` = 10,
  # Mosaic cropland
  `50` = 50,
  # broad leaved evergreen trees
  `60` = 60,
  # broad leaved deciduous trees
  `61` = 60,
  # broad leaved deciduous trees
  `62` = 60,
  # broad leaved deciduous trees
  `70` = 70,
  # Needleleaved evergreen trees
  `71` = 70,
  # Needleleaved evergreen
  `72` = 70,
  # Needleleaved deciduous
  `80` = 80,
  # Needleleaved evergreen
  `81` = 80,
  # Needleleaved deciduous
  `82` = 80,
  # Needleleaved deciduous
  `90` = 90,
  # Mixed forests
  `100` = 100,
  # Mosaic trees and shrub
  `110` = 100,
  # Mosaic trees and shrub
  `120` = 120,
  # Shrubland
  `121` = 120,
  # Evergreen shrubland
  `122` = 120,
  # Deciduous shrubland
  `130` = 130,
  # Grassland
  `140` = 140,
  # Lichen and Mosses
  `150` = 150,
  # Sparse vegetation
  `151` = 150,
  # Sparse tree
  `152` = 150,
  # sparse shrub
  `153` = 150,
  # sparse herbaceous
  `160` = 160,
  # tree cover flooded fresh
  `170` = 160,
  # tree cover flooded saline
  `180` = 180,
  # shrub or herbaceous cover flooded
  `190` = 190,
  # Urban
  `200` = 200,
  # bare areas
  `201` = 200,
  # consolidated bare areas
  `202` = 200,
  # unconsolidated bare areas
  `210` = 210,
  # water
  `220` = 220 # permanent snow and ice
)

parent_classes <- sort(unique(as.integer(merge_to_parent[merge_to_parent > 0])))

rcl <- cbind(from = as.integer(names(merge_to_parent)), to = as.integer(merge_to_parent))

for (year in years) {
  frac_file <- file.path(out_dir_frac, sprintf("lc025_fraction_%d.tif", year))

  majority_file <- file.path(out_dir_majority, sprintf("lc025_majority_%d.tif", year))

  if (file.exists(frac_file) &&
    file.exists(majority_file) &&
    !isTRUE(as.logical(Sys.getenv(
      "REMAKE_ALL", Sys.getenv("remake_all", "FALSE")
    )))) {
    message("Outputs exist, skipping year ", year)
    next
  }

  esacci_files <- list.files(esacci_dir,
    pattern = sprintf(".*%d.*\\.tif$", year),
    full.names = TRUE
  )

  if (!length(esacci_files)) {
    warning("No ESACCI file for year ", year, " — skipping")
    next
  }

  esacci_files <- sort(esacci_files)

  if (length(esacci_files) > 1L) {
    warning(
      "Multiple ESACCI files found for year ",
      year,
      "; using: ",
      basename(esacci_files[1])
    )
  }

  esacci_file <- esacci_files[1]
  message("Processing year ", year, ": ", basename(esacci_file))

  lc_native <- rast(esacci_file)

  # Merge original ESACCI classes to parent classes
  lc_parent <- classify(lc_native, rcl = rcl, others = 0)

  # Native-cell area in km2
  native_area <- cellSize(lc_parent, unit = "km")

  # Total valid land-cover area per 0.25° cell
  valid_area_native <- ifel(lc_parent > 0, native_area, NA)

  valid_area_025 <- resample(valid_area_native, ref_025, method = "sum")

  frac_layers <- vector("list", length(parent_classes))

  for (i in seq_along(parent_classes)) {
    cls <- parent_classes[i]

    message("  Class ", cls)

    class_area_native <- ifel(lc_parent == cls, native_area, NA)

    class_area_025 <- resample(class_area_native, ref_025, method = "sum")

    frac <- class_area_025 / valid_area_025
    frac <- ifel(is.finite(frac), frac, NA)

    names(frac) <- paste0("lc_", cls)
    frac_layers[[i]] <- frac

    rm(class_area_native, class_area_025, frac)
    gc()
  }

  frac_stack <- rast(frac_layers)
  names(frac_stack) <- paste0("lc_", parent_classes)

  # Write fractional composition
  writeRaster(
    frac_stack,
    filename = frac_file,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
    datatype = "FLT4S"
  )

  message("Wrote fractional composition: ", frac_file)

  idx_majority <- which.max(frac_stack)

  max_frac <- app(frac_stack, max, na.rm = TRUE)

  idx_majority <- ifel(is.finite(max_frac) &
    max_frac > 0, idx_majority, NA)

  rcl_idx <- cbind(from = seq_along(parent_classes), to = parent_classes)

  lc_majority_025 <- classify(idx_majority, rcl = rcl_idx, others = NA)

  names(lc_majority_025) <- "lc_id"

  lc_majority_025 <- align_to_template(lc_majority_025, ref_025, method = "near")

  writeRaster(
    lc_majority_025,
    filename = majority_file,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"),
    datatype = "INT2S"
  )

  message("Wrote majority class map: ", majority_file)

  rm(
    lc_native,
    lc_parent,
    native_area,
    valid_area_native,
    valid_area_025,
    frac_layers,
    frac_stack,
    idx_majority,
    max_frac,
    lc_majority_025
  )
  gc()
}
