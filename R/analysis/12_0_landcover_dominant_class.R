# ==============================================================================
# 12_0_landcover_dominant_class.R — Dominant 0.25-degree CCI land-cover class
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "landcover_classification.R"))

terraOptions(progress = 1, memfrac = 0.2, todisk = TRUE)

years <- 1992:2022
fraction_dir <- here(
  "analysis", "tmp", paste0("lc025_fraction_yearly_", lc_fraction_schema)
)
out_dir <- here(
  "analysis", "tmp", paste0("lc025_dominant_class_grouped19_", lc_fraction_schema)
)
out_file <- file.path(out_dir, "lc025_dominant_mean_1992-2022.tif")
files <- file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", years))

lc_assert_fraction_files(files)

if (file.exists(out_file) && file.mtime(out_file) >= max(file.mtime(files))) {
  message("Dominant land-cover classification is current; skipping calculation.")
  quit(save = "no", status = 0L)
}

group_names <- paste0("lc_", names(lc_group_members))
mean_fraction <- rast(lapply(lc_group_members, function(class_ids) {
  annual_group <- rast(lapply(files, function(path) {
    layers <- rast(path)[[paste0("lc_", class_ids)]]
    if (nlyr(layers) == 1L) {
      return(layers)
    }
    app(layers, sum, na.rm = FALSE)
  }))
  app(annual_group, mean, na.rm = TRUE)
}))
names(mean_fraction) <- group_names

dominant_index <- which.max(mean_fraction)
class_ids <- as.integer(names(lc_group_members))
dominant_class <- classify(
  dominant_index,
  cbind(seq_along(class_ids), class_ids),
  others = NA
)
names(dominant_class) <- "lc_id"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(out_file, ".tmp.tif")
writeRaster(
  dominant_class,
  tmp,
  overwrite = TRUE,
  datatype = "INT2S",
  NAflag = -9999,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER")
)
if (!file.rename(tmp, out_file)) {
  stop("Could not replace ", out_file, call. = FALSE)
}
message("Wrote dominant land-cover classification: ", out_file)
