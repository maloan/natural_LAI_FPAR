# ==============================================================================
# 12_0_landcover_dominant_class.R — Dominant 0.25-degree CCI land-cover class
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

years <- 1992:2022
fraction_dir <- here("analysis", "tmp", "lc025_fraction_yearly")
out_dir <- here("analysis", "tmp", "lc025_dominant_class")
out_file <- file.path(out_dir, "lc025_dominant_mean_1992-2022.tif")
files <- file.path(fraction_dir, sprintf("lc025_fraction_%d.tif", years))

missing_files <- files[!file.exists(files)]
if (length(missing_files)) {
  stop(
    "Missing annual 0.25-degree land-cover fractions. Run ",
    "R/12_make_lc025_fractions.R first:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (file.exists(out_file) && file.mtime(out_file) >= max(file.mtime(files))) {
  message("Dominant land-cover classification is current; skipping calculation.")
  quit(save = "no", status = 0L)
}

annual <- rast(files)
class_names <- names(rast(files[[1]]))
n_classes <- length(class_names)
mean_fraction <- rast(lapply(seq_len(n_classes), function(i) {
  app(annual[[seq(i, nlyr(annual), by = n_classes)]], mean, na.rm = TRUE)
}))
names(mean_fraction) <- class_names

dominant_index <- which.max(mean_fraction)
class_ids <- as.integer(sub("^lc_", "", class_names))
dominant_class <- classify(
  dominant_index,
  cbind(seq_along(class_ids), class_ids),
  others = NA
)
names(dominant_class) <- "lc_id"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeRaster(
  dominant_class,
  out_file,
  overwrite = TRUE,
  datatype = "INT2S",
  NAflag = -9999,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER")
)
message("Wrote dominant land-cover classification: ", out_file)
