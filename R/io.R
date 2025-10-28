# io.R

# --- Provenance helper for GeoTIFF DESCRIPTION tag ---
.area_def_str <- function(res_deg) {
  sprintf("global lon[-180,180], lat[-90,90], res[%.2f°]", res_deg)
}

gtiff_desc <- function(product,
                       grid_deg,
                       extras = list(),
                       crs = "EPSG:4326") {
  base <- c(
    sprintf("PRODUCT=%s", product),
    sprintf("GRID=%.2f°", grid_deg),
    sprintf("CRS=%s", crs),
    sprintf("AREA_DEF=%s", .area_def_str(grid_deg))
  )
  if (length(extras)) {
    add <- mapply(
      function(k, v)
        sprintf("%s=%s", toupper(k), as.character(v)),
      names(extras),
      extras,
      USE.NAMES = FALSE
    )
    base <- c(base, add)
  }
  paste(base, collapse = "; ")
}

cfg_sha1 <- function(cfg_path) {
  raw <- readChar(cfg_path, file.info(cfg_path)$size, useBytes = TRUE)
  as.character(openssl::sha1(charToRaw(raw)))
}

wopt_byte <- function(speed, na = 255L) {
  list(datatype = "INT1U",
       NAflag = as.integer(na),
       gdal = gdal_co_int(speed))
}
wopt_int  <- function(speed) {
  list(datatype = "INT2U", gdal = gdal_co_int(speed))
}
wopt_f32  <- function(speed) {
  list(datatype = "FLT4S", gdal = gdal_co_f32(speed))
}

gdal_co_int <- function(speed_over_size = FALSE) {
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (speed_over_size)
  {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(base, comp, "PREDICTOR=2")
}
gdal_co_f32 <- function(speed_over_size = FALSE) {
  base <- c("TILED=YES", "BIGTIFF=IF_SAFER")
  comp <- if (speed_over_size) {
    "COMPRESS=LZW"
  } else {
    "COMPRESS=DEFLATE"
  }
  c(base, comp, "PREDICTOR=3")
}
gdal_threads <- function() {
  "NUM_THREADS=ALL_CPUS"
}


write_session_info <- function(log_dir, tag) {
  dir.create(log_dir, TRUE, showWarnings = FALSE)
  fn <- file.path(log_dir, sprintf("session_%s.txt", tag))
  sink(fn)
  on.exit(sink(), add = TRUE)
  print(sessionInfo())
  cat("\nSys.getenv snapshot:\n")
  print(as.list(Sys.getenv()))
}

# small manifest helper (csv per stage)
write_manifest <- function(out_csv, rasters) {
  rows <- lapply(rasters, function(p) {
    r <- terra::rast(p)
    rng <- try(terra::global(r, fun = "range", na.rm = TRUE)[1, ], silent = TRUE)
    nafrac <- try({
      m <- terra::ifel(is.finite(r), 1, NA)
      1 - terra::global(m, "mean", na.rm = TRUE)[1, 1]
    }, silent = TRUE)
    data.frame(
      file = basename(p),
      full_path = normalizePath(p, winslash = "/", mustWork = FALSE),
      rows = terra::nrow(r),
      cols = terra::ncol(r),
      res_x = terra::res(r)[1],
      res_y = terra::res(r)[2],
      xmin = terra::xmin(r),
      xmax = terra::xmax(r),
      ymin = terra::ymin(r),
      ymax = terra::ymax(r),
      crs  = as.character(terra::crs(r)),
      vmin = suppressWarnings(ifelse(is.list(rng), NA, rng[1])),
      vmax = suppressWarnings(ifelse(is.list(rng), NA, rng[2])),
      na_fraction = suppressWarnings(as.numeric(nafrac)),
      stage = sub("^manifest_(.+)\\.csv$", "\\1", basename(out_csv)),
      generated_by = Sys.getenv("PIPE_SCRIPT", "<unknown>"),
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  dir.create(dirname(out_csv), TRUE, showWarnings = FALSE)
  utils::write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
}

mask_legend <- function(where = "bottomleft") {
  legend(where, fill = c("#f0f0f0", "#d73027"), legend = c("0 keep", "1 drop"), bty = "n")
}
