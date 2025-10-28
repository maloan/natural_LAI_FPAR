# ref/00_utils.R
suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(parallel)
})

# ---------- small helpers ----------
`%||%` <- function(a, b) if (is.null(a)) b else a
exp_ <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)

cfg_read <- function() {
  ROOT <- exp_(Sys.getenv("SNU_LAI_FPAR_ROOT", unset = "~/GitHub/natural_LAI_FPAR"))
  yaml::read_yaml(file.path(ROOT, "config", "config.yml"))
}

# ref/00_utils.R  — small, optional improvements

.with_timing <- function(label, expr, csv = NULL, enabled = TRUE, script_name = NULL) {
  if (!enabled) return(force(expr))
  t0 <- unclass(Sys.time())
  message(sprintf("[TIC] %s", label))
  on.exit({
    dt <- unclass(Sys.time()) - t0
    cat(sprintf("[TOC] %-28s %.2f s\n", label, dt))
    if (!is.null(csv)) {
      row <- data.frame(
        script  = if (is.null(script_name)) label else script_name,
        step    = label,
        seconds = dt,
        time    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      if (!file.exists(csv)) write.csv(row, csv, row.names = FALSE) else
        write.table(row, csv, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
  }, add = TRUE)
  force(expr)
}

gdal_wopt <- function(dtype = "FLT4S", compress = "DEFLATE", threads = "AUTO") {
  thr <- if (identical(threads, "AUTO")) "NUM_THREADS=ALL_CPUS" else sprintf("NUM_THREADS=%s", threads)
  list(
    datatype = dtype,
    # keep NAflag only for integer outputs; floats in terra generally don’t need a custom NAflag
    NAflag   = if (dtype %in% c("INT1U","INT2U","INT2S","INT4U","INT4S")) 255 else NULL,
    gdal     = c(
      sprintf("COMPRESS=%s", compress), "TILED=YES", "SPARSE_OK=TRUE",
      if (dtype == "FLT4S") "PREDICTOR=3" else "PREDICTOR=2",
      thr, "BIGTIFF=IF_SAFER"
    )
  )
}


align_to <- function(r, ref, method = c("near","bilinear")) {
  method <- match.arg(method)
  if (compareGeom(r, ref, stopOnError = FALSE)) return(r)
  resample(r, ref, method)
}


.is_snapped <- function(r, origin_x = -180, origin_y = -90, tol = 1e-10) {
  rx <- res(r)[1]; ry <- res(r)[2]
  qx <- (xmin(r) - origin_x) / rx
  qy <- (ymin(r) - origin_y) / ry
  .is_int_like(qx, tol) && .is_int_like(qy, tol)
}

.snap_to_grid <- function(r, origin_x = -180, origin_y = -90) {
  rx <- res(r)[1]; ry <- res(r)[2]
  nx <- ncol(r);   ny <- nrow(r)
  x0 <- origin_x + round((xmin(r) - origin_x) / rx) * rx
  y0 <- origin_y + round((ymin(r) - origin_y) / ry) * ry
  ext(r) <- ext(x0, x0 + nx * rx, y0, y0 + ny * ry)
  r
}



env_get_int <- function(key, def) {
  v <- Sys.getenv(key, unset = "")
  if (identical(v, "")) return(def)
  vv <- suppressWarnings(as.integer(v))
  if (is.na(vv)) def else vv
}

detect_workers <- function() max(1L, min(3L, max(1L, parallel::detectCores(TRUE) - 1L)))

runner <- function(tasks, FUN, n_workers = detect_workers()) {
  if (.Platform$OS.type == "windows" || n_workers == 1L) lapply(tasks, FUN) else
    parallel::mclapply(tasks, FUN, mc.cores = n_workers)
}
