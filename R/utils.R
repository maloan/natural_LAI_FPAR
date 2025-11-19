## =============================================================================
# utils.R — General-purpose helpers for config, timing, GDAL, and geometry
#
# Purpose
#   Shared low-level tools for configuration handling, timing, GDAL options,
#   raster alignment, and lightweight parallel execution.
#
# Dependencies
#   terra, yaml, parallel, here
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(parallel)
  library(here)
})

# ==============================================================================
# Basic utilities
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

tok <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))

exp_ <- function(p) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)

find_one <- function(dir, pattern) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(cand)) {
    stop("No file matching pattern '", pattern, "' in: ", dir, call. = FALSE)
  }
  cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
}

cfg_read <- function() {
  yaml::read_yaml(here("config", "config.yml"))
}


# ==============================================================================
# Timing helper
# ==============================================================================

.with_timing <- function(label,
                         expr,
                         csv = NULL,
                         enabled = TRUE,
                         script_name = NULL) {

  if (!enabled) return(force(expr))

  t0 <- unclass(Sys.time())
  message(sprintf("[TIC] %s", label))

  on.exit({
    dt <- unclass(Sys.time()) - t0
    cat(sprintf("[TOC] %-28s %.2f s\n", label, dt))
    if (!is.null(csv)) {
      row <- data.frame(
        script  = script_name %||% label,
        step    = label,
        seconds = dt,
        time    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      if (!file.exists(csv)) {
        write.csv(row, csv, row.names = FALSE)
      } else {
        write.table(
          row,
          csv,
          sep = ",",
          col.names = FALSE,
          row.names = FALSE,
          append = TRUE
        )
      }
    }
  }, add = TRUE)

  force(expr)
}


# ==============================================================================
# GDAL write options
# ==============================================================================

gdal_wopt <- function(dtype = "FLT4S",
                      compress = "DEFLATE",
                      threads = "AUTO") {

  thr <- if (identical(threads, "AUTO")) {
    "NUM_THREADS=ALL_CPUS"
  } else {
    sprintf("NUM_THREADS=%s", threads)
  }

  list(
    datatype = dtype,
    NAflag   = if (dtype %in% c("INT1U","INT2U","INT2S","INT4U","INT4S")) 255 else NULL,
    gdal = c(
      sprintf("COMPRESS=%s", compress),
      "TILED=YES",
      "SPARSE_OK=TRUE",
      if (dtype == "FLT4S") "PREDICTOR=3" else "PREDICTOR=2",
      thr,
      "BIGTIFF=IF_SAFER"
    )
  )
}


# ==============================================================================
# Geometry helpers
# ==============================================================================

align_to <- function(r, ref, method = c("near", "bilinear")) {
  method <- match.arg(method)
  if (compareGeom(r, ref, stopOnError = FALSE)) r else resample(r, ref, method)
}

.is_int_like <- function(x, tol = 1e-10) abs(x - round(x)) < tol

.is_snapped <- function(r, origin_x = -180, origin_y = -90, tol = 1e-10) {
  rx <- res(r)[1]
  ry <- res(r)[2]
  .is_int_like((xmin(r) - origin_x) / rx, tol) &&
    .is_int_like((ymin(r) - origin_y) / ry, tol)
}

.snap_to_grid <- function(r, origin_x = -180, origin_y = -90) {
  rx <- res(r)[1]
  ry <- res(r)[2]
  nx <- ncol(r)
  ny <- nrow(r)

  x0 <- origin_x + round((xmin(r) - origin_x) / rx) * rx
  y0 <- origin_y + round((ymin(r) - origin_y) / ry) * ry

  ext(r) <- ext(x0, x0 + nx * rx, y0, y0 + ny * ry)
  r
}


# ==============================================================================
# Environment handling
# ==============================================================================

env_get_int <- function(key, def) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) return(def)

  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) def else out
}


# ==============================================================================
# Parallel helpers
# ==============================================================================

detect_workers <- function() {
  max(1L, min(3L, max(1L, parallel::detectCores(TRUE) - 1L)))
}

runner <- function(tasks, FUN, n_workers = detect_workers()) {
  if (.Platform$OS.type == "windows" || n_workers == 1L) {
    lapply(tasks, FUN)
  } else {
    parallel::mclapply(tasks, FUN, mc.cores = n_workers)
  }
}
