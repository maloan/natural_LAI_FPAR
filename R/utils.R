## =============================================================================
# utils.R — General-purpose helpers for config, timing, GDAL, and geometry
#
# Purpose
#   Define shared low-level tools for configuration handling, timing, raster
#   alignment, GDAL write settings, and lightweight parallel execution.
#
# Functions
#   %||%()              — Null-coalescing operator
#   find_one()          — Find newest matching file in directory
#   tok()               — Encode numeric thresholds (e.g., 0.05 → "0p05")
#   exp_()              — Expand and normalize file paths
#   cfg_read()          — Load global config.yml
#   .with_timing()      — Time and log script steps to CSV
#   gdal_wopt()         — Standardized GDAL write options
#   align_to()          — Resample raster to reference grid
#   .is_snapped()       — Check raster alignment to lon/lat grid origin
#   .snap_to_grid()     — Correct minor grid misalignments
#   env_get_int()       — Read integer environment variable with fallback
#   detect_workers()    — Detect usable CPU cores
#   runner()            — Apply tasks serially or in parallel
#
# Inputs
#   - Environment variables and raster objects as needed by each helper.
#
# Outputs
#   - Utility return values, no direct file outputs.
#
# Dependencies
#   Packages: terra, yaml, parallel
#
# Processing overview
#   1) Configure and load workflow metadata.
#   2) Support standardized GDAL options and geometry checks.
#   3) Provide safe, cross-platform timing and parallel helpers.
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
  library(parallel)
})

`%||%` <- function(a, b) {
  if (is.null(a)) {
    return(b)
  } else {
    return(a)
  }
}

find_one <- function(dir, pattern) {
  cand <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(cand)) {
    stop("No file matching pattern '", pattern, "' in: ", dir, call. = FALSE)
  } else {
    cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
  }
}

tok <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

exp_ <- function(p) {
  normalizePath(path.expand(p), winslash = "/", mustWork = FALSE)
}

cfg_read <- function() {
  ROOT <- exp_(Sys.getenv("SNU_LAI_FPAR_ROOT", unset = "~/GitHub/natural_LAI_FPAR"))
  yaml::read_yaml(file.path(ROOT, "config", "config.yml"))
}

.with_timing <- function(label,
                         expr,
                         csv = NULL,
                         enabled = TRUE,
                         script_name = NULL) {
  if (!enabled) {
    return(force(expr))
  }

  t0 <- unclass(Sys.time())
  message(sprintf("[TIC] %s", label))

  on.exit({
    dt <- unclass(Sys.time()) - t0
    cat(sprintf("[TOC] %-28s %.2f s\n", label, dt))
    if (!is.null(csv)) {
      row <- data.frame(
        script  = if (is.null(script_name))
          label
        else
          script_name,
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
    NAflag   = if (dtype %in% c("INT1U", "INT2U", "INT2S", "INT4U", "INT4S")) {
      255
    } else {
      NULL
    },
    gdal     = c(
      sprintf("COMPRESS=%s", compress),
      "TILED=YES",
      "SPARSE_OK=TRUE",
      if (dtype == "FLT4S") {
        "PREDICTOR=3"
      } else {
        "PREDICTOR=2"
      },
      thr,
      "BIGTIFF=IF_SAFER"
    )
  )
}

align_to <- function(r, ref, method = c("near", "bilinear")) {
  method <- match.arg(method)
  if (compareGeom(r, ref, stopOnError = FALSE)) {
    return(r)
  } else {
    return(resample(r, ref, method))
  }
}

.is_snapped <- function(r,
                        origin_x = -180,
                        origin_y = -90,
                        tol = 1e-10) {
  rx <- res(r)[1]
  ry <- res(r)[2]
  qx <- (xmin(r) - origin_x) / rx
  qy <- (ymin(r) - origin_y) / ry
  .is_int_like(qx, tol) && .is_int_like(qy, tol)
}

.snap_to_grid <- function(r,
                          origin_x = -180,
                          origin_y = -90) {
  rx <- res(r)[1]
  ry <- res(r)[2]
  nx <- ncol(r)
  ny <- nrow(r)
  x0 <- origin_x + round((xmin(r) - origin_x) / rx) * rx
  y0 <- origin_y + round((ymin(r) - origin_y) / ry) * ry
  ext(r) <- ext(x0, x0 + nx * rx, y0, y0 + ny * ry)
  r
}

env_get_int <- function(key, def) {
  v <- Sys.getenv(key, unset = "")
  if (identical(v, "")) {
    return(def)
  } else {
    vv <- suppressWarnings(as.integer(v))
    if (is.na(vv)) {
      return(def)
    } else {
      return(vv)
    }
  }
}

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
