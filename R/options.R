## =============================================================================
# options.R — Environment-variable parsing and runtime option handling
#
# Purpose
#   Provide consistent parsing of environment variables controlling pipeline
#   behavior (e.g., overwrite, quicklooks, concurrency, compression modes).
#
# Functions
#   as_bool()       — Convert character input (e.g., "TRUE","1","yes") to logical
#   env_get_int()   — Retrieve integer environment variable with fallback default
#   opts_read()     — Assemble a list of runtime options for scripts
#
# Env vars
#   FORCE_REBUILD / FORCE, REMAKE_QL, SKIP_EXISTING, SPEED_OVER_SIZE,
#   USE_GDAL_MODE, N_WORKERS
## =============================================================================

suppressPackageStartupMessages({
  # parallel is only needed for detectCores; safe to require lazily
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' not available; N_WORKERS will default to 1.")
  }
})

as_bool <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return(FALSE)
  }
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y", "on")
}

env_get_int <- function(key, def) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) {
    return(as.integer(def))
  }
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) {
    as.integer(def)
  } else {
    as.integer(out)
  }
}

opts_read <- function() {
  # Detect core count robustly
  cores <- 1L
  if ("parallel" %in% .packages(all.available = TRUE)) {
    cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
    if (!is.finite(cores) || cores < 1L) {
      cores <- 1L
    }
  }

  # Defaults
  d_FORCE          <- FALSE
  d_REMAKE_QL      <- FALSE
  d_SKIP_EXISTING  <- TRUE
  d_SPEED_OVER_SIZE<- FALSE
  d_USE_GDAL_MODE  <- TRUE
  d_N_WORKERS      <- max(1L, cores - 1L)

  # Parse
  FORCE <- as_bool(Sys.getenv("FORCE_REBUILD", Sys.getenv("FORCE", ifelse(d_FORCE, "1", "0"))))
  REMAKE_QL <- as_bool(Sys.getenv("REMAKE_QL", ifelse(d_REMAKE_QL, "1", "0")))
  SKIP_EXISTING <- as_bool(Sys.getenv("SKIP_EXISTING", ifelse(d_SKIP_EXISTING, "1", "0")))
  SPEED_OVER_SIZE <- as_bool(Sys.getenv("SPEED_OVER_SIZE", ifelse(d_SPEED_OVER_SIZE, "1", "0")))
  USE_GDAL_MODE <- as_bool(Sys.getenv("USE_GDAL_MODE", ifelse(d_USE_GDAL_MODE, "1", "0")))
  N_WORKERS <- env_get_int("N_WORKERS", d_N_WORKERS)

  # Clamp workers to [1, cores]
  if (is.na(N_WORKERS) || N_WORKERS < 1L) {
    N_WORKERS <- 1L
  }
  if (is.finite(cores) && N_WORKERS > cores) {
    N_WORKERS <- cores
  }

  list(
    FORCE = FORCE,
    REMAKE_QL = REMAKE_QL,
    SKIP_EXISTING = SKIP_EXISTING,
    SPEED_OVER_SIZE = SPEED_OVER_SIZE,
    USE_GDAL_MODE = USE_GDAL_MODE,
    N_WORKERS = N_WORKERS
  )
}
