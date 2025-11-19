## =============================================================================
# options.R — Environment-variable parsing and runtime option handling
#
# Purpose
#   Unified parsing of environment variables controlling overwrite behavior,
#   quicklook generation, compression mode, and parallelism.
#
# Functions
#   as_bool()       — Robust logical conversion
#   env_get_int()   — Integer env var with fallback
#   opts_read()     — Consolidated runtime options
#
# Env vars used
#   FORCE_REBUILD / FORCE
#   REMAKE_QL
#   SKIP_EXISTING
#   SPEED_OVER_SIZE
#   USE_GDAL_MODE
#   N_WORKERS
## =============================================================================


# ------------------------------------------------------------------------------
# Basic converters
# ------------------------------------------------------------------------------

as_bool <- function(x) {
  if (length(x) == 0L || is.na(x)) return(FALSE)
  tolower(trimws(as.character(x))) %in% c("1","true","t","yes","y","on")
}

env_get_int <- function(key, def) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) return(as.integer(def))

  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) as.integer(def) else as.integer(out)
}


# ------------------------------------------------------------------------------
# Runtime-option collector
# ------------------------------------------------------------------------------

opts_read <- function() {

  # number of cores
  cores <- 1L
  if ("parallel" %in% .packages(all.available = TRUE)) {
    cores <- tryCatch(parallel::detectCores(logical = TRUE),
                      error = function(e) 1L)
    if (!is.finite(cores) || cores < 1L) cores <- 1L
  }

  # defaults
  d <- list(
    FORCE           = FALSE,
    REMAKE_QL       = FALSE,
    SKIP_EXISTING   = TRUE,
    SPEED_OVER_SIZE = FALSE,
    USE_GDAL_MODE   = TRUE,
    N_WORKERS       = max(1L, cores - 1L)
  )

  # parse with fallbacks
  FORCE <- as_bool(Sys.getenv("FORCE_REBUILD",
                              Sys.getenv("FORCE", ifelse(d$FORCE, "1","0"))))
  REMAKE_QL     <- as_bool(Sys.getenv("REMAKE_QL",     ifelse(d$REMAKE_QL, "1","0")))
  SKIP_EXISTING <- as_bool(Sys.getenv("SKIP_EXISTING", ifelse(d$SKIP_EXISTING, "1","0")))
  SPEED_OVER_SIZE <- as_bool(Sys.getenv("SPEED_OVER_SIZE",
                                        ifelse(d$SPEED_OVER_SIZE,"1","0")))
  USE_GDAL_MODE <- as_bool(Sys.getenv("USE_GDAL_MODE",
                                      ifelse(d$USE_GDAL_MODE,"1","0")))
  N_WORKERS <- env_get_int("N_WORKERS", d$N_WORKERS)

  # clamp parallelism
  if (is.na(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L
  if (is.finite(cores) && N_WORKERS > cores) N_WORKERS <- cores

  list(
    FORCE           = FORCE,
    REMAKE_QL       = REMAKE_QL,
    SKIP_EXISTING   = SKIP_EXISTING,
    SPEED_OVER_SIZE = SPEED_OVER_SIZE,
    USE_GDAL_MODE   = USE_GDAL_MODE,
    N_WORKERS       = N_WORKERS
  )
}
