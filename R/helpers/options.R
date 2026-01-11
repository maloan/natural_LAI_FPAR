## =============================================================================
## options.R — Environment-variable parsing and runtime option handling
## =============================================================================

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x)) {
    return(default)
  }
  x <- tolower(trimws(as.character(x)))
  if (!nzchar(x)) {
    return(default)
  }
  x %in% c("1", "true", "t", "yes", "y", "on")
}

env_get_int <- function(key, default = NA_integer_) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) {
    return(as.integer(default))
  }
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) as.integer(default) else out
}

detect_cores <- function() {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    return(1L)
  }
  cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) 1L
  )
  if (!is.finite(cores) || cores < 1L) 1L else as.integer(cores)
}

opts_read <- function() {
  cores <- detect_cores()

  FORCE <- as_bool(Sys.getenv(
    "FORCE_REBUILD",
    Sys.getenv("FORCE", "")
  ), default = FALSE)
  REMAKE_QL <- as_bool(Sys.getenv("REMAKE_QL", ""), default = FALSE)
  SKIP_EXISTING <- as_bool(Sys.getenv("SKIP_EXISTING", ""), default = TRUE)
  SPEED_OVER_SIZE <- as_bool(
    Sys.getenv("SPEED_OVER_SIZE", ""),
    default = FALSE
  )
  USE_GDAL_MODE <- as_bool(Sys.getenv("USE_GDAL_MODE", ""), default = TRUE)

  # default workers = cores - 1 (leave one core for OS), but at least 1
  default_workers <- max(1L, cores - 1L)
  N_WORKERS <- env_get_int("N_WORKERS", default_workers)
  if (is.na(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L
  if (N_WORKERS > cores) N_WORKERS <- cores

  list(
    FORCE = FORCE,
    REMAKE_QL = REMAKE_QL,
    SKIP_EXISTING = SKIP_EXISTING,
    SPEED_OVER_SIZE = SPEED_OVER_SIZE,
    USE_GDAL_MODE = USE_GDAL_MODE,
    N_WORKERS = N_WORKERS,
    N_CORES = cores
  )
}
