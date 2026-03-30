

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
  # Note: implicit return of last expression
}

detect_cores <- function() {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    return(1L)
  }
  cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) 1L
  )
  as.integer(max(1L, cores))
}

opts_read <- function() {
  cores <- detect_cores()

  force <- as_bool(Sys.getenv(
    "FORCE_REBUILD",
    Sys.getenv("force", "")
  ), default = FALSE)
  remake_ql <- as_bool(Sys.getenv("remake_ql", ""), default = FALSE)
  skip_existing <- as_bool(Sys.getenv("skip_existing", ""), default = TRUE)
  speed_over_size <- as_bool(
    Sys.getenv("speed_over_size", ""),
    default = FALSE
  )
  use_gdal_mode <- as_bool(Sys.getenv("use_gdal_mode", ""), default = TRUE)

  # default workers = cores - 1 (leave one core for OS), but at least 1
  default_workers <- max(1L, cores - 1L)
  n_workers <- env_get_int("n_workers", default_workers)
  if (is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  if (n_workers > cores) n_workers <- cores

  list(
    force = force,
    remake_ql = remake_ql,
    skip_existing = skip_existing,
    speed_over_size = speed_over_size,
    use_gdal_mode = use_gdal_mode,
    n_workers = n_workers,
    n_cores = cores
  )
}