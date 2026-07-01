# =============================================================================
# options.R — Helpers for reading and managing options from environment
# variables
# =============================================================================

as_bool <- function(x, default = FALSE) {
  # Convert a string to a boolean value, returning a default if the input is
  # empty or NA.
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
  # Retrieve an integer value from an environment variable, returning a default
  # if the variable is not set or cannot be converted to an integer.
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) {
    return(as.integer(default))
  }
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) {
    as.integer(default)
  } else {
    out
  }
  # Note: implicit return of last expression
}


opts_read <- function() {
  # Read and return a list of options from environment variables, with defaults
  # for each option.
  cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) {
      1L
    }
  )
  cores <- as.integer(max(1L, cores))
  force <- as_bool(Sys.getenv("FORCE_REBUILD", Sys.getenv("force", "")), default = FALSE)
  remake_ql <- as_bool(Sys.getenv("remake_ql", ""), default = FALSE)
  skip_existing <- as_bool(Sys.getenv("skip_existing", ""), default = TRUE)
  speed_over_size <- as_bool(Sys.getenv("speed_over_size", ""), default = FALSE)
  use_gdal_mode <- as_bool(Sys.getenv("use_gdal_mode", ""), default = TRUE)

  # default workers = cores - 1 (leave one core for OS), but at least 1
  default_workers <- max(1L, cores - 1L)
  n_workers <- env_get_int("n_workers", default_workers)
  if (is.na(n_workers) || n_workers < 1L) {
    n_workers <- 1L
  }
  if (n_workers > cores) {
    n_workers <- cores
  }

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
