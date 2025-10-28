# ===============================
# File: ref/options.R
# ===============================

as_bool <- function(x) {
  tolower(trimws(x)) %in% c("1", "true", "t", "yes", "y", "on")
}
env_get_int <- function(k, def) {
  v <- Sys.getenv(k, "")
  if (!nzchar(v)) {
    def
  } else {
    x <- suppressWarnings(as.integer(v))
    if (is.na(x)) {
      def
    } else{
      x
    }
  }
}

opts_read <- function() {
  list(
    FORCE          = as_bool(Sys.getenv(
      "FORCE_REBUILD", Sys.getenv("FORCE", "0")
    )),
    REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "FALSE")),
    SKIP_EXISTING  = as_bool(Sys.getenv("SKIP_EXISTING", "1")),
    SPEED_OVER_SIZE = as_bool(Sys.getenv("SPEED_OVER_SIZE", "0")),
    USE_GDAL_MODE  = as_bool(Sys.getenv("USE_GDAL_MODE", "1")),
    N_WORKERS      = env_get_int("N_WORKERS", max(
      1L, parallel::detectCores(logical = TRUE) - 1L
    ))
  )
}