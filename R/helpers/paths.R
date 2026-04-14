## =============================================================================
## paths.R — Path and config helpers
## =============================================================================

exp_ <- function(p) {
  normalizePath(path.expand(p),
    winslash = "/",
    mustWork = FALSE
  )
}

cfg_read <- function() {
  root <- exp_(Sys.getenv(
    "SNU_LAI_FPAR_ROOT",
    unset = "~/GitHub/natural_LAI_FPAR"
  ))
  yaml::read_yaml(file.path(root, "config", "config.yml"))
}
