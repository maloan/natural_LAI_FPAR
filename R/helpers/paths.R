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
  run_tag <- Sys.getenv("run_tag", "tau_0.1")
  cfg_file <- file.path(root, "config", sprintf("config_%s.yml", run_tag))

  # Fall back to generic config if run_tag specific one doesn't exist
  if (!file.exists(cfg_file)) {
    cfg_file <- file.path(root, "config", "config.yml")
  }

  yaml::read_yaml(cfg_file)
}
