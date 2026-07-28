# =============================================================================
# cli_args.R — Command-line argument parsing and scenario specification
# =============================================================================

parse_cli_args <- function(defaults) {
  # Parse command-line arguments in key=value format and return a named list of
  # values.
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    return(defaults)
  }

  kv <- strsplit(args, "=", fixed = TRUE)
  ok <- lengths(kv) == 2
  if (!all(ok)) {
    bad <- args[!ok]
    stop(
      "Invalid argument(s); use key=value format: ",
      paste(bad, collapse = ", ")
    )
  }

  vals <- defaults
  for (pair in kv) {
    k <- pair[[1]]
    v <- pair[[2]]
    if (!nzchar(v)) {
      next
    }
    if (!k %in% names(vals)) {
      warning("Unknown argument ignored: ", k, call. = FALSE)
      next
    }
    vals[[k]] <- v
  }

  vals$use_relative <- tolower(as.character(vals$use_relative)) %in% c("true", "1", "yes", "y")
  vals$make_preview_map <- tolower(as.character(vals$make_preview_map)) %in% c("true", "1", "yes", "y")
  vals$lc_year_start <- suppressWarnings(as.integer(vals$lc_year_start))
  vals$lc_year_end <- suppressWarnings(as.integer(vals$lc_year_end))
  vals$stability_min <- suppressWarnings(as.numeric(vals$stability_min))
  vals$chunk_size <- suppressWarnings(as.integer(vals$chunk_size))
  vals$top_n <- suppressWarnings(as.numeric(vals$top_n))
  vals$top_n3 <- suppressWarnings(as.numeric(vals$top_n3))
  vals
}

create_scenario_spec <- function(cci_alphas = c("alpha_0.05", "alpha_0.1", "alpha_0.2"),
                                 glc_run_tag = "alpha_0.1") {
  # Create a tibble of scenario specifications for the analysis.
  tibble::tibble(
    scenario = c("Unmasked", sprintf("CCI alpha=%s", sub(
      "alpha_", "", cci_alphas
    )), "GLC"),
    source = c("unmasked", rep("CCI", length(cci_alphas)), "GLC"),
    run_tag = c(NA_character_, cci_alphas, glc_run_tag)
  )
}

scenario_order <- function(scenario_spec) {
  # Return the order of scenarios for plotting and analysis.
  as.character(scenario_spec$scenario)
}
