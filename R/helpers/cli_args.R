# =============================================================================
# cli_args.R — Parse command-line arguments in key=value format
# =============================================================================

#' Parse and validate command-line arguments
#'
#' Converts command-line arguments in key=value format to a named list,
#' validating against a provided defaults list. Applies type conversions
#' based on known parameter names (logical, integer, numeric).
#'
#' @param defaults Named list with default values and target types.
#'                 Keys must match argument names; values define types.
#'
#' @return Updated list with command-line values overriding defaults.
#'
#' @examples
#' defaults <- list(var = "LAI", mask = "CCI", overwrite = FALSE)
#' cfg <- parse_cli_args(defaults) # Called with: Rscript script.R var=FPAR overwrite=TRUE
#'
parse_cli_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    return(defaults)
  }

  kv <- strsplit(args, "=", fixed = TRUE)
  ok <- lengths(kv) == 2
  if (!all(ok)) {
    bad <- args[!ok]
    stop("Invalid argument(s); use key=value format: ", paste(bad, collapse = ", "))
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
