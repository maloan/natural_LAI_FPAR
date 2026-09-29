# ==============================================================================
# 00_area_validdomain_after_nonvegetated.R — Create base and scenario-specific
# retained-area rasters (0.05° and 0.25°)
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(tibble)
  library(readr)
  library(here)
})

source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "cli_args.R"))

terraOptions(progress = 0, memfrac = 0.4)

outdir <- here("analysis", "results", "tables", "masks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_cdf_if_stale <- function(x, path, dependencies) {
  dependencies <- dependencies[file.exists(dependencies)]
  is_current <- file.exists(path) &&
    length(dependencies) > 0 &&
    file.info(path)$mtime >= max(file.info(dependencies)$mtime)

  if (!is_current) {
    writeCDF(x, path, overwrite = TRUE)
  }
  invisible(path)
}

area_005 <- here("src", "area_0p05_km2.nc")
mask_nonveg <- here(
  "output",
  "alpha_0.1",
  "masks",
  "mask_nonvegetated",
  "mask_nonvegetated_CCI_2007_alphaW0p05_alphaI0p05_0p05.tif"
)

area <- rast(area_005)[[1]]
m_nonveg <- rast(mask_nonveg)[[1]]
stopifnot(compareGeom(m_nonveg, area, stopOnError = TRUE))

area_sum <- function(cond) {
  global(ifel(cond, area, NA), "sum", na.rm = TRUE)[1, 1] |> as.numeric()
}
support_dom <- is.finite(area) & (area > 0)
nonveg_excl <- support_dom & (m_nonveg == 1)
land_dom <- support_dom & !nonveg_excl & is.finite(m_nonveg)

area_support <- area_sum(support_dom)
area_nonveg_excl <- area_sum(support_dom & nonveg_excl)
area_land_after_nonveg <- area_sum(land_dom)

tbl_nonveg <- tibble(
  denom_support_km2 = area_support,
  nonvegetated_removed_km2 = area_nonveg_excl,
  nonvegetated_removed_pct_of_support = 100 * (area_nonveg_excl / area_support),
  land_after_nonvegetated_km2 = area_land_after_nonveg,
  land_after_nonvegetated_pct_of_support = 100 * (area_land_after_nonveg / area_support)
)

area_valid_005 <- area
area_valid_005[!land_dom] <- NA
area_valid_025 <- aggregate(area_valid_005,
  fact = 5,
  fun = "sum",
  na.rm = TRUE
)
area_valid_025[area_valid_025 == 0] <- NA

# Write summary and refresh rasters only when their inputs are newer.
write_csv(
  round_numeric(tbl_nonveg, 5),
  file.path(outdir, "domain_nonvegetated_0p05.csv")
)
area_valid_005_file <- here("src", "area_0p05_validdomain_km2.nc")
area_valid_025_file <- here("src", "area_0p25_validdomain_km2.nc")
write_cdf_if_stale(
  area_valid_005,
  area_valid_005_file,
  c(area_005, mask_nonveg)
)
write_cdf_if_stale(
  area_valid_025,
  area_valid_025_file,
  c(area_005, mask_nonveg)
)

# Aggregate the exact 0.05-degree area retained by each managed-land mask.
scenario_spec <- create_scenario_spec(
  cci_alphas = c("alpha_0.05", "alpha_0.1", "alpha_0.2"),
  glc_run_tag = "alpha_0.1"
)
scenario_spec <- scenario_spec[scenario_spec$source != "unmasked", ]
area_025_vals <- values(area_valid_025, dataframe = FALSE)
base_valid <- is.finite(area_025_vals) & area_025_vals > 0
base_area_total <- sum(area_025_vals[base_valid])
qa_rows <- vector("list", nrow(scenario_spec))

for (i in seq_len(nrow(scenario_spec))) {
  sc <- scenario_spec[i, ]
  mask_file <- combined_mask_path(sc$source, sc$run_tag)
  if (!file.exists(mask_file)) {
    stop("Missing combined mask: ", mask_file)
  }
  drop_mask <- rast(mask_file)[[1]]
  stopifnot(compareGeom(drop_mask, area_valid_005, stopOnError = TRUE))

  retained_005 <- ifel(
    is.finite(area_valid_005) &
      area_valid_005 > 0 &
      (is.na(drop_mask) | drop_mask != 1),
    area_valid_005,
    NA
  )
  retained_025 <- aggregate(
    retained_005,
    fact = 5,
    fun = "sum",
    na.rm = TRUE
  )
  retained_025[retained_025 <= 0] <- NA
  retained_025 <- ifel(
    is.finite(retained_025) & is.finite(area_valid_025),
    ifel(
      retained_025 > area_valid_025,
      area_valid_025,
      retained_025
    ),
    NA
  )
  stopifnot(compareGeom(retained_025, area_valid_025, stopOnError = TRUE))

  output_file <- scenario_area_path(sc$source, sc$run_tag)
  write_cdf_if_stale(
    retained_025,
    output_file,
    c(area_valid_005_file, mask_file)
  )

  retained_vals <- values(retained_025, dataframe = FALSE)
  retained <- base_valid & is.finite(retained_vals) & retained_vals > 0
  retained_fraction <- numeric(length(area_025_vals))
  retained_fraction[retained] <- retained_vals[retained] /
    area_025_vals[retained]
  partial <- retained & retained_fraction < (1 - 1e-6)
  retained_area <- sum(retained_vals[retained])
  full_cell_area <- sum(area_025_vals[retained])

  qa_rows[[i]] <- tibble(
    scenario = sc$scenario,
    run_tag = sc$run_tag,
    retained_area_km2 = retained_area,
    excluded_area_km2 = base_area_total - retained_area,
    retained_cells = sum(retained),
    partially_retained_cells = sum(partial),
    partially_retained_cells_pct = 100 * sum(partial) / sum(retained),
    retained_area_in_partial_cells_pct = 100 *
      sum(retained_vals[partial]) / retained_area,
    old_full_cell_area_assigned_km2 = full_cell_area,
    old_to_correct_area_ratio = full_cell_area / retained_area,
    output_file = output_file
  )
}

qa <- dplyr::bind_rows(qa_rows)
write_csv(
  round_numeric(qa, 6),
  file.path(outdir, "scenario_retained_area_0p25_qa.csv")
)
print(qa)
