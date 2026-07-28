# =============================================================================
# 07_nonveg_snapshot_sensitivity.R
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "io.R"))
cfg <- cfg_read()

# locate files
run_tag <- cfg$project$run_tag
nonveg_dir <- here("output", run_tag, "masks", "mask_nonvegetated")

# Read all matching mask files
mask_files <- list.files(nonveg_dir, pattern = "^mask_nonvegetated_CCI_\\d{4}_alphaW0p05_alphaI0p05_0p05\\.tif$", full.names = TRUE)
if (length(mask_files) == 0) {
  stop("No non-vegetated mask files found in: ", nonveg_dir)
}
# Extract and sort by year
get_year <- function(p) {
  as.integer(sub("^.*_CCI_(\\d{4})_.*$", "\\1", basename(p)))
}
yrs <- vapply(mask_files, get_year, integer(1))
ord <- order(yrs)
mask_files <- mask_files[ord]
yrs <- yrs[ord]
# area raster (km²)
area005 <- rast(cfg$grids$grid_005$area_raster)
support_dom <- is.finite(area005) & area005 > 0
support_area_km2 <- as.numeric(global(ifel(support_dom, area005, NA), "sum", na.rm = TRUE)[1, 1])
if (!is.finite(support_area_km2) || support_area_km2 <= 0) {
  stop("Invalid support-domain area.")
}
# read masks and compute removed area
read_mask <- function(p) {
  m <- rast(p)
  compareGeom(area005, m, stopOnError = TRUE)
  m
}

masks <- lapply(mask_files, read_mask)
names(masks) <- as.character(yrs)
removed_km2 <- vapply(masks, function(m) {
  as.numeric(global(ifel(support_dom &
    (m == 1), area005, NA), "sum", na.rm = TRUE)[1, 1])
}, numeric(1))
removed_pct_support <- 100 * removed_km2 / support_area_km2
summary_tbl <- data.frame(
  year = yrs,
  file = basename(mask_files),
  removed_km2 = round(removed_km2, 0),
  removed_pct_valid_domain = round(removed_pct_support, 4)
)

# pairwise disagreement
if (length(masks) < 2) {
  pair_tbl <- data.frame(
    year_a = integer(0),
    year_b = integer(0),
    disagree_km2 = numeric(0),
    disagree_pct_valid_domain = numeric(0),
    area_agreement_pct = numeric(0)
  )
} else {
  # Generate all pairs (i, j) where i < j
  idx_pairs <- combn(seq_along(masks), 2)
  pair_list <- list()
  for (k in seq_len(ncol(idx_pairs))) {
    i <- idx_pairs[1, k]
    j <- idx_pairs[2, k]
    mi <- masks[[i]]
    mj <- masks[[j]]
    disagree <- support_dom & (mi != mj)
    disagree_km2 <- as.numeric(global(ifel(disagree, area005, NA), "sum", na.rm = TRUE)[1, 1])
    disagree_pct <- 100 * disagree_km2 / support_area_km2
    agree_rate <- 100 - disagree_pct
    pair_list[[k]] <- data.frame(
      year_a = yrs[i],
      year_b = yrs[j],
      disagree_km2 = round(disagree_km2, 0),
      disagree_pct_valid_domain = round(disagree_pct, 5),
      area_agreement_pct = round(agree_rate, 4)
    )
  }
  pair_tbl <- do.call(rbind, pair_list)
}

# nearly identical check
max_disagree_pct <- if (nrow(pair_tbl) == 0) {
  NA_real_
} else {
  max(pair_tbl$disagree_pct_valid_domain, na.rm = TRUE)
}
max_removed_diff_pct <- diff(range(summary_tbl$removed_pct_valid_domain, na.rm = TRUE))

# save
outdir <- here("analysis", "results", "tables", "masks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  round_numeric(summary_tbl, 5),
  file.path(
    outdir,
    "nonvegetated_snapshot_removed_area_sensitivity.csv"
  ),
  row.names = FALSE
)
write.csv(
  round_numeric(pair_tbl, 5),
  file.path(outdir, "nonvegetated_snapshot_pairwise_disagreement.csv"),
  row.names = FALSE
)
