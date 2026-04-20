# =============================================================================
# 07_nonveg_snapshot_sensitivity.R
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "paths.R"))
source(here("R", "helpers", "files.R"))
source(here("R", "helpers", "netcdf.R"))
cfg <- cfg_read()

# --- locate files -------------------------------------------------------------
run_tag <- cfg$project$run_tag
nonveg_dir <- here("output", run_tag, "masks", "mask_nonvegetated")

# Read all matching mask files
mask_files <- list.files(
  nonveg_dir,
  pattern = "^mask_nonvegetated_CCI_\\d{4}_tauW0p05_tauI0p05_0p05\\.tif$",
  full.names = TRUE
)
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

# --- area raster (km²) --------------------------------------------------------
area005 <- rast(cfg$grids$grid_005$area_raster)
global_area_km2 <- as.numeric(global(area005, "sum", na.rm = TRUE)[1, 1])

# --- read masks and compute removed area --------------------------------------
read_mask <- function(p) {
  m <- rast(p)
  compareGeom(area005, m, stopOnError = TRUE)
  m
}

masks <- lapply(mask_files, read_mask)
names(masks) <- as.character(yrs)

removed_km2 <- vapply(masks, function(m) {
  as.numeric(global(area005 * (m == 1), "sum", na.rm = TRUE)[1, 1])
}, numeric(1))

removed_pct_global <- 100 * removed_km2 / global_area_km2
summary_tbl <- data.frame(
  year = yrs,
  file = basename(mask_files),
  removed_km2 = round(removed_km2, 0),
  removed_pct_global = round(removed_pct_global, 4)
)

cat("\n=== Non-vegetated snapshot masks: removed area ===\n")
print(summary_tbl, row.names = FALSE)

# --- pairwise disagreement ----------------------------------------------------
if (length(masks) < 2) {
  pair_tbl <- data.frame(
    year_a = integer(0),
    year_b = integer(0),
    disagree_km2 = numeric(0),
    disagree_pct_global = numeric(0),
    cell_agreement_pct = numeric(0)
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

    disagree <- (mi != mj)
    disagree_km2 <- as.numeric(
      global(area005 * disagree, "sum", na.rm = TRUE)[1, 1]
    )
    agree_rate <- 100 * (1 - global(disagree, "mean", na.rm = TRUE)[1, 1])

    pair_list[[k]] <- data.frame(
      year_a = yrs[i],
      year_b = yrs[j],
      disagree_km2 = round(disagree_km2, 0),
      disagree_pct_global = round(100 * disagree_km2 / global_area_km2, 5),
      cell_agreement_pct = round(agree_rate, 4)
    )
  }
  pair_tbl <- do.call(rbind, pair_list)
}

cat("\n=== Pairwise disagreement (XOR) ===\n")
print(pair_tbl, row.names = FALSE)

# ---  nearly identical check ------------------------------------------
max_disagree_pct <- if (nrow(pair_tbl) == 0) {
  NA_real_
} else {
  max(pair_tbl$disagree_pct_global, na.rm = TRUE)
}

max_removed_diff_pct <- if (nrow(summary_tbl) < 2) {
  NA_real_
} else {
  max(abs(diff(summary_tbl$removed_pct_global)), na.rm = TRUE)
}

cat("\n=== Quick diagnostics ===\n")
cat(sprintf(
  "Max pairwise disagreement (%% of global surface): %.6f\n", max_disagree_pct
))
cat(sprintf(
  "Max adjacent-year difference in removed area (percentage points of global): %.6f\n",
  max_removed_diff_pct
))
