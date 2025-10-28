## =============================================================================
# 02_cci_frac_0p05.R — Aggregate ESA-CCI/C3S land cover to 0.05° fractions

# Load packages
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(stringr)
  library(glue)
  library(tibble)
})

# --- config & refs -------------------------------------------------------------
# Root dir
ROOT <- normalizePath(path.expand(
  Sys.getenv("SNU_LAI_FPAR_ROOT", "~/GitHub/natural_LAI_FPAR")
),
winslash = "/",
mustWork = FALSE)
# Other source files
source(file.path(ROOT, "R", "00_utils.R"))
source(file.path(ROOT, "R", "io.R"))
source(file.path(ROOT, "R", "geom.R"))
source(file.path(ROOT, "R", "viz.R"))
cfg <- cfg_read()
opts <- opts_read()
terraOptions(progress = 1, memfrac = 0.25)

# --- Directories --------------------------------------------------------------
cci_dir <- cfg$paths$cci_dir
out_dir <- path.expand(cfg$paths$cci_out_dir)
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, TRUE, showWarnings = FALSE)
dir.create(ql_dir, TRUE, showWarnings = FALSE)

# --- Settings --------------------------------------------------------------
tmpl    <- rast(cfg$grids$grid_005$ref_raster) # raster
# --- rebuild policy (env toggles) ---------------------------------------------
REMAKE_ALL <- as.logical(Sys.getenv("REMAKE_ALL", "FALSE"))
REMAKE_QL  <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))

# --- class sets ---------------------------------------------------------------
# load ESACCI classes
ESACCI <- cfg$esa_cci$classes
# define no data value
nodata_vals <- unique(c(ESACCI$nodata, 255))
# create groups with the different classes =
groups <- list(
  cropland  = ESACCI$cropland,
  urban     = ESACCI$urban,
  cls30     = ESACCI$cls30,
  cls40     = ESACCI$cls40,
  grass     = ESACCI$grassland
)
# remove entries that are NULL
groups <- Filter(Negate(is.null), groups)

# --- choose one file per year (prefer C3S), restrict to cfg window ------------
# list all tif files in the cci directory
files <- list.files(cci_dir, "\\.tif$", full.names = TRUE) |>
  tibble(path = _) |>
  # Parse year from filename and Tag source as "C3S" if basename starts with
  # C3S, else "ESACCI".
  mutate(year   = as.integer(str_extract(basename(path), "(19|20)\\d{2}")),
         source = if_else(str_detect(basename(path), "^C3S"), "C3S", "ESACCI")) |>
  # only keep rows with valid year
  filter(
    !is.na(year),
    dplyr::between(year, cfg$project$years$cci_start, cfg$project$years$cci_end)
  ) |>
  # Order by source descending (so C3S before ESACCI)
  arrange(desc(source)) |>
  # Group by year; take the first row per group → one file/year.
  group_by(year) |>
  slice_head(n = 1) |>
  ungroup()

# --- 300 m → 0.05° fractions ---------------------------------------------------
# build a plan table from files, so loop knows what to (re)build per year
plan <- files |>
  mutate(
    out_tif   = file.path(out_dir, glue("ESACCI_frac_{year}_0p05.tif")),
    ql_global = file.path(ql_dir, sprintf("quicklook_global_%d.png", year)),
    ql_aoi    = file.path(ql_dir, sprintf("quicklook_aoi_%d.png", year)),
    have_tif  = file.exists(out_tif),
    have_qg   = file.exists(ql_global),
    have_qa   = file.exists(ql_aoi)
  )

for (i in seq_len(nrow(plan))) {
  # Iterate over each planned year, pull the row's paths and flags and define what to rebuild
  f   <- plan$path[i] # input path
  yr  <- plan$year[i] # year
  ot  <- plan$out_tif[i] # output path
  qg  <- plan$ql_global[i] # quicklook path
  qa  <- plan$ql_aoi[i] # quicklook paths

  # build the GeoTIFF if REMAKE_ALL is TRUE or we’re not in “skip existing” mode or the file doesn’t exist
  need_tif <- (REMAKE_ALL || !(SKIP_EXISTING && file.exists(ot)))
  # same logic for the global and AOI quicklooks, but they also rebuild when REMAKE_QL is TRUE
  need_qg  <- (REMAKE_ALL ||
                 REMAKE_QL || !(SKIP_EXISTING && file.exists(qg)))
  need_qa  <- (REMAKE_ALL ||
                 REMAKE_QL || !(SKIP_EXISTING && file.exists(qa)))

  if (!need_tif && !need_qg && !need_qa) {
    # Do not need anything -> do nothing
    message("✓ Year ", yr, " already complete — skipping.")
    next
  }

  if (need_tif) {
    # Need tif -> build fractions
    message("→ Building fractions for ", yr, "  ←  ", basename(f))
    # Load input
    r <- rast(f)
    # warn if unexpected classes present (ignores NA)
    seen <- unique(values(r))
    seen <- seen[is.finite(seen)]
    allowed <- unique(unlist(ESACCI[c("cropland",
                                      "urban",
                                      "cls30",
                                      "cls40",
                                      "grassland",
                                      "bare",
                                      "water",
                                      "snow_ice")], use.names = FALSE))
    extra <- setdiff(seen, allowed)
    if (length(extra)) {
      warning("CCI unexpected codes in ",
              basename(f),
              ": ",
              paste(head(extra, 20), collapse = ", "))
    }

    if (is.na(crs(r))) {
      # set CRS to EPSG:4326 if missing
      crs(r) <- "EPSG:4326"
    }
    # turn nodata codes into NA
    r[r %in% nodata_vals] <- NA

    # --- build fused fraction then strict mask ---
    frac <- lapply(groups, function(cls) {
      # For each class group in groups
      m <- classify(r, cbind(cls, 1), others = 0) # makes a binary mask (1 = class in group, 0 = other)
      resample(m, tmpl, method = "average") # aggregates 300 m → 0.05°; the mean of 0/1 equals fractional cover
    }) |> rast()
    # Stack results
    names(frac) <- paste0("frac_", names(groups))
    # Pull weights (w30, w40, defaults 0.75/0.25).
    w30 <- cfg$esa_cci$weights$cls30 %||% 0.75
    w40 <- cfg$esa_cci$weights$cls40 %||% 0.25
    # Extract available layers (fc cropland, fu urban, f30 cls30, f40 cls40; else 0).
    fc  <- if ("frac_cropland" %in% names(frac)) {
      frac$frac_cropland
    } else {
      0
    }
    fu  <- if ("frac_urban"    %in% names(frac)) {
      frac$frac_urban
    }   else {
      0
    }
    f30 <- if ("frac_cls30"    %in% names(frac)) {
      frac$frac_cls30
    }    else{
      0
    }
    f40 <- if ("frac_cls40"    %in% names(frac)) {
      frac$frac_cls40
    }   else {
      0
    }

    # Define frac mask
    frac_fused <- clamp(fc + fu + w30 * f30 + w40 * f40, 0, 1)
    names(frac_fused) <- "frac_fused"
    writeRaster(
      frac_fused,
      ot,
      overwrite = TRUE,
      gdal   = gdal_wopt("FLT4S")$gdal,
      # Float32 GeoTIFF
      NAflag = -9999
    )

  } else {
    # else-branch (reopen existing mask)
    frac_fused <- rast(ot)
  }

  # (Re)draw quicklooks if needed (Global + all AOIs from cfg$aois)
  if (need_qg || need_qa) {
    quicklook_all_aois(
      frac    = frac_fused,
      year    = yr,
      cfg     = cfg,
      ql_root = ql_dir,
      down    = 4L,
      include_global   = TRUE,
      drop_global_key  = FALSE
    )
  }

}
message("Done.")
