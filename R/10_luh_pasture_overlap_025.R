## =============================================================================
# 10_luh_pasture_overlap_025.R — Build LUH pasture-overlap mask (0.25° + 0.05°)
## =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

source(here("R", "helpers", "utils.R"))
source(here("R", "helpers", "io.R"))
source(here("R", "helpers", "geom.R"))
source(here("R", "helpers", "viz.R"))
source(here("R", "helpers", "options.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

# --- params -------------------------------------------------------------------
GRASS_SOURCE <- toupper(Sys.getenv("GRASS_SOURCE", "CCI")) # CCI | GLC_TEMP
REMAKE_QL <- as.logical(Sys.getenv("REMAKE_QL", "TRUE"))

G_MIN <- as.numeric(Sys.getenv("G_MIN", "0.1"))
P_MIN <- as.numeric(Sys.getenv("P_MIN", "0.1"))
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.5"))

Y0 <- env_get_int("LUH_AVG_START", cfg$project$years$cci_start)
Y1 <- env_get_int("LUH_AVG_END", cfg$project$years$cci_end)

stopifnot(
  is.finite(G_MIN), is.finite(P_MIN),
  is.finite(ALPHA), is.finite(Y0), is.finite(Y1)
)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

tok2 <- function(x) gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
tag <- sprintf(
  "%s_Gmin%s_Pmin%s_alpha%s_%d-%d",
  GRASS_SOURCE, tok2(G_MIN), tok2(P_MIN), tok2(ALPHA), Y0, Y1
)

out025 <- file.path(out_dir, sprintf("mask_luh_overlap_%s_0p25.tif", tag))
out005 <- file.path(out_dir, sprintf("mask_luh_overlap_%s_0p05_rep.tif", tag))

# --- 1) grass fraction @0.05° -------------------------------------------------
grass_005 <- switch(GRASS_SOURCE,
  CCI = {
    frac_dir <- cfg$paths$cci_out_dir
    ff <- list.files(frac_dir,
      pattern = "ESACCI_frac_\\d{4}_0p05\\.tif$",
      full.names = TRUE
    )
    if (!length(ff)) stop("No CCI fraction files in: ", frac_dir)
    yrs <- as.integer(substr(
      basename(ff), regexpr("\\d{4}", basename(ff)),
      regexpr("\\d{4}", basename(ff)) + 3
    ))
    keep <- which(yrs >= Y0 & yrs <= Y1)
    if (!length(keep)) stop("No CCI fraction years in window ", Y0, "-", Y1)
    stk <- rast(lapply(ff[keep], function(x) rast(x)[["frac_grass"]]))
    mean(stk, na.rm = TRUE)
  },
  GLC_TEMP = {
    p <- file.path(cfg$paths$glc_out_dir, "glc_cat_yearstack_0p05.tif")
    if (!file.exists(p)) stop("GLC yearstack not found: ", p)
    s <- rast(p)
    yrs <- suppressWarnings(as.integer(substr(names(s), 2, 5)))
    keep <- which(yrs >= Y0 & yrs <= Y1)
    if (!length(keep)) stop("No GLC years in window ", Y0, "-", Y1)
    grass_vals <- as.integer(
      unlist(cfg$glc$classes$grassland, use.names = FALSE)
    )
    is_grass <- classify(s[[keep]], cbind(grass_vals, 1), others = 0)
    app(is_grass, mean, na.rm = TRUE)
  },
  stop("Unknown GRASS_SOURCE: ", GRASS_SOURCE)
)

if (!compareGeom(grass_005, ref005, stopOnError = FALSE)) {
  grass_005 <- resample(grass_005, ref005, method = "bilinear")
}
grass_005 <- clamp(grass_005, 0, 1)
names(grass_005) <- "grass_005"

# --- 2) area-weighted aggregate grass 0.05° → 0.25° ---------------------------
grass_025 <- agg005_to_025_aw(grass_005, area005, ref025)
grass_025 <- clamp(grass_025, 0, 1)
names(grass_025) <- "grass_025"

# --- 3) LUH pasture @0.25° (mean over window) ---------------------------------
luh_nc <- cfg$luh2$states_nc
if (!file.exists(luh_nc)) stop("LUH file not found: ", luh_nc)

v_pas <- cfg$luh2$variables$pasture
pas <- rast(luh_nc, subds = v_pas)

ty <- suppressWarnings(as.integer(time(pas)))
keep <- which(ty >= Y0 & ty <= Y1)
if (!length(keep)) stop("No LUH timesteps in window ", Y0, "-", Y1)

pasture_025 <- clamp(mean(pas[[keep]], na.rm = TRUE), 0, 1)
if (!compareGeom(pasture_025, ref025, stopOnError = FALSE)) {
  pasture_025 <- resample(pasture_025, ref025, method = "bilinear")
}
names(pasture_025) <- "pasture_025"

# --- 4) decision at 0.25° ------------------------------------------------------
ratio <- clamp(pasture_025 / (grass_025 + 1e-9), 0, 1)
drop_025 <- (grass_025 >= G_MIN) & (pasture_025 >= P_MIN) & (ratio >= ALPHA)
mask_025 <- ifel(drop_025, 1L, 0L)

# --- 5) write 0.25° + 0.05° replica -------------------------------------------
wopt <- wopt_byte(speed = TRUE, na = 255L)

writeRaster(mask_025, out025,
  overwrite = TRUE, datatype = wopt$datatype,
  gdal = wopt$gdal, NAflag = wopt$NAflag
)

mask_005 <- disagg(mask_025, fact = 5, method = "near")
if (!compareGeom(mask_005, ref005, stopOnError = FALSE)) {
  mask_005 <- resample(mask_005, ref005, method = "near")
}
writeRaster(mask_005, out005,
  overwrite = TRUE,
  datatype = wopt$datatype, gdal = wopt$gdal, NAflag = wopt$NAflag
)

# --- 6) quicklooks (one function; global + AOIs) ------------------------------
write_3panel <- function(g, p, m, out_png, main) {
  png(out_png, 1600, 750, res = 120)
  op <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 4), oma = c(2, 0, 2, 0))
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  plot(g, main = "Grass (0.25°)", col = pal_green(64), zlim = c(0, 1))
  plot(p, main = "Pasture (0.25°)", col = pal_green(64), zlim = c(0, 1))
  plot(m,
    main = "Mask (1=drop)", col = c("#f0f0f0", "#d73027"),
    breaks = c(-0.5, 0.5, 1.5), legend = FALSE, axes = TRUE, box = TRUE
  )
  legend("bottomleft",
    fill = c("#f0f0f0", "#d73027"),
    legend = c("0 keep", "1 drop"), bty = "n"
  )
  mtext(main, outer = TRUE, line = 0, cex = 1.2)
}

ql_global <- file.path(ql_dir, sprintf("quicklook_global_%s.png", tag))
if (REMAKE_QL || !file.exists(ql_global)) {
  write_3panel(
    grass_025, pasture_025, mask_025, ql_global,
    sprintf(
      "LUH overlap (α=%s, G≥%s, P≥%s), %s, %d–%d",
      tok2(ALPHA), tok2(G_MIN), tok2(P_MIN), GRASS_SOURCE, Y0, Y1
    )
  )
}

aoi_root <- file.path(ql_dir, "aois")
dir.create(aoi_root, recursive = TRUE, showWarnings = FALSE)

for (nm in names(cfg$aois)) {
  a <- cfg$aois[[nm]]
  ext_aoi <- ext(a$lon_min, a$lon_max, a$lat_min, a$lat_max)
  g_aoi <- try(crop(grass_025, ext_aoi, snap = "near"), silent = TRUE)
  p_aoi <- try(crop(pasture_025, ext_aoi, snap = "near"), silent = TRUE)
  m_aoi <- try(crop(mask_025, ext_aoi, snap = "near"), silent = TRUE)
  if (inherits(g_aoi, "try-error") || inherits(p_aoi, "try-error") ||
    inherits(m_aoi, "try-error")) {
    next
  }

  ql_aoi <- file.path(aoi_root, sprintf("quicklook_%s_%s.png", nm, tag))
  if (REMAKE_QL || !file.exists(ql_aoi)) {
    write_3panel(g_aoi, p_aoi, m_aoi, ql_aoi, sprintf(
      "AOI: %s — %d–%d",
      nm, Y0, Y1
    ))
  }
}

gc()
cat(sprintf(
  "Wrote:\n  - %s\n  - %s\n
  Rule: drop if grass≥%.3f & pasture≥%.3f & pasture/grass≥%.3f\n",
  out025, out005, G_MIN, P_MIN, ALPHA
))
