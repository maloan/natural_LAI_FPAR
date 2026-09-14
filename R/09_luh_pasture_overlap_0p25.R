# =============================================================================
# 09_luh_pasture_overlap_0p25.R — Build LUH pasture-overlap mask (0.25° + 0.05°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

source(here("R", "helpers", "plotting.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

#  params
grass_source <- toupper(Sys.getenv("GRASS_SOURCE", "GLC"))
g_min <- cfg$luh2$pasture_mask$grass_min
p_min <- cfg$luh2$pasture_mask$pasture_min
beta <- cfg$luh2$pasture_mask$pasture_grass_ratio_min

# LUH2 v2h historical averaging window used in the manuscript.
year_0 <- as.integer(cfg$luh2$pasture_mask$start_year)
year_1 <- as.integer(cfg$luh2$pasture_mask$end_year)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)
area005 <- rast(cfg$grids$grid_005$area_raster)

out_dir <- file.path(cfg$paths$masks_root_dir, "mask_luh_overlap")
ql_dir <- file.path(out_dir, "quicklooks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ql_dir, recursive = TRUE, showWarnings = FALSE)

year_span <- function(y) {
  # Helper to summarize a vector of years as a range string
  y <- as.integer(y)
  y <- y[!is.na(y)]
  if (!length(y)) {
    return("unknown")
  }
  paste(range(y), collapse = "-")
}
tag <- sprintf(
  "%s_Gmin%s_Pmin%s_beta%s_%d-%d",
  grass_source,
  tok(g_min),
  tok(p_min),
  tok(beta),
  year_0,
  year_1
)

out025 <- file.path(out_dir, sprintf("mask_luh_overlap_%s_0p25.tif", tag))
out005 <- file.path(out_dir, sprintf("mask_luh_overlap_%s_0p05_rep.tif", tag))

#  1) grass fraction @0.05°
grass_005 <- switch(grass_source,
  CCI = {
    frac_dir <- cfg$paths$cci_out_dir
    yrs <- year_0:year_1
    ff <- file.path(frac_dir, sprintf("ESACCI_frac_%d_0p05.tif", yrs))
    missing_files <- ff[!file.exists(ff)]
    if (length(missing_files)) {
      stop("Missing expected CCI fraction files:\n", paste(missing_files, collapse = "\n"))
    }
    stk <- rast(lapply(ff, function(x) {
      rast(x)[["frac_grass"]]
    }))
    mean(stk, na.rm = TRUE)
  },
  GLC = {
    p <- file.path(cfg$paths$glc_out_dir, "glc_grass_yearstack_0p05.tif")
    if (!file.exists(p)) {
      stop_msg("Missing GLC grass-fraction stack: ", p)
    }
    s <- rast(p)
    yrs <- as.integer(substr(names(s), 2, 5))
    keep <- filter_by_year_range(yrs, year_0, year_1)
    if (!length(keep)) {
      stop_msg(
        "No GLC layers in requested window ",
        year_0,
        "-",
        year_1,
        ". Available years: ",
        year_span(yrs)
      )
    }
    mean(s[[keep]], na.rm = TRUE)
  },
  stop_msg("Unknown grass_source: ", grass_source)
)

grass_005 <- align_to_template(grass_005, ref005, method = "bilinear")
grass_005 <- clamp(grass_005, 0, 1)
names(grass_005) <- "grass_005"

#  2) area-weighted aggregate grass 0.05° → 0.25°
grass_025 <- agg005_to_025_aw(grass_005, area005, ref025)
grass_025 <- clamp(grass_025, 0, 1)
names(grass_025) <- "grass_025"

#  3) LUH pasture @0.25° (mean over window)
luh_nc <- cfg$luh2$states_nc
v_pas <- cfg$luh2$variables$pasture
pas <- rast(luh_nc, subds = v_pas)
times <- suppressWarnings(as.integer(time(pas)))
if (!length(times) || all(is.na(times))) {
  stop_msg("LUH pasture time axis missing or not parseable as integer years")
}
keep <- filter_by_year_range(times, year_0, year_1)
if (!length(keep)) {
  stop_msg(
    "No LUH pasture layers in requested window ",
    year_0,
    "-",
    year_1,
    ". Available years: ",
    year_span(times)
  )
}
pasture_025 <- clamp(mean(pas[[keep]], na.rm = TRUE), 0, 1)
pasture_025 <- align_to_template(pasture_025, ref025, method = "bilinear")
names(pasture_025) <- "pasture_025"

#  4) decision at 0.25°
ratio <- clamp(pasture_025 / (grass_025 + 1e-9), 0, 1)
drop_025 <- (grass_025 >= g_min) &
  (pasture_025 >= p_min) & (ratio >= beta)
mask_025 <- ifel(drop_025, 1L, 0L)

#  5) write 0.25° + 0.05° replica
wopt <- wopt_byte(speed_over_size = TRUE, na = 255L)

writeRaster(mask_025, out025, overwrite = TRUE, wopt = wopt)

mask_005 <- disagg(mask_025, fact = 5, method = "near")
mask_005 <- align_to_template(mask_005, ref005, method = "near")
writeRaster(mask_005, out005, overwrite = TRUE, wopt = wopt)

#  6) quicklooks
write_3panel <- function(g, p, m, out_png, main) {
  # Helper to write a 3-panel quicklook of grass, pasture, and mask
  png(out_png, 1600, 750, res = 120)
  op <- par(
    mfrow = c(1, 3),
    mar = c(3, 3, 3, 4),
    oma = c(2, 0, 2, 0)
  )
  on.exit(
    {
      par(op)
      dev.off()
    },
    add = TRUE
  )

  plot(g,
    main = "Grass (0.25°)",
    col = pal_green(64),
    zlim = c(0, 1)
  )
  plot(p,
    main = "Pasture (0.25°)",
    col = pal_green(64),
    zlim = c(0, 1)
  )
  plot(
    m,
    main = "Mask (1=drop)",
    col = c("#f0f0f0", "#d73027"),
    breaks = c(-0.5, 0.5, 1.5),
    legend = FALSE,
    axes = TRUE,
    box = TRUE
  )
  legend(
    "bottomleft",
    fill = c("#f0f0f0", "#d73027"),
    legend = c("0 keep", "1 drop"),
    bty = "n"
  )
  mtext(main,
    outer = TRUE,
    line = 0,
    cex = 1.2
  )
}

ql_global <- file.path(ql_dir, sprintf("quicklook_global_%s.png", tag))
if (!file.exists(ql_global)) {
  write_3panel(
    grass_025,
    pasture_025,
    mask_025,
    ql_global,
    sprintf(
      "LUH overlap (α=%s, G≥%s, P≥%s), %s, %d–%d",
      tok(beta),
      tok(g_min),
      tok(p_min),
      grass_source,
      year_0,
      year_1
    )
  )
}

gc()
