# =============================================================================
# 08_luh_use_masks.R — Build LUH pasture and rangeland share maps (0.25° +
# 0.05°)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})


source(here("R", "helpers", "netcdf.R"))
source(here("R", "helpers", "io.R"))

cfg <- cfg_read()

terraOptions(progress = 1, memfrac = 0.25)

year_0 <- as.integer(Sys.getenv("LUH_AVG_START", cfg$project$years$cci_start))
year_1 <- as.integer(Sys.getenv("LUH_AVG_END", cfg$project$years$cci_end))
tag <- sprintf("%d-%d", year_0, year_1)

ref005 <- rast(cfg$grids$grid_005$ref_raster)
ref025 <- rast(cfg$grids$grid_025$ref_raster)

out_dir <- cfg$paths$masks_luh_dir
wopt <- wopt_f32(FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

luh_nc <- cfg$luh2$states_nc
v_pas <- cfg$luh2$variables$pasture
v_rng <- cfg$luh2$variables$rangeland %||% cfg$luh2$variables$range

pas <- rast(luh_nc, subds = v_pas)
rng <- rast(luh_nc, subds = v_rng)
times <- suppressWarnings(as.integer(time(pas)))
idx <- filter_by_year_range(times, year_0, year_1)
if (!length(idx)) {
  stop("No LUH steps in ", year_0, "-", year_1)
}

m025_pas <- clamp(mean(pas[[idx]], na.rm = TRUE), 0, 1)
m025_rng <- clamp(mean(rng[[idx]], na.rm = TRUE), 0, 1)
names(m025_pas) <- "pasture_share"
names(m025_rng) <- "rangeland_share"

# align once to canonical 0.25° grid
m025_pas <- align_to_template(m025_pas, ref025, method = "bilinear")
m025_rng <- align_to_template(m025_rng, ref025, method = "bilinear")

f_pas_025 <- file.path(out_dir, sprintf("m_LUH_pasture_%s_0p25.tif", tag))
f_rng_025 <- file.path(out_dir, sprintf("m_LUH_rangeland_%s_0p25.tif", tag))

writeRaster(m025_pas, f_pas_025, overwrite = TRUE, wopt = wopt)
writeRaster(m025_rng, f_rng_025, overwrite = TRUE, wopt = wopt)

# 0.05° replicas (nearest-neighbour replication)
m005_pas <- disagg(m025_pas, fact = 5, method = "near")
m005_rng <- disagg(m025_rng, fact = 5, method = "near")

m005_pas <- align_to_template(m005_pas, ref005, method = "near")
m005_rng <- align_to_template(m005_rng, ref005, method = "near")

f_pas_005 <- file.path(out_dir, sprintf("m_LUH_pasture_%s_0p05_rep.tif", tag))
f_rng_005 <- file.path(out_dir, sprintf("m_LUH_rangeland_%s_0p05_rep.tif", tag))

writeRaster(m005_pas, f_pas_005, overwrite = TRUE, wopt = wopt)
writeRaster(m005_rng, f_rng_005, overwrite = TRUE, wopt = wopt)

gc()
