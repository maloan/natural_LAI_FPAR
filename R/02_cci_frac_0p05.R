# ## =============================================================================
# # 02_cci_frac_0p05.R — Aggregate ESA-CCI/C3S land cover to 0.05° fractional cover
# ## =============================================================================
#
# suppressPackageStartupMessages({
#   library(terra)
#   library(dplyr)
#   library(stringr)
#   library(glue)
#   library(tibble)
#   library(here)
# })
#
# # --- config & refs -------------------------------------------------------------
# ROOT <- here()
#
# source(here("R", "utils.R"))
# source(here("R", "io.R"))
# source(here("R", "geom.R"))
# source(here("R", "viz.R"))
# source(here("R", "options.R"))
#
# cfg  <- cfg_read()
# opts <- opts_read()
#
# terraOptions(progress = 1, memfrac = 0.25)
#
# # --- Directories --------------------------------------------------------------
# cci_dir <- cfg$paths$cci_dir
# out_dir <- cfg$paths$cci_out_dir
# ql_dir  <- file.path(out_dir, "quicklooks")
#
# dir.create(out_dir, TRUE, showWarnings = FALSE)
# dir.create(ql_dir,  TRUE, showWarnings = FALSE)
#
# # --- Settings --------------------------------------------------------------
# tmpl <- rast(cfg$grids$grid_005$ref_raster)
#
# # toggles
# REMAKE_ALL    <- as.logical(Sys.getenv("REMAKE_ALL", "FALSE"))
# REMAKE_QL     <- as.logical(Sys.getenv("REMAKE_QL", "FALSE"))
# SKIP_EXISTING <- as.logical(Sys.getenv("SKIP_EXISTING", "TRUE"))
#
# # --- class sets ---------------------------------------------------------------
# ESACCI <- cfg$esa_cci$classes
# nodata_vals <- unique(c(ESACCI$nodata, 255))
#
# groups <- list(
#   cropland = ESACCI$cropland,
#   urban    = ESACCI$urban,
#   cls30    = ESACCI$cls30,
#   cls40    = ESACCI$cls40,
#   grass    = ESACCI$grassland
# )
# groups <- Filter(Negate(is.null), groups)
#
# # --- choose one file per year -------------------------------------------------
# files <- list.files(cci_dir, "\\.tif$", full.names = TRUE) |>
#   tibble(path = _) |>
#   mutate(year   = as.integer(str_extract(basename(path), "(19|20)\\d{2}")),
#          source = if_else(str_detect(basename(path), "^C3S"), "C3S", "ESACCI")) |>
#   filter(
#     !is.na(year),
#     dplyr::between(year, cfg$project$years$cci_start, cfg$project$years$cci_end)
#   ) |>
#   arrange(desc(source)) |>
#   group_by(year) |>
#   slice_head(n = 1) |>
#   ungroup()
#
# # --- build plan table ----------------------------------------------------------
# plan <- files |>
#   mutate(
#     out_tif   = file.path(out_dir, glue("ESACCI_frac_{year}_0p05.tif")),
#     ql_global = file.path(ql_dir, sprintf("quicklook_global_%d.png", year)),
#     ql_aoi    = file.path(ql_dir, sprintf("quicklook_aoi_%d.png", year)),
#     have_tif  = file.exists(out_tif),
#     have_qg   = file.exists(ql_global),
#     have_qa   = file.exists(ql_aoi)
#   )
#
# # --- main loop -----------------------------------------------------------------
# for (i in seq_len(nrow(plan))) {
#
#   f  <- plan$path[i]
#   yr <- plan$year[i]
#
#   ot <- plan$out_tif[i]
#   qg <- plan$ql_global[i]
#   qa <- plan$ql_aoi[i]
#
#   need_tif <- (REMAKE_ALL || !(SKIP_EXISTING && file.exists(ot)))
#   need_qg  <- (REMAKE_ALL || REMAKE_QL || !(SKIP_EXISTING && file.exists(qg)))
#   need_qa  <- (REMAKE_ALL || REMAKE_QL || !(SKIP_EXISTING && file.exists(qa)))
#
#   if (!need_tif && !need_qg && !need_qa) {
#     message("✓ Year ", yr, " already complete — skipping.")
#     next
#   }
#
#   if (need_tif) {
#
#     message("→ Building fractions for ", yr, "  ←  ", basename(f))
#
#     r <- rast(f)
#
#     seen <- unique(values(r))
#     seen <- seen[is.finite(seen)]
#
#     allowed <- unique(unlist(ESACCI[c(
#       "cropland","urban","cls30","cls40","grassland","bare","water","snow_ice"
#     )], use.names = FALSE))
#
#     extra <- setdiff(seen, allowed)
#     if (length(extra))
#       warning("CCI unexpected codes in ", basename(f), ": ",
#               paste(head(extra, 20), collapse = ", "))
#
#     if (is.na(crs(r)))
#       crs(r) <- "EPSG:4326"
#
#     r[r %in% nodata_vals] <- NA
#
#     frac <- lapply(groups, function(cls) {
#       m <- classify(r, cbind(cls, 1), others = 0)
#       resample(m, tmpl, method = "average")
#     }) |> rast()
#
#     names(frac) <- paste0("frac_", names(groups))
#
#     w30 <- cfg$esa_cci$weights$cls30 %||% 0.75
#     w40 <- cfg$esa_cci$weights$cls40 %||% 0.25
#
#     fc  <- if ("frac_cropland" %in% names(frac)) frac$frac_cropland else 0
#     fu  <- if ("frac_urban"    %in% names(frac)) frac$frac_urban    else 0
#     f30 <- if ("frac_cls30"    %in% names(frac)) frac$frac_cls30    else 0
#     f40 <- if ("frac_cls40"    %in% names(frac)) frac$frac_cls40    else 0
#
#     frac_fused <- clamp(fc + fu + w30*f30 + w40*f40, 0, 1)
#     names(frac_fused) <- "frac_fused"
#
#     writeRaster(
#       frac_fused,
#       ot,
#       overwrite = TRUE,
#       gdal   = gdal_wopt("FLT4S")$gdal,
#       NAflag = -9999
#     )
#
#   } else {
#     frac_fused <- rast(ot)
#   }
#
#   if (need_qg || need_qa) {
#     quicklook_all_aois(
#       frac    = frac_fused,
#       year    = yr,
#       cfg     = cfg,
#       ql_root = ql_dir,
#       down    = 4L,
#       include_global   = TRUE,
#       drop_global_key  = FALSE
#     )
#   }
# }
#
# gc()
# message("Done.")
