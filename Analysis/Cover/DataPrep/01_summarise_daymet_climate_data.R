# //////////////////////////////////////////////////////////////////////////
# 01b_summarise_daymet_climate_data.R
#
# Derive historical climate "normals" (the _CLIM variable set) from Daymet v4
# monthly data, on the CONUS Daymet 1 km grid, using exactly the same metric
# definitions and across-year reductions as
# 01_summarise_projected_climate_data.R (MACA projections).
#
# The output grid is the CONUS mask/snap raster returned by read_mask()
# (created in Analysis/Cover/00_create_snap_raster.R). All outputs are cropped
# and masked to that grid.
#
# Inputs: monthly Daymet v4 (R1) North America GeoTIFFs, 12 bands per file
#   (Jan..Dec), under paths$large/Data_raw/daymet/rawMonthlyData:
#     daymet_v4_prcp_monttl_na_<year>.tif  (mm, monthly total)
#     daymet_v4_tmax_monavg_na_<year>.tif  (deg C, monthly mean of daily max)
#     daymet_v4_tmin_monavg_na_<year>.tif  (deg C, monthly mean of daily min)
#
#
# started August 2026
# //////////////////////////////////////////////////////////////////////////

library(terra)
library(stringr)

source("Functions/init.R")
source_functions() # Functions/data/climate.R contains the key functions used here


# terra scratch space -----------------------------------------------------
# The per-year metrics and across-year reductions spill large intermediate
# tiles to terra's temp dir. Point that at paths$large (2 TB) instead of the
# system C: drive, which otherwise fills and triggers
# "_tiffWriteProc: No space left on device".
terra_tmp <- file.path(paths$large, "terra_tmp")
dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
# todisk = TRUE forces block-wise processing to disk instead of allocating
# whole rasters in RAM; this is what prevents the std::bad_alloc failures on
# the full CONUS 1 km grid. memfrac caps the share of RAM terra will use.
terra::terraOptions(tempdir = terra_tmp, todisk = TRUE, memfrac = 0.5)


# Parameters --------------------------------------------------------------

year_start <- 1991        # first year of the climate-normal window
year_end   <- 2020        # last year
years      <- year_start:year_end

# Daymet variable file tag -> our short name.
daymet_vars <- c(prcp_monttl = "prcp", tmax_monavg = "tmax", 
                 tmin_monavg = "tmin", vp_monavg = "vp")

daymet_dir <- file.path(paths$large, "Data_raw/daymet/rawMonthlyData")

intermediate_dir <- file.path(paths$large, "Data_processed/WallToWallClimateData",
                              "DaymetClimate_intermediate")

out_dir <- file.path(paths$large, "Data_processed/WallToWallClimateData")

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Target grid/mask (CONUS-clipped Daymet 1 km grid).
mask_r <- read_mask()


# Helper: locate a Daymet monthly file ------------------------------------

# All .tif files under the raw monthly directory (the order-hash subdirectory
# name is not hard-coded).
daymet_tifs <- list.files(daymet_dir, pattern = "\\.tif$",
                          recursive = TRUE, full.names = TRUE)

#' Path to the Daymet monthly GeoTIFF for one variable and year
#'
#' @param var_tag Daymet file tag (e.g. "prcp_monttl").
#' @param yr Year (integer).
#' @return Single file path.
daymet_path <- function(var_tag, yr) {
  pattern <- paste0("daymet_v4_", var_tag, "_na_", yr, "\\.tif$")
  hits <- str_subset(daymet_tifs, pattern)
  if (length(hits) != 1) {
    stop("Expected 1 file for ", var_tag, " ", yr, ", found ", length(hits))
  }
  hits
}

# Fail early if any year is missing.
invisible(lapply(years, function(yr) {
  lapply(names(daymet_vars), daymet_path, yr = yr)
}))


# Helper: read one variable-year and put it on the mask grid ---------------

#' Read a 12-band Daymet monthly raster and align it to the CONUS mask
#'
#' Crops (and masks) if the source already shares CRS, resolution and origin
#' with the mask, otherwise projects. Projection should not be needed for
#' native Daymet files; a message is emitted if it is.
#'
#' @param var_tag Daymet file tag (e.g. "tmin_monavg").
#' @param yr Year (integer).
#' @return 12-layer SpatRaster on the mask grid, masked to CONUS.
read_daymet_year <- function(var_tag, yr) {
  r <- terra::rast(daymet_path(var_tag, yr))
  stopifnot(terra::nlyr(r) == 12)

  same_grid <- terra::compareGeom(r, mask_r, crs = TRUE, res = TRUE,
                                  ext = FALSE, rowcol = FALSE,
                                  stopOnError = FALSE) &&
    all(abs(terra::origin(r) - terra::origin(mask_r)) < 1e-6)

  if (same_grid) {
    r <- terra::crop(r, mask_r)
  } else {
    message("  ", var_tag, " ", yr, ": grid differs from mask, projecting")
    r <- terra::project(r, mask_r, method = "bilinear")
  }
  r <- terra::mask(r, mask_r)

  # Robustness: every variable-year must land on the mask grid cell-for-cell,
  # otherwise the per-year metrics and across-year stacking would silently
  # misalign. crs/res/ext/rowcol must all match; error out if they don't.
  if (!terra::compareGeom(r, mask_r, crs = TRUE, res = TRUE,
                          ext = TRUE, rowcol = TRUE, stopOnError = FALSE)) {
    stop("Grid mismatch after alignment for ", var_tag, " ", yr,
         ": result does not match the mask grid.")
  }
  r
}


# 1. Per-year annual metrics ----------------------------------------------
# For each year, read the 12 monthly layers of each variable, compute annual
# metrics, and write to an intermediate file. Re-running skips finished years.

for (yr in years) {
  out_file <- file.path(intermediate_dir, paste0("annualMetrics_", yr, ".tif"))
  if (file.exists(out_file)) {
    message("  year ", yr, " already done, skipping")
    next
  }
  message("  computing annual metrics for ", yr)

  metrics <- calc_annual_metrics(
    tmin12 = read_daymet_year("tmin_monavg", yr),
    tmax12 = read_daymet_year("tmax_monavg", yr),
    prcp12 = read_daymet_year("prcp_monttl", yr),
    vp12 = read_daymet_year("vp_monavg", yr)
  )
  
  terra::writeRaster(metrics, out_file, overwrite = TRUE)
  # Free per-year rasters so peak RAM does not grow across iterations.
  rm(metrics)
  gc()
}


# 2. Reduce across years to the _CLIM normals -----------------------------
# Most metrics: mean across years. A few use percentiles, matching training.


message("Loading intermediate rasters and reducing across years ...")
year_files  <- file.path(intermediate_dir, paste0("annualMetrics_", years, ".tif"))
year_stacks <- lapply(year_files, terra::rast)

# For a given per-year metric, build a multi-year stack (one layer per year).
stack_metric <- function(metric_name) {
  layers <- lapply(year_stacks, function(s) s[[metric_name]])
  do.call(c, layers)
}

clim_layers <- lapply(climate_reductions, function(spec) {
  metric_name <- spec[[1]]
  out_name    <- spec[[2]]
  fun         <- spec[[3]]
  r <- fun(stack_metric(metric_name))
  names(r) <- out_name
  r
})

clim <- do.call(c, clim_layers)

# Already on the mask grid; mask again so the percentile/mean reductions
# cannot introduce values outside CONUS.
clim <- terra::mask(clim, mask_r)


# 3. Write final raster ---------------------------------------------------

out_file <- file.path(out_dir,
                      paste0("DaymetClimateData_", year_start, "-", year_end,
                             "_CLIM.tif"))
terra::writeRaster(clim, out_file, overwrite = TRUE)
message("Wrote ", out_file)


# 4. Clean up terra scratch -----------------------------------------------
# The final output is now a standalone file on disk, so the in-memory rasters
# (which still reference tiles in terra_tmp) are no longer needed. Drop them
# and force garbage collection so no SpatRaster holds an open handle, then
# delete the scratch tiles. We delete the *contents* of terra_tmp rather than
# the folder itself, so terra's tempdir stays valid if the session continues.
rm(clim, clim_layers, year_stacks)
gc()
terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
unlink(list.files(terra_tmp, full.names = TRUE, recursive = TRUE), recursive = TRUE)
message("Cleared terra scratch in ", terra_tmp)
