# //////////////////////////////////////////////////////////////////////////
#
# Derive 1991-2020 climate "normals" (the _CLIM variable set) from monthly
# gridMET, on the native gridMET ~4 km grid, using the same metric
# definitions and across-year reductions as
# 01_summarise_daymet_climate_data.R and 01_summarise_projected_climate_data.R.
#
# VPD uses the tmin dew-point approximation (vp12 = NULL in
# calc_annual_metrics()), as for MACA, so gridMET and MACA normals are
# directly comparable.
#
# Input: monthly gridMET GeoTIFFs from 00_gridmet_daily_to_monthly.R, in
#   paths$large/Data_processed/WallToWallClimateData/gridmet/monthly
#
# Output:
#   paths$large/Data_processed/WallToWallClimateData/gridmet/
#     gridmet_1991-2020_CLIM.tif
#
# Run from the repository root. On the HPC run via hpc/summarise_gridmet.sh
#
# started November 2026
# //////////////////////////////////////////////////////////////////////////

library(terra)
library(stringr)

source("Functions/init.R")
source_functions() # Functions/data/climate.R contains the key functions used here

# Cap terra's RAM use (GB); keep well below --mem in hpc/summarise_gridmet.sh
terra::terraOptions(memmax = 8)


# Parameters --------------------------------------------------------------

year_start <- 1991
year_end   <- 2020
years      <- year_start:year_end
rerun <- FALSE # if TRUE, recompute per-year intermediates

gridmet_dir <- file.path(paths$large,
                         "Data_processed/WallToWallClimateData/gridmet")
monthly_dir <- file.path(gridmet_dir, "monthly")
intermediate_dir <- file.path(paths$large,
                              "Data_processed/WallToWallClimateData",
                              "GridmetClimate_intermediate")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)


# Helper: read one variable-year ------------------------------------------

#' Read a 12-layer monthly gridMET GeoTIFF
#'
#' @param var_tag File tag ("prcp_monttl", "tmax_monavg" or "tmin_monavg").
#' @param yr Year (integer).
#' @return 12-layer SpatRaster (Jan..Dec).
read_gridmet_year <- function(var_tag, yr) {
  f <- file.path(monthly_dir, paste0("gridmet_", var_tag, "_", yr, ".tif"))
  if (!file.exists(f)) stop("Missing ", f, "; run 00_gridmet_daily_to_monthly.R")
  r <- terra::rast(f)
  stopifnot(terra::nlyr(r) == 12)
  r
}


# 1. Per-year annual metrics ----------------------------------------------
# Re-running skips finished years.

for (yr in years) {
  out_file <- file.path(intermediate_dir, paste0("annualMetrics_", yr, ".tif"))
  if (file.exists(out_file) && !rerun) {
    message("  year ", yr, " already done, skipping")
    next
  }
  message("  computing annual metrics for ", yr)

  tmin12 <- read_gridmet_year("tmin_monavg", yr)
  tmax12 <- read_gridmet_year("tmax_monavg", yr)
  prcp12 <- read_gridmet_year("prcp_monttl", yr)
  terra::compareGeom(tmin12, tmax12, prcp12, stopOnError = TRUE)

  metrics <- calc_annual_metrics(tmin12 = tmin12, tmax12 = tmax12,
                                 prcp12 = prcp12)
  terra::writeRaster(metrics, out_file, overwrite = TRUE)
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
  fun  <- .reduction_raster[[spec[[3]]]]
  r <- fun(stack_metric(metric_name))
  names(r) <- out_name
  r
})

clim <- do.call(c, clim_layers)
# Left on the native gridMET grid; cells outside CONUS land are already NA.


# 3. Write final raster ---------------------------------------------------

out_file <- file.path(gridmet_dir, paste0("gridmet_", year_start, "-",
                                          year_end, "_CLIM.tif"))
terra::writeRaster(clim, out_file, overwrite = TRUE)
message("Wrote ", out_file)
