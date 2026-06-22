# //////////////////////////////////////////////////////////////////////////
# 07_GetForecastedClimateData.R
#
# Derive end-of-century climate "normals" (the _CLIM variable set) from MACA
# climate projections, as a multi-band raster on the Daymet 1 km grid, for use
# as inputs to the vegetation-climate models.
#
# Run once per climate model by setting `climate_model` below. Two models are
# expected: "BNU-ESM" (cool/wet-ish) and "IPSL-CM5A-MR" (warm/dry), both
# RCP 8.5, end of century.
#
# Anomalies and the 3-year window are intentionally NOT calculated: anomalies
# are treated as 0 at prediction time. Soils are added downstream and are not
# duplicated here.
#
# Assumed environment (initialise the repo before running):
#   - source("init.R") has been run, providing:
#       * paths$large   : directory for large outputs and per-year intermediates
#       * paths$large0   : directory holding the raw MACA NetCDF files
#       * daymet_grid    : SpatRaster used as the snap/target grid (1 km LCC)
#   - source_functions() (or equivalent) has sourced climate_functions.R, providing:
#       weighted_annual_mean(), calc_annual_metrics(), and the metric helpers.
#   - packages: terra, stringr (via tidyverse), and the repo's standard set.
#
# Consistency note (weighted means / VPD mean):
#   This script uses the mathematically correct day-length denominator
#   (sum of month weights, 365.5) via the default of weighted_annual_mean().
#   The Daymet training pipeline effectively used a denominator of 372
#   (= 31 * 12), because it divided each weight by 31 and then called mean()
#   (which divides by 12). This makes the four affected inputs
#   (tmin/tmax/tmean annual means and annVPD_mean) ~1.8% larger in training
#   than here (372 / 365.5 ≈ 1.018). This is a known, accepted offset for
#   the current fitted model. To reproduce training exactly, pass
#   denom = 31 * 12 (i.e. 372) to calc_annual_metrics() below.

# VPD--currently vpd isn't calculated correctly (for consistency with other dataset)
# it is closely related to saturation vapor pressure 
# //////////////////////////////////////////////////////////////////////////

library(terra)
library(stringr)

source("Functions/init.R")          # provides paths$large, paths$large0, daymet_grid
source_functions()         # sources climate_functions.R


# Parameters --------------------------------------------------------------

# Which climate model to process: "BNU-ESM" or "IPSL-CM5A-MR".
climate_model <- "IPSL-CM5A-MR" # "BNU-ESM"

rcp        <- "rcp85"
year_start <- 2068        # first year of the climate-normal window
year_end   <- 2099        # last year (full available span: 32 years)

# MACA files come in two time chunks; we read both and subset to the window.
maca_chunks <- c("2066_2085", "2086_2099")

# Directory holding the raw MACA NetCDFs (a subdirectory of paths$large0).
maca_dir <- file.path(paths$large0, "Data_raw/macaClimateProjections/Data")

# MACA variable name -> our short name. tasmin/tasmax are in Kelvin; pr in mm.
maca_vars <- c(tasmin = "tmin", tasmax = "tmax", pr = "prcp")

# Where per-year intermediate rasters are written (under paths$large).
intermediate_dir <- file.path(paths$large, 'Data_processed/WallToWallClimateData', 
                              "ForecastedClimate_intermediate",
                              paste0(climate_model, "_", rcp))

out_dir <-  file.path(paths$large, 'Data_processed/WallToWallClimateData')
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)


# Helper: build a MACA file path -----------------------------------------

#' Construct the path to a MACA monthly NetCDF for a variable and time chunk.
#'
#' @param maca_var MACA variable string (e.g. "tasmin", "tasmax", "pr").
#' @param chunk Time-chunk string (e.g. "2066_2085").
#' @return Full file path under paths$large0.
maca_path <- function(maca_var, chunk) {
  fname <- paste0("macav2livneh_", maca_var, "_", climate_model,
                  "_r1i1p1_", rcp, "_", chunk, "_CONUS_monthly.nc")
  file.path(maca_dir, fname)
}


# 1. Read and prepare monthly stacks --------------------------------------

#' Read both MACA chunks for one variable, project to the Daymet grid, label
#' layers by date, and (for temperature) convert Kelvin to Celsius.
#'
#' @param maca_var MACA variable string ("tasmin", "tasmax", "pr").
#' @return SpatRaster of monthly layers on daymet_grid, named by date.
read_maca_var <- function(maca_var) {
  parts <- lapply(maca_chunks, function(chunk) {
    r <- terra::rast(maca_path(maca_var, chunk))
    names(r) <- terra::time(r)
    r
  })
  r <- do.call(c, parts)

  # Kelvin -> Celsius for temperature variables.
  if (maca_var %in% c("tasmin", "tasmax")) {
    r <- r - 273.15
  }
  r
}

message("Reading MACA data for ", climate_model, " ...")
tmin_all <- read_maca_var("tasmin")
tmax_all <- read_maca_var("tasmax")
prcp_all <- read_maca_var("pr")

# Year of each layer, used to subset to the window and to pull each year's
# 12 monthly layers.
layer_years <- as.integer(str_sub(names(tmin_all), 1, 4))
years <- year_start:year_end
stopifnot(all(years %in% layer_years))


# 2. Per-year annual metrics ----------------------------------------------
# For each year, pull its 12 monthly layers, compute annual metrics, and write
# the result to an intermediate file. Re-running can skip years already done.

for (yr in years) {
  out_file <- file.path(intermediate_dir, paste0("annualMetrics_", yr, ".tif"))
  if (file.exists(out_file)) {
    message("  year ", yr, " already done, skipping")
    next
  }
  message("  computing annual metrics for ", yr)

  idx <- which(layer_years == yr)
  stopifnot(length(idx) == 12)   # expect exactly 12 monthly layers per year

  metrics <- calc_annual_metrics(
    tmin12 = tmin_all[[idx]],
    tmax12 = tmax_all[[idx]],
    prcp12 = prcp_all[[idx]]
  )
  terra::writeRaster(metrics, out_file, overwrite = TRUE)
}


# 3. Reduce across years to the _CLIM normals -----------------------------
# Most metrics: mean across years. A few use percentiles, matching training.

# Reductions keyed by the per-year metric name. Each entry is a function taking
# a multi-year SpatRaster (one layer per year) and returning a single layer.
mean_fun  <- function(r) terra::app(r, fun = mean, na.rm = TRUE)
q95_fun   <- function(r) terra::app(r, fun = function(x) stats::quantile(x, 0.95, na.rm = TRUE))
q05_fun   <- function(r) terra::app(r, fun = function(x) stats::quantile(x, 0.05, na.rm = TRUE))

# Map: per-year metric name -> list(output _CLIM name, reduction function).
# Note durationFrostFreeDays and a few others appear in BOTH a mean and a
# percentile output, matching the training variable set.
reductions <- list(
  list("tmin_annAvg",            "tmin_meanAnnAvg_CLIM",                  mean_fun),
  list("tmax_annAvg",            "tmax_meanAnnAvg_CLIM",                  mean_fun),
  list("tmean",                  "tmean_meanAnnAvg_CLIM",                 mean_fun),
  list("totalAnnPrecip",         "prcp_meanAnnTotal_CLIM",                mean_fun),
  list("T_warmestMonth",         "T_warmestMonth_meanAnnAvg_CLIM",        mean_fun),
  list("T_coldestMonth",         "T_coldestMonth_meanAnnAvg_CLIM",        mean_fun),
  list("precip_wettestMonth",    "precip_wettestMonth_meanAnnAvg_CLIM",   mean_fun),
  list("precip_driestMonth",     "precip_driestMonth_meanAnnAvg_CLIM",    mean_fun),
  list("precip_Seasonality",     "precip_Seasonality_meanAnnAvg_CLIM",    mean_fun),
  list("PrecipTempCorr",         "PrecipTempCorr_meanAnnAvg_CLIM",        mean_fun),
  list("aboveFreezing_month",    "aboveFreezing_month_meanAnnAvg_CLIM",   mean_fun),
  list("isothermality",          "isothermality_meanAnnAvg_CLIM",         mean_fun),
  list("annWaterDeficit",        "annWaterDeficit_meanAnnAvg_CLIM",       mean_fun),
  list("annWetDegDays",          "annWetDegDays_meanAnnAvg_CLIM",         mean_fun),
  list("annVPD_mean",            "annVPD_mean_meanAnnAvg_CLIM",           mean_fun),
  list("annVPD_max",             "annVPD_max_meanAnnAvg_CLIM",            mean_fun),
  list("annVPD_min",             "annVPD_min_meanAnnAvg_CLIM",            mean_fun),
  list("annVPD_max",             "annVPD_max_95percentile_CLIM",          q95_fun),
  list("annWaterDeficit",        "annWaterDeficit_95percentile_CLIM",     q95_fun),
  list("annWetDegDays",          "annWetDegDays_5percentile_CLIM",        q05_fun),
  list("durationFrostFreeDays",  "durationFrostFreeDays_5percentile_CLIM", q05_fun),
  list("durationFrostFreeDays",  "durationFrostFreeDays_meanAnnAvg_CLIM", mean_fun)
)

message("Loading intermediate rasters and reducing across years ...")
year_files <- file.path(intermediate_dir, paste0("annualMetrics_", years, ".tif"))
year_stacks <- lapply(year_files, terra::rast)

# For a given per-year metric, build a multi-year stack (one layer per year).
stack_metric <- function(metric_name) {
  layers <- lapply(year_stacks, function(s) s[[metric_name]])
  do.call(c, layers)
}

clim_layers <- lapply(reductions, function(spec) {
  metric_name <- spec[[1]]
  out_name    <- spec[[2]]
  fun         <- spec[[3]]
  r <- fun(stack_metric(metric_name))
  names(r) <- out_name
  r
})

clim <- do.call(c, clim_layers)


# * project to daymet grid --------------------------------------------------

clim <- terra::project(clim, daymet_grid)   # bilinear (default for continuous)

# 4. Write final per-model raster -----------------------------------------

out_file <- file.path(out_dir,
                      paste0("ForecastedClimateData_", climate_model, "_", rcp, "_CLIM.tif"))
terra::writeRaster(clim, out_file, overwrite = TRUE)
message("Wrote ", out_file)
