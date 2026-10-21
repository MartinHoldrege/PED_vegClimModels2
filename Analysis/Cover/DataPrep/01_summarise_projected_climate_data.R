# //////////////////////////////////////////////////////////////////////////
# 07_GetForecastedClimateData.R
#
# Derive end-of-century climate "normals" (the _CLIM variable set) from MACA
# climate projections, as a multi-band raster on the Daymet 1 km grid, for use
# as inputs to the vegetation-climate models.
#
# Run once per climate model (GCM), scenario, and period. Periods match the
# Climate Toolbox: historical 1971-2000; rcp45 and rcp85 2040-2069 and
# 2070-2099. From the command line:
#   Rscript 01_summarise_projected_climate_data.R --gcm=<gcm> \
#     --scenario=<scenario> --period=<period>
# e.g. --gcm=BNU-ESM --scenario=rcp85 --period=2070-2099. The options are
# defined in Functions/params.R. On the HPC, slurm/summarise_maca.sh
# runs all 20 GCMs x 5 scenario-periods as a job array. Run interactively, the
# defaults below are used.
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

# //////////////////////////////////////////////////////////////////////////

library(terra)
library(stringr)

source("Functions/init.R")          # provides paths$large, paths$large0, daymet_grid
source_functions()         # sources climate_functions.R

# Cap terra's RAM use (GB) so it writes to temp files instead of exceeding the
# Slurm memory request; terra otherwise sizes itself from the node's total RAM.
# Keep this well below --mem in slurm/summarise_maca.sh.
terra::terraOptions(memmax = 8)


# Parameters --------------------------------------------------------------

# Passed as --gcm, --scenario, --period (defined in Functions/params.R);
# defaults below are used when run interactively.
climate_model <- if (is.null(opt$gcm)) "IPSL-CM5A-MR" else opt$gcm
scenario <- if (is.null(opt$scenario)) "rcp85" else opt$scenario
period <- if (is.null(opt$period)) "2070-2099" else opt$period
rerun <- FALSE

# Climate-normal windows (matching the Climate Toolbox), by scenario.
windows <- list(
  historical = "1971-2000",
  rcp45 = c("2040-2069", "2070-2099"),
  rcp85 = c("2040-2069", "2070-2099")
)
stopifnot(scenario %in% names(windows), period %in% windows[[scenario]])
year_start <- as.integer(str_sub(period, 1, 4))
year_end   <- as.integer(str_sub(period, 6, 9))

# MACA files come in time chunks; read every chunk that overlaps the window,
# then subset to the window.
all_chunks <- if (scenario == "historical") {
  c("1950_1969", "1970_1989", "1990_2005")
} else {
  c("2006_2025", "2026_2045", "2046_2065", "2066_2085", "2086_2099")
}
chunk_start <- as.integer(str_sub(all_chunks, 1, 4))
chunk_end   <- as.integer(str_sub(all_chunks, 6, 9))
maca_chunks <- all_chunks[chunk_start <= year_end & chunk_end >= year_start]

message("GCM: ", climate_model, "; scenario: ", scenario,
        "; years: ", year_start, "-", year_end,
        "; MACA chunks: ", paste(maca_chunks, collapse = ", "))

# Directory holding the raw MACA NetCDFs (a subdirectory of paths$large0).
maca_dir <- file.path(paths$large, "Data_raw/macaClimateProjections/Data")

# MACA variable name -> our short name. tasmin/tasmax are in Kelvin; pr in mm.
maca_vars <- c(tasmin = "tmin", tasmax = "tmax", pr = "prcp")

# Where per-year intermediate rasters are written (under paths$large).
intermediate_dir <- file.path(paths$large, 'Data_processed/WallToWallClimateData', 
                              "ForecastedClimate_intermediate",
                              paste0(climate_model, "_", scenario))

out_dir <-  file.path(paths$large, 'Data_processed/WallToWallClimateData')
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)


# Helper: build a MACA file path -----------------------------------------

#' Construct the path to a MACA monthly NetCDF for a variable and time chunk.
#'
#' @param maca_var MACA variable string (e.g. "tasmin", "tasmax", "pr").
#' @param chunk Time-chunk string (e.g. "2066_2085").
#' The ensemble member is matched by pattern because it is not always r1i1p1
#' (CCSM4 is r6i1p1).
#'
#' @return Full file path under maca_dir.
#' @examples
#' \dontrun{
#' climate_model <- "IPSL-CM5A-MR"; scenario <- "rcp85"
#' maca_path("pr", "2066_2085")
#' }
maca_path <- function(maca_var, chunk) {
  pattern <- paste0("^macav2livneh_", maca_var, "_", climate_model,
                    "_r\\d+i\\d+p\\d+_", scenario, "_", chunk,
                    "_CONUS_monthly\\.nc$")
  f <- list.files(maca_dir, pattern = pattern, full.names = TRUE)
  if (length(f) != 1) {
    stop("Expected 1 file matching ", pattern, ", found ", length(f))
  }
  f
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
  if (file.exists(out_file) & !rerun) {
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


# Map: per-year metric name -> list(output _CLIM name, reduction function).
# Note durationFrostFreeDays and a few others appear in BOTH a mean and a
# percentile output, matching the training variable set.


message("Loading intermediate rasters and reducing across years ...")
year_files <- file.path(intermediate_dir, paste0("annualMetrics_", years, ".tif"))
year_stacks <- lapply(year_files, terra::rast)

# For a given per-year metric, build a multi-year stack (one layer per year).
stack_metric <- function(metric_name) {
  layers <- lapply(year_stacks, function(s) s[[metric_name]])
  do.call(c, layers)
}

clim_layers <- lapply(climate_reductions, function(spec) {
  metric_name <- spec[[1]]
  out_name    <- spec[[2]]
  fun         <- .reduction_raster[[spec[[3]]]]
  r <- fun(stack_metric(metric_name))
  names(r) <- out_name
  r
})

clim <- do.call(c, clim_layers)


# * project to daymet grid --------------------------------------------------

clim <- terra::project(clim, read_mask())   # bilinear (default for continuous)

# 4. Write final per-model raster -----------------------------------------

out_file <- file.path(out_dir,
                      paste0("maca_", climate_model, "_",
                             scenario, "_", period, "_CLIM.tif"))
terra::writeRaster(clim, out_file, overwrite = TRUE)
message("Wrote ", out_file)
