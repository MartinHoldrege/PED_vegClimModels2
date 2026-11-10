# //////////////////////////////////////////////////////////////////////////
#
# Aggregate daily gridMET (pr, tmmx, tmmn) to monthly values, 1991-2020,
# on the native gridMET ~4 km grid (same grid family as MACAv2-METDATA).
#
# Output is one 12-layer GeoTIFF (Jan..Dec) per variable-year, in the same
# form as the Daymet monthly files read by 01_summarise_daymet_climate_data.R:
#   gridmet_prcp_monttl_<year>.tif  (mm, monthly total)
#   gridmet_tmax_monavg_<year>.tif  (deg C, monthly mean of daily max)
#   gridmet_tmin_monavg_<year>.tif  (deg C, monthly mean of daily min)
# These are read by 01_summarise_gridmet_climate_data.R.
#
# Input: daily gridMET NetCDFs (pr_<year>.nc, tmmx_<year>.nc, tmmn_<year>.nc;
# temperatures in K, precip in mm) in a colleague's folder next to this repo
# on Hovenweep (see gridmet_daily_dir).
#
# Days with NA give an NA month (na.rm = FALSE), so cells outside CONUS land
# stay NA. Re-running skips files that already exist.
#
# Run from the repository root. On the HPC run via hpc/summarise_gridmet.sh
#
# started November 2026
# //////////////////////////////////////////////////////////////////////////

source("config/paths.R") # provides paths
library(terra)

# Cap terra's RAM use (GB); keep well below --mem in hpc/summarise_gridmet.sh
terra::terraOptions(memmax = 8)


# Parameters --------------------------------------------------------------

years <- 1991:2020
rerun <- FALSE # if TRUE, overwrite existing monthly files

# Daily gridMET folder, relative to the repository root (the gridMET folder
# sits next to PED_vegClimModels2 on Hovenweep)
gridmet_daily_dir <- "../gridMET/data-raw"

out_dir <- file.path(paths$large,
                     "Data_processed/WallToWallClimateData/gridmet/monthly")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# gridMET variable -> output file tag, aggregation, expected units, and
# offset added to the monthly value (K -> deg C for temperatures)
gridmet_vars <- list(
  pr   = list(tag = "prcp_monttl", fun = "sum",  unit = "mm", offset = 0),
  tmmx = list(tag = "tmax_monavg", fun = "mean", unit = "K",  offset = -273.15),
  tmmn = list(tag = "tmin_monavg", fun = "mean", unit = "K",  offset = -273.15)
)

# Fail early if any input is missing
if (!dir.exists(gridmet_daily_dir)) {
  stop("gridMET folder not found: ", normalizePath(gridmet_daily_dir, mustWork = FALSE))
}
in_files <- expand.grid(var = names(gridmet_vars), yr = years,
                        stringsAsFactors = FALSE)
in_files$path <- file.path(gridmet_daily_dir,
                           paste0(in_files$var, "_", in_files$yr, ".nc"))
missing <- in_files$path[!file.exists(in_files$path)]
if (length(missing) > 0) {
  stop("Missing gridMET files:\n", paste(missing, collapse = "\n"))
}

# Every file is checked against the grid of the first one
ref_grid <- terra::rast(in_files$path[1])


# Aggregate one variable-year ---------------------------------------------

#' Aggregate one year of daily gridMET to monthly values
#'
#' @param var gridMET variable ("pr", "tmmx" or "tmmn").
#' @param yr Year (integer).
#' @return 12-layer SpatRaster (Jan..Dec), layers named "YYYY-MM", in mm
#'   (precip) or deg C (temperature).
#' @examples
#' \dontrun{
#' r <- daily_to_monthly("tmmx", 2020)
#' terra::global(r[[7]], "range", na.rm = TRUE) # July, deg C
#' }
daily_to_monthly <- function(var, yr) {
  spec <- gridmet_vars[[var]]
  r <- terra::rast(file.path(gridmet_daily_dir, paste0(var, "_", yr, ".nc")))

  # checks: same grid, expected units, one layer per day of the year
  terra::compareGeom(r, ref_grid, stopOnError = TRUE)
  stopifnot(all(terra::units(r) == spec$unit))
  dates <- as.Date(terra::time(r))
  all_days <- seq(as.Date(paste0(yr, "-01-01")), as.Date(paste0(yr, "-12-31")),
                  by = "day")
  # compare values (not identical(), which also compares storage type)
  if (length(dates) != length(all_days) || !isTRUE(all(dates == all_days))) {
    stop(var, " ", yr, ": layer dates (terra::time) are not one per day of ",
         "the year; got ", length(dates), " layers, ", sum(is.na(dates)),
         " NA dates, range ", paste(range(dates, na.rm = TRUE), collapse = " to "))
  }

  month <- as.integer(format(dates, "%m"))
  out <- terra::tapp(r, index = month, fun = spec$fun) + spec$offset
  names(out) <- paste0(yr, "-", sprintf("%02d", 1:12))
  out
}


# Loop over variable-years ------------------------------------------------

for (i in seq_len(nrow(in_files))) {
  var <- in_files$var[i]
  yr <- in_files$yr[i]
  out_file <- file.path(out_dir, paste0("gridmet_", gridmet_vars[[var]]$tag,
                                        "_", yr, ".tif"))
  if (file.exists(out_file) && !rerun) next
  message(var, " ", yr)

  out <- daily_to_monthly(var, yr)

  # Catches unapplied scale factors or a missed K -> deg C conversion
  rng <- unlist(terra::global(out, "range", na.rm = TRUE))
  ok <- if (var == "pr") all(rng >= 0 & rng < 3000) else all(rng > -60 & rng < 60)
  if (!ok) stop(var, " ", yr, ": implausible monthly values (",
                paste(round(range(rng), 1), collapse = " to "), ")")

  # write to a temporary name first so an interrupted run never leaves a
  # partial file under the final name
  tmp <- sub("\\.tif$", "_partial.tif", out_file)
  terra::writeRaster(out, tmp, overwrite = TRUE)
  file.rename(tmp, out_file)
}

message("Monthly gridMET files are in ", out_dir)
