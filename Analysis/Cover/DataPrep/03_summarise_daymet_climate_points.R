# //////////////////////////////////////////////////////////////////////////
# 02_summarise_daymet_climate_points.R
#
# Point-location analogue of 01_summarise_daymet_climate_data.R.
#
# For each observation (a location and a year), calculate the same climate
# variables as the gridded product, over two trailing windows:
#   _CLIM : the 30 years preceding the observation year (truncated where the
#           Daymet record does not reach back that far)
#   _3yr  : the 3 years preceding the observation year
# Neither window includes the observation year itself.
#
# This is done at points rather than wall-to-wall because a gridded product
# would need a separate 30-year normal for every observation year.
#
# Metrics are calculated with calc_annual_metrics_df() and reduced with
# climate_reductions, so the definitions match the gridded pipeline exactly.
#
# started August 2026
# //////////////////////////////////////////////////////////////////////////

library(terra)
library(sf)
library(tidyverse)

source("Functions/init.R")
source_functions()


# Parameters --------------------------------------------------------------

rerun <-  TRUE # recreate intermediate files

# Daymet years available on disk. Windows are truncated to this range.
daymet_years <- 1980:2023

n_years_clim <- 30
n_years_lag  <- 3

# Shortest acceptable trailing window. Early observation years cannot have a
# full 30 years of preceding Daymet data (e.g. 2000 has only 1980-1999).
min_years_clim <- 20
min_years_lag  <- 3

daymet_dir <- file.path(paths$large, "Data_raw/daymet/rawMonthlyData")

intermediate_dir <- file.path(paths$large, "Data_processed/CoverData",
                              "DaymetPoints_intermediate")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(paths$large, "Data_processed/CoverData",
                      "daymetClimateAtPoints.rds")


# Input observations ------------------------------------------------------
# DUMMY DATA -- replace with the real observation sf object. Only two things
# are required: a `year` column and point geometry.

set.seed(1)
obs_sf <- tibble(
  year = sample(2000:2024, 500, replace = TRUE),
  lon  = runif(500, -120, -75),
  lat  = runif(500, 30, 48)
) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

stopifnot("year" %in% names(obs_sf), inherits(obs_sf, "sf"))


# Locate Daymet files -----------------------------------------------------

daymet_tifs <- list.files(daymet_dir, pattern = "\\.tif$",
                          recursive = TRUE, full.names = TRUE)

#' Path to the Daymet monthly GeoTIFF for one variable and year
#'
#' @param var_tag Daymet file tag (e.g. "prcp_monttl").
#' @param yr Year (integer).
#' @return Single file path.
daymet_path <- function(var_tag, yr) {
  hits <- str_subset(daymet_tifs,
                     paste0("daymet_v4_", var_tag, "_na_", yr, "\\.tif$"))
  if (length(hits) != 1) {
    stop("Expected 1 file for ", var_tag, " ", yr, ", found ", length(hits))
  }
  hits
}

var_tags <- c(tmin = "tmin_monavg", tmax = "tmax_monavg",
              prcp = "prcp_monttl",  vp = "vp_monavg")

# Fail early if a year is missing.
invisible(lapply(daymet_years, function(yr) lapply(var_tags, daymet_path, yr = yr)))


# Assign observations to Daymet cells -------------------------------------
# Cell identity comes from the CONUS mask returned by read_mask(), so the `cell`
# column is directly comparable to the gridded product and can be used to join
# back to it later. Extraction happens once per unique cell, not once per
# observation: many observations share a cell, and repeat visits share a cell
# across years.

mask_r <- read_mask()

# cellFromXY() returns the cell the point falls in, or NA if the point is
# outside the raster. There is no nearest-cell search: a point off the grid
# stays NA rather than borrowing a distant cell's values.
obs_sf <- obs_sf |>
  st_transform(crs(mask_r)) |>
  mutate(cell = cellFromXY(mask_r, st_coordinates(geometry)))

n_off_grid <- sum(is.na(obs_sf$cell))
if (n_off_grid > 0) {
  message(n_off_grid, " observations fall outside the mask extent; dropping")
  obs_sf <- filter(obs_sf, !is.na(cell))
}

# read_mask() holds the cell number at unmasked cells and NA elsewhere, so the
# value is both the CONUS test and the identifier. Confirm the stored value
# really is the terra cell number; if it is a different index the join key
# below would be wrong.
mask_vals <- mask_r[obs_sf$cell][[1]]
outside <- is.na(mask_vals)
if (any(outside)) {
  message(sum(outside), " observations fall outside the CONUS mask; dropping")
  obs_sf <- filter(obs_sf, !outside)
  mask_vals <- mask_vals[!outside]
}
stopifnot(all(mask_vals == obs_sf$cell))

cells <- sort(unique(obs_sf$cell))

if (any(obs_sf$year - 1 > max(daymet_years))) {
  stop("Observation years require Daymet data past ", max(daymet_years), ": ",
       paste(sort(unique(obs_sf$year[obs_sf$year - 1 > max(daymet_years)])),
             collapse = ", "))
}

# Years of Daymet needed to cover the longest window for any observation.
years_needed <- seq(max(min(daymet_years), min(obs_sf$year) - n_years_clim),
                    min(max(daymet_years), max(obs_sf$year) - 1))


#' Read monthly Daymet values at cells of the mask grid
#'
#' Crops the North America file to the mask so that cell numbers refer to the
#' mask grid, and errors if the two do not align cell-for-cell.
#'
#' @param var_tag Daymet file tag (e.g. "tmin_monavg").
#' @param yr Year (integer).
#' @param cells Integer vector of mask cell numbers.
#' @return Numeric matrix, length(cells) x 12.
read_daymet_cells <- function(var_tag, yr, cells) {
  r <- crop(rast(daymet_path(var_tag, yr)), mask_r)
  if (!compareGeom(r, mask_r, crs = TRUE, res = TRUE, ext = TRUE,
                   rowcol = TRUE, stopOnError = FALSE)) {
    stop("Grid mismatch for ", var_tag, " ", yr, ": does not match the mask grid")
  }
  m <- as.matrix(r[cells])
  stopifnot(ncol(m) == 12,
            nrow(m) == length(cells))
  m
}

# Cells inside the mask but with no Daymet data hold NA in every year.
# Identify them once and drop those observations here, so that no window is
# silently built from fewer years than it appears to use.
ref_vals <- read_daymet_cells("tmin_monavg", daymet_years[1], cells)
cells_no_data <- cells[is.na(ref_vals[, 1])]

# given the mask being used is based on daymet, this shouldn't yield 
# any cells with no data
if (length(cells_no_data) > 0) {
  n_drop <- sum(obs_sf$cell %in% cells_no_data)
  message(n_drop, " observations at ", length(cells_no_data),
          " cells with no Daymet data; dropping")
  obs_sf <- filter(obs_sf, !cell %in% cells_no_data)
  cells <- sort(unique(obs_sf$cell))
}
rm(ref_vals)

message(nrow(obs_sf), " observations at ", length(cells), " unique Daymet cells")



# 1. Per-year annual metrics at each cell ---------------------------------
# Written per year so a crash mid-run can be resumed. Only the 18 annual
# metrics are kept, not the 48 monthly values they come from.

for (yr in years_needed) {
  yr_file <- file.path(intermediate_dir, paste0("annualMetrics_", yr, ".rds"))
  if (file.exists(yr_file) & !rerun) next
  message("  extracting and summarising ", yr)

  monthly <- map(var_tags, read_daymet_cells, yr = yr, cells = cells)

  metrics <- calc_annual_metrics_df(tmin12 = monthly$tmin,
                                    tmax12 = monthly$tmax,
                                    prcp12 = monthly$prcp,
                                    vp12   = monthly$vp)
  metrics$cell <- cells
  metrics$year <- yr

  saveRDS(metrics, yr_file)
  rm(monthly, metrics)
  gc()
}

annual <- map_dfr(years_needed, function(yr) {
  readRDS(file.path(intermediate_dir, paste0("annualMetrics_", yr, ".rds")))
})

# No-data cells were dropped above, so any NA remaining here is unexpected
# (e.g. a year with incomplete coverage) and would quietly shorten a window.
if (anyNA(annual$tmin_annAvg)) {
  stop("Unexpected NA in annual metrics for ",
       length(unique(annual$cell[is.na(annual$tmin_annAvg)])), " cells")
}


# 2. Trailing-window reductions -------------------------------------------

targets <- obs_sf |>
  st_drop_geometry() |>
  distinct(cell, year)

message("Reducing over the ", n_years_clim, "-year window ...")
clim <- roll_point_normals(annual, targets, n_years = n_years_clim,
                           out_suffix = "_CLIM", min_years = min_years_clim)

message("Reducing over the ", n_years_lag, "-year window ...")
lag3 <- roll_point_normals(annual, targets, n_years = n_years_lag,
                           out_suffix = "_3yr", min_years = min_years_lag)

clim <- rename(clim, n_years_CLIM = n_years_used)
lag3 <- rename(lag3, n_years_3yr  = n_years_used)


# 3. Anomalies and join back to the observations --------------------------

# Cell-centre coordinates (Daymet LCC), stored alongside `cell` so the cell
# identifier stays interpretable if the snap raster is ever rebuilt.
cell_xy <- xyFromCell(mask_r, cells) |>
  as_tibble() |>
  rename(cell_x = x, cell_y = y) |>
  mutate(cell = cells)

out <- obs_sf |>
  left_join(cell_xy, by = "cell") |>
  left_join(clim, by = c("cell", "year")) |>
  left_join(lag3, by = c("cell", "year"))

out <- bind_cols(out, calc_anomalies(st_drop_geometry(out), short_suffix = "_3yr"))

n_missing <- sum(is.na(out$tmean_meanAnnAvg_CLIM))
if (n_missing > 0) {
  message(n_missing, " observations have no ", n_years_clim,
          "-year normal (window shorter than min_years_clim = ",
          min_years_clim, ")")
}

saveRDS(out, out_file)
message("Wrote ", out_file)
