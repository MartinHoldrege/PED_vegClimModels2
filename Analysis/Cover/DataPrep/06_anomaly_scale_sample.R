# 06_anomaly_scale_sample.R
#
# A CONUS-wide random sample of cell-years, used only to compute the scaling
# constants (sd) of the climate anomaly predictors (_3yrAnom) in
# 09_compute_scale_params.R. Like the climate normals, which are scaled with
# all CONUS cells, the sample isn't masked by LCMAP or fire, so the constants
# don't depend on how the training data are sampled. Made once (fixed seed),
# and not changed when the training data are rebuilt.
#
# n_cells random cells of the snap grid, each with every year in `years` (the
# RAP sample's years, all with full 30-year normals).
#
# Input:
#   daymet_conus_snap_1000m.tif - 00_create_snap_raster.R (read_mask())
# Output:
#   cover_cell_years_anom-sample.csv (cell, year, x, y), read by
#   07_summarise_daymet_climate_points.R with vc = "anom-sample"
#
# Note number if this script isn't strictly 'correct'.
# it is numbered '06_...' b/ it creates a file for '07_'
# but it itself doesn't have an '05_...' dependency. 
# November 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

n_cells <- 2e4
years   <- 2010:2024
seed    <- 1234

cov_dir <- file.path(paths$large, "Data_processed/cover_combined")

# sample ------------------------------------------------------------------

snap <- read_mask()
pool <- which(!is.na(terra::values(snap, mat = FALSE)))

set.seed(seed)
cells <- sort(sample(pool, n_cells))
xy <- terra::xyFromCell(snap, cells)

# every cell crossed with every year
sample_cell_years <- tibble(cell = cells, x = xy[, "x"], y = xy[, "y"]) |>
  expand_grid(year = years) |>
  select(cell, year, x, y)

stopifnot(nrow(sample_cell_years) == n_cells * length(years),
          !anyDuplicated(sample_cell_years[c("cell", "year")]))

write_csv(sample_cell_years,
          file.path(cov_dir, "cover_cell_years_anom-sample.csv"))
