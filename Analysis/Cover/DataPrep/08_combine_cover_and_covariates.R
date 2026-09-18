# 08_combine_cover_and_covariates.R
#
# Joins the cover pixel-years (field + RAP) to the climate at each pixel-year
# and the soil covariates at each pixel, giving one table for the cover models.
# Joins are on cell (and year), with checks that x, y agree with the cell.
# Rows missing soils are kept.
#
# Inputs:
#   cover_by_pixel_year_all-sources_<vc>.csv - 06_cover_add-rap.R
#   daymet_climate-at-points_<vc>.csv        - 07_summarise_daymet_climate_points.R
#   soil_covariates_solus100_1000m.tif       - 02_soils_calculate_variables.R
#   daymet_conus_snap_1000m.tif              - 00_create_snap_raster.R
#                                              (read_mask())
#
# Output:
#   cover_clim_soils_<vc>.csv (read with read_cover_training())
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

vc <- opt$vc  # cover data version

# params ------------------------------------------------------------------

cov_dir  <- file.path(paths$large, "Data_processed/cover_combined")
soil_dir <- file.path(paths$large, "Data_processed/soils")
out_file <- file.path(cov_dir, paste0("cover_clim_soils_", vc, ".csv"))

# read --------------------------------------------------------------------

cover <- read_csv(
  file.path(cov_dir, paste0("cover_by_pixel_year_all-sources_", vc, ".csv")),
  show_col_types = FALSE)

clim <- read_csv(
  file.path(cov_dir, paste0("daymet_climate-at-points_", vc, ".csv")),
  show_col_types = FALSE)

mask_r <- read_mask()
soil_r <- terra::rast(file.path(soil_dir, "soil_covariates_solus100_1000m.tif"))

# checks ------------------------------------------------------------------

stopifnot(
  !anyDuplicated(cover[c("cell", "year")]),
  !anyDuplicated(clim[c("cell", "year")]),
  # x, y lie in the cell they're labelled with
  all(terra::cellFromXY(mask_r, as.matrix(cover[c("x", "y")])) == cover$cell),
  # 07 is built from the 06 file; a mismatch means one is stale (rerun 07)
  nrow(clim) == nrow(cover),
  nrow(anti_join(cover, clim, by = c("cell", "year"))) == 0,
  # soils are looked up by cell, so the grids must match
  terra::compareGeom(soil_r, mask_r)
)

# combine -----------------------------------------------------------------

out <- cover |>
  left_join(clim, by = c("cell", "year"), suffix = c("", "_clim"))

# joined on cell; the x, y in each file must agree too
stopifnot(all(out$x == out$x_clim), all(out$y == out$y_clim))

out <- out |>
  bind_cols(as_tibble(soil_r[out$cell]))

write_csv(out, out_file)
