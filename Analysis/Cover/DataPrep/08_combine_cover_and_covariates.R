# 08_combine_cover_and_covariates.R
#
# Joins the cover pixel-years (field + RAP) to the climate at each pixel-year
# and the soil covariates at each pixel, giving one table for the cover models.
# Joins are on cell (and year), with checks that x, y agree with the cell.
# Rows missing soils are kept.
#
# c_needle and c_broad (needleleaf and broadleaf tree cover) are replaced by
# the field fractions (frac_needle, frac_broad) times 2021 RAP tree cover at
# the cell, set to 0 where the RAP zero-tree layer is 1. RAP is treated as the
# better estimate of tree cover over the 1 km cell, and it is what the tree
# cover models are fit to, so these are comparable with the model's
# needleleaf share x predicted tree cover. 0 where the RAP tree cover is 0,
# even where the field fraction is NA. NA where the cell has no 2021 RAP
# tree cover (e.g. masked by LCMAP or fire in that export). Note: c_needle and
# c_broad are not used for model fitting just in validation
#
# Inputs:
#   cover_by_pixel_year_all-sources_<vc>.csv - 06_cover_add-rap.R
#   daymet_climate-at-points_<vc>.csv        - 07_summarise_daymet_climate_points.R
#   soil_covariates_solus100_1000m.tif       - 02_soils_calculate_variables.R
#   daymet_conus_snap_1000m.tif              - 00_create_snap_raster.R
#                                              (read_mask())
#   RAP_v3_cover_2021_1000m.tif              - 03_rap_sample.js
#   RAP_v3_fracZeroTree_lt3_2021-2021_1000m.tif
#                                            - 03_rap_frac-zero-tree.js
#                                              (read_zero_tree_raster())
#
# Output:
#   cover_clim_soils_<vc>.csv (read with read_cover_training())
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

vc <- opt$vc

if(vc == 'c02') {
  stop("cover data version 'c02' is created in 
       Analysis/Cover/DataPrep/05_rap_training_sample.R
       that version can't be created in this script")
} else if (!vc %in% c('c01', 'c03')) {
  stop('code not setup for other vc yet')
}

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
  select(-x_clim, -y_clim) |>
  bind_cols(as_tibble(soil_r[out$cell]))

# needleleaf and broadleaf cover from RAP tree cover ----------------------

rap_tree <- terra::rast(file.path(paths$large, "Data_processed/CoverData/rap",
                                  "RAP_v3_cover_2021_1000m.tif"))[["tree"]] |>
  # GEE exports come back one cell larger in each direction than the snap grid
  align_raster(mask_r)
zero_tree <- read_zero_tree_raster("zeroTree") |>
  align_raster(mask_r)

# RAP tree cover, 0 where the zero-tree layer is 1 (left as is where that
# layer is NA)
rap_tree_cell <- rap_tree[out$cell][[1]]
zero_cell <- zero_tree[out$cell][[1]]
rap_tree_cell[!is.na(zero_cell) & zero_cell == 1] <- 0

out <- out |>
  mutate(c_needle_field = c_needle,
         c_broad_field = c_broad,
         # 0 where RAP has no trees, even if the field fraction is NA (no
         # trees in the plot); this adds observations only where RAP is 0, so
         # the sample isn't random, but these columns are only for validation
         c_needle = if_else(!is.na(rap_tree_cell) & rap_tree_cell == 0, 0,
                            frac_needle * rap_tree_cell),
         c_broad  = if_else(!is.na(rap_tree_cell) & rap_tree_cell == 0, 0,
                            frac_broad * rap_tree_cell))

message("rows with frac_needle: ", sum(!is.na(out$frac_needle)),
        " | with c_needle (RAP tree cover available): ",
        sum(!is.na(out$c_needle)))

write_csv(out, out_file)
