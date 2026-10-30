# 08_combine_cover_and_covariates.R
#
# Joins the cover pixel-years of one cover version (c04, c05 or c06) to the
# climate at each pixel-year, the zero-tree status at each pixel-year and the
# soil covariates at each pixel, giving one table for the cover models. Joins
# are on cell (and year), with checks that x, y agree with the cell. Rows
# missing soils are kept.
#
# Zero tree (response of the zero-tree model): pct_zero_tree, the % of the
# cell's natural-land, unburned 30 m pixels with < 3% RAP tree cover that
# year, and zero_tree (1 where that is above 90%). NA only where a cell has
# no unmasked 30 m pixels that year.
#
# c_needle and c_broad (needleleaf and broadleaf tree cover) are replaced by
# the field fractions (frac_needle, frac_broad) times same-year RAP tree cover
# at the cell, set to 0 where zero_tree is 1. RAP is treated as the better
# estimate of tree cover over the 1 km cell, and it is what the tree cover
# models are fit to, so these are comparable with the model's needleleaf
# share x predicted tree cover. 0 where the RAP tree cover is 0, even where
# the field fraction is NA. NA where RAP has no tree cover for that cell-year,
# and NA in RAP rows (they have no field fractions). The original field
# values are kept as c_needle_field and c_broad_field. Note: c_needle and
# c_broad are not used for model fitting just in validation
#
# Inputs:
#   cover_by_pixel_year_<vc>.csv              - 06_build_cover_datasets.R
#   daymet_climate-at-points_c04-c06.csv      - 07_summarise_daymet_climate_points.R
#   RAP_v3_zeroTree_lt3_cell-years_c04-c06.csv
#                                             - 07_rap_zero-tree_cell-years.js,
#                                               via 04_download_gee_output.R
#   soil_covariates_solus100_1000m.tif        - 02_soils_calculate_variables.R
#   daymet_conus_snap_1000m.tif               - 00_create_snap_raster.R
#                                               (read_mask())
#   RAP_v3_cover-tree_2000-2024_1000m.tif     - 01_rap_cover.js
#   RAP_v3_fracZeroTree_lt3_2021-2021_1000m.tif
#                                             - 03_rap_frac-zero-tree.js
#                                               (read_zero_tree_raster(); only
#                                               to check the 2021 rows)
#
# Output:
#   cover_clim_soils_<vc>.csv (read with read_cover_training())
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

vc <- opt$vc

if (!vc %in% c("c04", "c05", "c06")) {
  stop("this script only builds c04-c06")
}

# params ------------------------------------------------------------------

cov_dir  <- file.path(paths$large, "Data_processed/cover_combined")
rap_dir  <- file.path(paths$large, "Data_processed/CoverData/rap")
soil_dir <- file.path(paths$large, "Data_processed/soils")
cell_years_vc <- "c04-c06"  # 07 and the zero-tree extraction cover all three
out_file <- file.path(cov_dir, paste0("cover_clim_soils_", vc, ".csv"))

# read --------------------------------------------------------------------

cover <- read_csv(
  file.path(cov_dir, paste0("cover_by_pixel_year_", vc, ".csv")),
  show_col_types = FALSE)

clim <- read_csv(
  file.path(cov_dir, paste0("daymet_climate-at-points_", cell_years_vc,
                            ".csv")),
  show_col_types = FALSE)

zero_tree <- read_csv(
  file.path(rap_dir, paste0("RAP_v3_zeroTree_lt3_cell-years_", cell_years_vc,
                            ".csv")),
  # set explicitly: masked cells are written as empty values, and readr
  # could guess the wrong type if the first rows are all empty
  col_types = cols(cell = "d", year = "d", fracZeroTree = "d",
                   zeroTree = "d")) |>
  transmute(cell, year,
            # % (not 0-1), like cover, so cover_threshold and the
            # quasibinomial fit (response / 100) work on it unchanged
            pct_zero_tree = 100 * fracZeroTree,
            zero_tree = as.integer(zeroTree))

mask_r <- read_mask()
soil_r <- terra::rast(file.path(soil_dir, "soil_covariates_solus100_1000m.tif"))

# checks ------------------------------------------------------------------

stopifnot(
  !anyDuplicated(cover[c("cell", "year")]),
  !anyDuplicated(clim[c("cell", "year")]),
  !anyDuplicated(zero_tree[c("cell", "year")]),
  # x, y lie in the cell they're labelled with
  all(terra::cellFromXY(mask_r, as.matrix(cover[c("x", "y")])) == cover$cell),
  # 07 and the zero-tree extraction are built from the cell-years of all three
  # 06 tables, so every cover row should be in both; a missing one means one
  # is stale (rerun it)
  nrow(anti_join(cover, clim, by = c("cell", "year"))) == 0,
  nrow(anti_join(cover, zero_tree, by = c("cell", "year"))) == 0,
  all(zero_tree$zero_tree %in% c(0, 1, NA)),
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
  left_join(zero_tree, by = c("cell", "year")) |>
  bind_cols(as_tibble(soil_r[out$cell]))

message("rows missing zero-tree data: ", sum(is.na(out$zero_tree)))

# the extraction uses the 2021 layer's definition, so its 2021 rows should
# match that layer (small differences can come from resampling)
zero_2021 <- read_zero_tree_raster("zeroTree") |>
  align_raster(mask_r)
is_2021 <- out$year == 2021
match_2021 <- out$zero_tree[is_2021] == zero_2021[out$cell[is_2021]][[1]]
message("2021 rows where zero_tree matches the 2021 layer: ",
        round(100 * mean(match_2021, na.rm = TRUE), 2), "% of ",
        sum(!is.na(match_2021)))

# needleleaf and broadleaf cover from RAP tree cover ----------------------

rap_tree <- terra::rast(file.path(rap_dir,
                                  "RAP_v3_cover-tree_2000-2024_1000m.tif")) |>
  # GEE exports come back one cell larger in each direction than the snap grid
  align_raster(mask_r)

# same-year RAP tree cover, 0 where zero_tree is 1 (left as is where zero_tree
# is NA)
rap_tree_cell <- rap_at(rap_tree, out)
rap_tree_cell[!is.na(out$zero_tree) & out$zero_tree == 1] <- 0

# field rows where RAP has no trees. Field rows only: RAP rows have no
# needle/broad fractions, so without this they'd get c_needle = 0 at zero-tree
# cells and NA everywhere else
field_no_tree <- out$sources != "rap" & !is.na(rap_tree_cell) &
  rap_tree_cell == 0

out <- out |>
  mutate(c_needle_field = c_needle,
         c_broad_field = c_broad,
         # 0 where RAP has no trees, even if the field fraction is NA (no
         # trees in the plot); this adds observations only where RAP is 0, so
         # the sample isn't random, but these columns are only for validation
         c_needle = if_else(field_no_tree, 0, frac_needle * rap_tree_cell),
         c_broad  = if_else(field_no_tree, 0, frac_broad * rap_tree_cell))

message("rows with frac_needle: ", sum(!is.na(out$frac_needle)),
        " | with c_needle (RAP tree cover available): ",
        sum(!is.na(out$c_needle)))

write_csv(out, out_file)
