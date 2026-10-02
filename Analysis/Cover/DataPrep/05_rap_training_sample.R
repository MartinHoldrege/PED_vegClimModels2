# 05_rap_training_sample.R
#
# Training data from RAP alone: a simple random sample of cells from the
# wall-to-wall 2021 RAP cover (LCMAP and fire masked in GEE), with the
# 1991-2020 climate normals and soils at each cell. For trying the forest
# classification on a large, spatially even dataset. No anomalies.
#
# Written in the same format as 08_combine_cover_and_covariates.R (long
# _CLIM climate names, NA columns for the field-only variables), except soil
# AWC is written as awc, so it can be read with read_cover_training("c02").
#
# Inputs:
#   RAP_v3_cover_2021_1000m.tif          - 03_rap_sample.js (04_download)
#   DaymetClimateData_1991-2020_CLIM.tif - 01_summarise_daymet_climate_data.R
#   soil_covariates_solus100_1000m.tif   - 02_soils_calculate_variables.R
#   EPA_L3_ecoregion_daymet_1000m.tif    - 01_rasterize_ecoregions.R
#   daymet_conus_snap_1000m.tif          - 00_create_snap_raster.R
#
# Output:
#   cover_clim_soils_<vc>.csv (read with read_cover_training())
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

vc          <- "c02"  # cover data version this script creates
year        <- 2021
n_sample    <- 3e5
tree_cutoff <- 10     # percent; RAP understorey invalid at or above (as in 06)
seed        <- 5817

rap_file  <- file.path(paths$large, "Data_processed/CoverData/rap",
                       paste0("RAP_v3_cover_", year, "_1000m.tif"))
clim_file <- file.path(paths$large, "Data_processed/WallToWallClimateData",
                       "DaymetClimateData_1991-2020_CLIM.tif")
out_file  <- file.path(paths$large, "Data_processed/cover_combined",
                       paste0("cover_clim_soils_", vc, ".csv"))

# read --------------------------------------------------------------------

snap <- read_mask()

rap <- terra::rast(rap_file) |>
  # GEE exports come back one cell larger in each direction than the snap grid
  align_raster(snap)
stopifnot(names(rap) == c("tree", "shrub", "herbaceous", "bare_ground"))
names(rap) <- paste0("cov_", names(rap))

# climate (short names) and soils
clim_r <- read_climate_raster(path = clim_file)
# long names, as written by 08; read_cover_training() shortens them again
clim_long <- names(terra::rast(clim_file))
clim_short <- climate_name_lookup(clim_long)
soil_vars <- setdiff(names(clim_r), clim_short)  # includes awc

eco        <- load_ecoregion_raster("L3")
eco_lookup <- load_ecoregion_lookup("L3")

stopifnot(
  terra::compareGeom(clim_r, snap),
  terra::compareGeom(eco, snap),
  all(clim_short %in% names(clim_r)),
  "awc" %in% soil_vars
)

# sample ------------------------------------------------------------------

# cells with RAP tree cover (i.e. passing the masks), climate and an ecoregion;
# rows missing soils are kept, as in 08
ok <- !is.na(snap) & !is.na(rap[["cov_tree"]]) & !is.na(clim_r[[1]]) & 
  !is.na(clim_r$awc) & !is.na(eco)

pool <- which(terra::values(ok, mat = FALSE) == 1)

message("cells available to sample: ", length(pool))
stopifnot(length(pool) >= n_sample)

set.seed(seed)
cell <- sort(sample(pool, n_sample))
xy <- terra::xyFromCell(snap, cell)

vals <- c(rap, clim_r, eco)[cell] |>
  as_tibble()

# combine -----------------------------------------------------------------

# field-only columns, NA for RAP rows (as in 08's output)
na_cols <- c("c_forb", "c_c3", "c_c4", "c_gram", "c_needle", "c_broad",
             "frac_forb", "frac_c3", "frac_c4", "frac_needle", "frac_broad")

out <- vals |>
  mutate(cell = cell, year = year, x = xy[, "x"], y = xy[, "y"],
         sources = "rap", n_plots = NA_integer_) |>
  left_join(eco_lookup, by = c("region" = "eco_id")) |>
  # RAP understorey is not usable under tree canopy (as in 06)
  mutate(under_tree     = cov_tree >= tree_cutoff,
         cov_shrub      = if_else(under_tree, NA_real_, cov_shrub),
         cov_herbaceous = if_else(under_tree, NA_real_, cov_herbaceous)) |>
  rename(any_of(setNames(clim_short, clim_long))) |>
  select(cell, year, cov_tree, cov_shrub, cov_herbaceous, cov_bare_ground,
         n_plots, sources, x, y, eco_code, eco_name, all_of(clim_long), 
         all_of(soil_vars))

# checks ------------------------------------------------------------------

stopifnot(
  nrow(out) == n_sample,
  !anyDuplicated(out$cell),
  all(terra::cellFromXY(snap, as.matrix(out[c("x", "y")])) == out$cell),
  !anyNA(out[c("cell", "year", "x", "y", "eco_code", "eco_name", "cov_tree")])
)

message("rows with tree cover > ", tree_cutoff, "%: ",
        round(100 * mean(out$cov_tree > tree_cutoff), 1), "%")
message("rows missing soils: ", sum(!complete.cases(out[soil_vars])))

write_csv(out, out_file)
