# 06_build_cover_datasets.R
#
# Builds the cover tables for three training-data versions:
#   c04  field only: the field pixel-years from 05, no RAP
#   c05  RAP only: a random sample of RAP pixel-years
#   c06  blended, for shrub and herbaceous: the c05 rows with same-year RAP
#        tree cover < tree_cutoff, plus field rows where it is
#        >= tree_cutoff (where RAP understorey isn't usable). The same tree
#        cover value decides both, so a cell-year can't be in both
#
# RAP sample (c05):
#   - pool: snap cells with an ecoregion, kept by the LCMAP mask, and with
#     at least min_years usable years in rap_years. A year is usable if the
#     cell passes that year's fire mask (< 10% burned in the preceding 20
#     years) and has RAP cover. Cells with field data can be drawn too.
#   - n_cells cells drawn at random from the pool, then n_years_cell years
#     drawn at random from each cell's usable years
#   - shrub and herbaceous are NA where same-year RAP tree cover is
#     >= tree_cutoff; bare ground is kept at all tree cover
#
# Field rows (c04) use the same site-year filters as before: the per-year
# fire mask for all, plus the LCMAP mask except for FIA-only pixel-years (FIA
# data was already filtered to exclude developed and agricultural conditions
# at the plot site).
#
# All three tables have the same columns; variables a version doesn't have
# are NA. Also writes the combined list of cell-years across the three
# versions, so climate at points (07) only has to be computed once.
#
# Inputs:
#   field_cover_by_pixel_year.csv        - 05_combine_field_cover.R
#   RAP_v3_cover-<group>_2000-2024_1000m.tif
#                                        - 01_rap_cover.js, via
#                                          04_download_gee_output.R
#   LCMAP_fracKeep_gte90_1000m.tif       - 03_export_masks.js, via 04
#   MTBS_fracUnburned_gte90_20yr_2000-2024_1000m.tif
#                                        - 03_export_masks.js, via 04
#   EPA_L3_ecoregion_daymet_1000m.tif, EPA_L3_ecoregion_lookup.csv
#                                        - 01_rasterize_ecoregions.R
#   daymet_conus_snap_1000m.tif          - 00_create_snap_raster.R
#                                          (read_mask())
#
# Outputs:
#   cover_by_pixel_year_<vc>.csv, for c04, c05 and c06
#   cover_cell_years_c04-c06.csv
#
# Replaces 06_cover_add-rap.R.
#
# October 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

vc_field <- "c04"
vc_rap   <- "c05"
vc_blend <- "c06"

n_cells      <- 1e5
n_years_cell <- 3
# 2010 is the first year with a full 30-yr normal. Field rows keep all years
# (2000 on, >= 20-yr normals) to retain as much field data as possible
rap_years    <- 2010:2024
min_years    <- 7          # usable years a cell needs to enter the pool
tree_cutoff  <- 10         # percent; RAP understorey not usable at or above
seed         <- 1234

cov_dir  <- file.path(paths$large, "Data_processed/cover_combined")
rap_dir  <- file.path(paths$large, "Data_processed/CoverData/rap")
mask_dir <- file.path(paths$large, "Data_processed/masks")

groups <- c("tree", "shrub", "herbaceous", "bare_ground")
rap_files <- file.path(rap_dir, paste0("RAP_v3_cover-", groups,
                                       "_2000-2024_1000m.tif"))
names(rap_files) <- groups

# columns only field data carry; NA in RAP rows
field_only_cols <- c("c_forb", "c_c3", "c_c4", "c_gram", "c_needle",
                     "c_broad", "frac_forb", "frac_c3", "frac_c4",
                     "frac_needle", "frac_broad")


# functions ---------------------------------------------------------------

#' RAP cover at given cell-years
#'
#' @param r SpatRaster with one band per year, named 'year_YYYY'.
#' @param cell_years Data frame with `cell` and `year` columns.
#' @return Numeric vector of cover, one value per row of `cell_years`.
#' @examples
#' r <- rap$tree
#' cell_years <- tibble(cell = c(4190, 4191), year = c(2015L, 2020L))
#' rap_at(r = r, cell_years = cell_years)
rap_at <- function(r, cell_years) {
  # one row per cell-year, one column per year
  vals <- as.matrix(r[cell_years$cell])
  col <- match(paste0("year_", cell_years$year), names(r))
  stopifnot(!anyNA(col))
  vals[cbind(seq_len(nrow(cell_years)), col)]
}

#' Final formatting shared by the three cover tables
#'
#' Adds NA columns for any field-only variables the table lacks (so all
#' versions have the same columns), adds ecoregion names and sorts.
#'
#' @param d Cover table with a `region` column (ecoregion id).
#' @param eco_lookup Ecoregion lookup, from load_ecoregion_lookup().
#' @param field_only_cols Columns that only field data carry.
#' @examples
#' d <- rap_rows
#' eco_lookup <- load_ecoregion_lookup("L3")
#' field_only_cols <- c("c_forb", "frac_forb")
#' format_cover(d = d, eco_lookup = eco_lookup,
#'              field_only_cols = field_only_cols)
format_cover <- function(d, eco_lookup, field_only_cols) {
  for (cc in setdiff(field_only_cols, names(d))) {
    d[[cc]] <- NA_real_
  }
  d |>
    left_join(eco_lookup, by = c("region" = "eco_id")) |>
    select(-region) |>
    arrange(eco_code, cell, year)
}


# read --------------------------------------------------------------------

snap <- read_mask()

field <- read_csv(file.path(cov_dir, "field_cover_by_pixel_year.csv"),
                  show_col_types = FALSE)

eco        <- load_ecoregion_raster("L3")
eco_lookup <- load_ecoregion_lookup("L3")

# GEE exports come back one cell larger in each direction than the snap grid
lcmap_mask <- terra::rast(file.path(mask_dir,
                                    "LCMAP_fracKeep_gte90_1000m.tif")) |>
  align_raster(snap)
fire_mask <- terra::rast(file.path(
  mask_dir, "MTBS_fracUnburned_gte90_20yr_2000-2024_1000m.tif")) |>
  align_raster(snap)
rap <- map(rap_files, \(f) align_raster(terra::rast(f), snap))

year_bands <- paste0("year_", rap_years)

stopifnot(
  all(terra::cellFromXY(snap, as.matrix(field[c("x", "y")])) == field$cell),
  terra::compareGeom(eco, snap),
  all(map_lgl(rap, \(r) identical(names(r), paste0("year_", 2000:2024)))),
  all(year_bands %in% names(fire_mask))
)


# field rows (c04) --------------------------------------------------------

field$region <- eco[field$cell][[1]]

# checking for NA's in masks that shouldn't be present
stopifnot(all(!is.na(field$masked_by_fire)),
          all(!is.na(field$masked_by_lcmap)))

field_final <- field |>
  filter(!is.na(region)) |> # this filter not actually removing any rows
  # FIA-only pixel-years are masked by fire only; everything else by both.
  filter(if_else(sources == "fia",
                 !masked_by_fire,
                 !masked_by_fire & !masked_by_lcmap)) |>
  select(-masked_by_fire, -masked_by_lcmap)

# RAP sample (c05) --------------------------------------------------------

# usable cell-years: pass that year's fire mask and have RAP cover
usable <- fire_mask[[year_bands]] == 1 & !is.na(rap$tree[[year_bands]])
names(usable) <- year_bands
n_usable <- sum(usable, na.rm = TRUE)

pool_ok <- !is.na(snap) & !is.na(eco) & lcmap_mask == 1 &
  n_usable >= min_years
pool <- which(terra::values(pool_ok, mat = FALSE) == 1)

message("cells in the RAP pool: ", length(pool))
stopifnot(length(pool) >= n_cells)

set.seed(seed)
rap_cells <- sort(sample(pool, n_cells))

rap_cell_years <- as_tibble(usable[rap_cells]) |>
  mutate(cell = rap_cells) |>
  pivot_longer(-cell, names_to = "year", values_to = "is_usable") |>
  filter(is_usable == 1) |>
  mutate(year = as.integer(str_remove(year, "^year_"))) |>
  slice_sample(n = n_years_cell, by = cell) |>
  select(cell, year)

stopifnot(all(count(rap_cell_years, cell)$n == n_years_cell))

xy <- terra::xyFromCell(snap, rap_cell_years$cell)

rap_rows <- rap_cell_years |>
  mutate(cov_tree        = rap_at(rap$tree,        pick(cell, year)),
         cov_shrub       = rap_at(rap$shrub,       pick(cell, year)),
         # AFG and PFG are predicted separately, so their sum can exceed 100
         cov_herbaceous  = pmin(rap_at(rap$herbaceous, pick(cell, year)), 100),
         cov_bare_ground = rap_at(rap$bare_ground, pick(cell, year)),
         x = xy[, "x"], y = xy[, "y"],
         region = eco[cell][[1]],
         sources = "rap", n_plots = NA_integer_)

stopifnot(!anyNA(select(rap_rows, -n_plots)))

message("RAP pixel-years with tree cover >= ", tree_cutoff, "%: ",
        round(100 * mean(rap_rows$cov_tree >= tree_cutoff), 1), "%")

# RAP understorey is not usable under tree canopy
rap_rows <- rap_rows |>
  mutate(cov_shrub      = if_else(cov_tree >= tree_cutoff, NA_real_,
                                  cov_shrub),
         cov_herbaceous = if_else(cov_tree >= tree_cutoff, NA_real_,
                                  cov_herbaceous))


# field shrub and herbaceous under tree canopy (c06) ----------------------

field_final$rap_tree <- rap_at(rap$tree, field_final)

# pixel-years with no RAP tree cover can't be classed, so drop out here
field_forest <- field_final |>
  filter(rap_tree >= tree_cutoff) |>
  filter(!is.na(cov_shrub) | !is.na(cov_herbaceous)) |>
  select(cell, year, x, y, region, sources, n_plots,
         cov_shrub, cov_herbaceous)

message("field pixel-years added to c06, by source:")
count(field_forest, sources, sort = TRUE) |> as.data.frame() |> print()


# * combine -----------------------------------------------------------------

# c06: RAP below the tree cutoff, field at or above it
blend <- bind_rows(filter(rap_rows, cov_tree < tree_cutoff), field_forest)

cover <- list(
  field = format_cover(select(field_final, -rap_tree), eco_lookup,
                       field_only_cols),
  rap   = format_cover(rap_rows, eco_lookup, field_only_cols),
  blend = format_cover(blend, eco_lookup, field_only_cols)
)
names(cover) <- c(vc_field, vc_rap, vc_blend)

# combined cell-years for 07 and 07_rap_zero-tree_cell-years.js; c06's are a
# subset of c04's and c05's. longitude, latitude (cell centres) are for the
# GEE table upload
cell_years <- bind_rows(cover) |>
  distinct(cell, year, x, y) |>
  arrange(cell, year)

lonlat <- terra::vect(as.matrix(cell_years[c("x", "y")]),
                      crs = terra::crs(snap)) |>
  terra::project("EPSG:4326") |>
  terra::crds()

cell_years$longitude <- lonlat[, "x"]
cell_years$latitude  <- lonlat[, "y"]

# checks ------------------------------------------------------------------

walk(cover, \(d) {
  stopifnot(
    !anyDuplicated(d[c("cell", "year")]),
    !anyNA(d[c("cell", "year", "x", "y", "sources", "eco_code", "eco_name")]),
    setequal(names(d), names(cover[[vc_field]])),
    all(map_lgl(d[str_subset(names(d), "^cov_")],
                \(z) all(is.na(z) | between(z, 0, 100))))
  )
})

stopifnot(!anyDuplicated(cell_years[c("cell", "year")]))

message("pixel-years by version:")
map_int(cover, nrow) |> print()


# write -------------------------------------------------------------------

iwalk(cover, \(d, vc) {
  write_csv(d, file.path(cov_dir,
                         paste0("cover_by_pixel_year_", vc, ".csv")))
})

write_csv(cell_years, file.path(cov_dir, "cover_cell_years_c04-c06.csv"))
