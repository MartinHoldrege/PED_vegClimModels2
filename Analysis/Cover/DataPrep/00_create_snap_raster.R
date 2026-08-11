# Purpose:
# Create a snap/mask raster on the daymet 1 km grid, cropped to CONUS.
# Cells are 1 where daymet has data (and fall within CONUS), NA elsewhere.
# The output defines the analysis grid (crs, resolution, origin, extent) for
# downstream rasters.
#
# Author: Martin Holdrege
# August 2026
# params ------------------------------------------------------------------

# keep cells whose centers fall outside the CONUS polygon but which the
# polygon touches (permissive coastal edge)
touches <- TRUE

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# paths -------------------------------------------------------------------

p_daymet <- file.path(paths$large0, # this can be changed to paths$large once files are moved
                      "Data_raw/dayMet/yearly/daymet_v4_prcp_annttl_na_1980.tif")

# level 3 ecoregions used as the CONUS boundary (level 2 covers all of
# North America and has no country attribute)
p_shp <- file.path(paths$large, "Data_raw/Level3Ecoregions/us_eco_l3.shp")

stopifnot(file.exists(p_daymet), file.exists(p_shp))

out_dir <- file.path(paths$large, "Data_processed", 'masks')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_out <- file.path(out_dir, "daymet_conus_snap_1000m.tif")

# load daymet template ----------------------------------------------------

r_daymet <- terra::rast(p_daymet)

# CONUS boundary ----------------------------------------------------------

conus <- sf::st_read(p_shp, quiet = TRUE) |>
  sf::st_transform(crs = terra::crs(r_daymet)) |>
  sf::st_union() |>          # dissolve ecoregions to a single outline
  sf::st_make_valid()

conus_v <- terra::vect(conus)

# crop and mask -----------------------------------------------------------

# crop first (cheap) so the mask operates on a much smaller raster;
# snap = "out" keeps the cropped extent aligned to the daymet grid
r_crop <- terra::crop(r_daymet, conus_v, snap = "out")
r_mask <- terra::mask(r_crop, conus_v, touches = touches)

# 1 where daymet has data, NA elsewhere
r_snap <- terra::ifel(is.na(r_mask), NA, 1L)

r_snap <- terra::trim(r_snap)

# cell numbers, so points can be matched to cells later.
# must come after trim(): cell numbers depend on the extent/dim
r_cell <- terra::init(terra::rast(r_snap), "cell")
r_snap <- terra::mask(r_cell, r_snap)

names(r_snap) <- "cell_id"

# checks ------------------------------------------------------------------

# grid alignment must be preserved (same res, and origin offset by whole cells)
stopifnot(
  all(terra::res(r_snap) == terra::res(r_daymet)),
  all(abs(terra::origin(r_snap) - terra::origin(r_daymet)) < 1e-6)
)

n_cell <- terra::global(!is.na(r_snap), "sum", na.rm = TRUE)[[1]]

cat("snap raster:\n")
cat("  dim:   ", dim(r_snap)[1:2], "\n")
cat("  extent:", as.vector(terra::ext(r_snap)), "\n")
cat("  n cells with data:", n_cell, "\n")


# create sparse grid ------------------------------------------------------
# this is for extracting rap data at a lower
# thinning factor: keep every 5th row and every 5th column
thin_by <- 5

r_row <- terra::init(terra::rast(r_snap), "row")
r_col <- terra::init(terra::rast(r_snap), "col")

keep <- (r_row %% thin_by == 1) & (r_col %% thin_by == 1)

r_sparse <- terra::mask(r_snap, keep, maskvalues = c(0, NA))
names(r_sparse) <- "cell_id"

p_sparse <- file.path(out_dir,
                      paste0("daymet_conus_snap_1000m_thin", thin_by, ".tif"))

if(FALSE) {
  # zoom in to a section of the sparse raster to see what it actually looks like
  e <- terra::ext(r_sparse)
  xc <- mean(e[1:2]); yc <- mean(e[3:4])
  e_zoom <- terra::ext(xc, xc + 100000, yc, yc + 100000)
  terra::plot(terra::crop(r_sparse, e_zoom))
}

# save --------------------------------------------------------------------

terra::writeRaster(r_sparse, p_sparse, overwrite = TRUE, datatype = "INT4S")

terra::writeRaster(r_snap, p_out, overwrite = TRUE, datatype = "INT4S")
