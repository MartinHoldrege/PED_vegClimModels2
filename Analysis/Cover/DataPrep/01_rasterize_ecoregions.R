# Purpose:
# Rasterize EPA Level 2 or Level 3 ecoregions to the daymet 1km grid.
# Each cell gets an integer ecoregion ID.
# Output: GeoTIFF with ecoregion IDs on the daymet grid.
#
# Author: Martin Holdrege


# params ------------------------------------------------------------------

# rasterize either the L2 or L3 ecoregions
epa_lev <- 'L3' # 'L2' # 

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# paths --------------------------------------------------------------------


p_shp <- if(epa_lev == 'L3')  {
  file.path( paths$large, "Data_raw/Level3Ecoregions/us_eco_l3.shp")
} else if(epa_lev == 'L2') {
  file.path(paths$large0, "Data_raw/Level2Ecoregions/NA_CEC_eco_Level2.shp")
} else {
  stop('selected region level not available')
}

stopifnot(file.exists(p_shp))

out_dir <- file.path(paths$large, "Data_processed/regions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
p_out <- file.path(out_dir, 
                   paste0("EPA_", epa_lev, "_ecoregion_daymet_1000m.tif"))

# column name mapping by level ----------------------------------------------

if (epa_lev == "L3") {
  code_col <- "US_L3CODE"
  name_col <- "US_L3NAME"
} else if (epa_lev == "L2") {
  code_col <- "NA_L2CODE"
  name_col <- "NA_L2NAME"
}

# load ecoregion shapefile -------------------------------------------------

eco <- sf::st_read(p_shp, quiet = TRUE)

# load a daymet template raster --------------------------------------------

r_template <- read_mask()

# project ecoregions to daymet CRS ----------------------------------------

eco_proj <- sf::st_transform(eco, crs = terra::crs(r_template))

# create numeric ecoregion ID ----------------------------------------------

# code column is character; create an integer mapping
eco_proj$eco_id <- as.integer(as.factor(eco_proj[[code_col]]))

# save lookup table
eco_lookup <- eco_proj |>
  sf::st_drop_geometry() |>
  dplyr::distinct(eco_id, .data[[code_col]], .data[[name_col]]) |>
  dplyr::arrange(eco_id)

cat("  Ecoregion ID range:", range(eco_lookup$eco_id), "\n")

p_lookup <- file.path(out_dir, paste0("EPA_", epa_lev, "_ecoregion_lookup.csv"))
utils::write.csv(eco_lookup, p_lookup, row.names = FALSE)

# rasterize ----------------------------------------------------------------

eco_vect <- terra::vect(eco_proj)
r_eco <- terra::rasterize(eco_vect, r_template, field = "eco_id")
names(r_eco) <- "ecoregion"

cat("  N unique ecoregions in raster:", length(unique(terra::values(r_eco, na.rm = TRUE))), "\n")

# check for NAs where climate data exists ----------------------------------

has_climate <- !is.na(r_template)
has_eco <- !is.na(r_eco)

n_climate <- terra::global(has_climate, "sum", na.rm = TRUE)[[1]]
n_eco <- terra::global(has_eco, "sum", na.rm = TRUE)[[1]]
n_gap <- terra::global(has_climate & !has_eco, "sum", na.rm = TRUE)[[1]]

cat("  Pixels with climate data:", n_climate, "\n")
cat("  Pixels with ecoregion:   ", n_eco, "\n")
cat("  Gap (climate but no eco):", n_gap, "\n")

if (n_gap > 0) {
  pct_gap <- round(n_gap / n_climate * 100, 2)
  cat("  Gap is", pct_gap, "% of climate pixels\n")
  cat("  Consider filling gaps with nearest-neighbor or buffering polygons.\n")
} 

# fill border gaps with nearest ecoregion ----------------------------------
# Gaps are coastal and border fringe: snap cells whose centre falls just
# outside an ecoregion polygon. Filled with the modal ecoregion of the
# surrounding cells, iterating so that gaps more than one cell wide close.
# Only cells with climate data are filled, and only from existing values.
if (n_gap > 0) {
  max_iter <- 10
  
  for (i in seq_len(max_iter)) {
    gap <- has_climate & is.na(r_eco)
    n_left <- terra::global(gap, "sum", na.rm = TRUE)[[1]]
    if (n_left == 0) break
    
    filled <- terra::focal(r_eco, w = 3, fun = "modal", na.policy = "only",
                           na.rm = TRUE)
    r_eco <- terra::cover(r_eco, terra::mask(filled, has_climate,
                                             maskvalue = FALSE))
    cat("  fill pass", i, ": ", n_left, " gap cell(s) remaining\n", sep = "")
  }
  
  n_gap_final <- terra::global(has_climate & is.na(r_eco), "sum",
                               na.rm = TRUE)[[1]]
  cat("  Gap after filling:", n_gap_final, "\n")
}


# save ---------------------------------------------------------------------

terra::writeRaster(r_eco, p_out, overwrite = TRUE, datatype = "INT2S")

