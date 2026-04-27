# calculated available water content 
# adapted from Alice Stear's, for calculating across a Raster

# started 4/27/26, by Martin Holdrege

library(terra)
library(rSOILWAT2)
source('Functions/init.R')
source_functions()

# ---- inputs ----
soil_dir <- "E:/USGS/large_files/PED_vegClimModels/Data_raw/soilsDB_new"

# point depths (cm) — from Daniel's preprocessing
depths <- c(2, 7, 15, 25, 35, 50, 70, 90, 125, 176)



# ---- load only the bands needed ----

clay   <- rast(file.path(soil_dir, "claytotal_PED-CONUS4km_SOLUS100.nc"))
sand   <- rast(file.path(soil_dir, "sandtotal_PED-CONUS4km_SOLUS100.nc"))
coarse <- rast(file.path(soil_dir, "fragvol_PED-CONUS4km_SOLUS100.nc"))

# stack: 10 clay bands, 10 sand bands, 10 coarse bands (30 total, in that order)
soil_stack <- c(clay, sand, coarse)

# ---- per-pixel AWC function ----

#' Compute total profile available water-holding capacity (cm) for one pixel
#'
#' @param x Numeric vector of length 30: clay (10), sand (10), coarse (10),
#'   each in percent, ordered shallow to deep.
#' @return Total profile AWC in cm, or NA if surface layer is missing.
compute_awc <- function(x) {
  clay   <- x[1:10]  / 100
  sand   <- x[11:20] / 100
  coarse <- x[21:30] / 100
  
  # layer thicknesses (cm) corresponding to each point depth
  # layer boundaries: 0-3, 3-10, 10-20, 20-30, 30-40, 40-60, 60-80, 80-100, 100-150, 150-201
  thickness <- c(3, 7, 10, 10, 10, 20, 20, 20, 50, 51)
  
  
  # require surface layer; drop deeper NA layers (e.g., bedrock truncation)
  if (is.na(sand[1]) || is.na(clay[1])) return(NA_real_)
  keep <- !is.na(sand) & !is.na(clay) & !is.na(coarse)
  if (sum(keep) == 0) return(NA_real_)
  
  p <- rSOILWAT2::ptf_estimate(
    sand = sand[keep], clay = clay[keep], fcoarse = coarse[keep],
    swrc_name = "Campbell1974", ptf_name = "Cosby1984"
  )
  vwc <- rSOILWAT2::swrc_swp_to_vwc(
    c(-1.5, -0.033),  # FC, WP — order chosen so diff is negative, abs gives AWC
    fcoarse = coarse[keep],
    swrc = list(name = "Campbell1974", swrcp = p)
  )
  awc_per_layer <- abs(as.vector(diff(vwc)))  # cm³/cm³
  sum(thickness[keep] * awc_per_layer)        # cm of water
}

# ---- apply to raster ----
awc_rast <- app(soil_stack, fun = compute_awc, cores = 4)
names(awc_rast) <- "awc_cm"


# project AWC to the Daymet grid
awc_1km <- project(awc_rast, daymet_grid, method = "bilinear")

writeRaster(awc_1km, file.path(paths$large, 
                               "./Data_processed/soils/", 
                               "awc_SOLUS100_1000m.tif"), 
                               overwrite = TRUE)
