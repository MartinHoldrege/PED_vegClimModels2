# calculated available water content 
# adapted from Alice Stear's, for calculating across a Raster

# started 4/27/26, by Martin Holdrege

library(terra)
library(rSOILWAT2)
source('Functions/init.R')
source_functions()

# ---- inputs ----
soil_dir <- file.path(paths$large, 'Data_processed/soils')

# ---- load only the bands needed ----

# files created in 01_soils_aggregate.R; order here defines band order below
properties <- c("sandtotal", "claytotal", "fragvol", "soc")

soil_stack <- c(
  rast(file.path(soil_dir, paste0(properties, "_solus100_1000m.tif"))),
  rast(file.path(soil_dir, "resdept_solus100_1000m.tif"))
)
# layer depths ------------

solus_depths <- names(soil_stack) |> 
  stringr::str_subset('sandtotal') |> 
  stringr::str_extract('(?<=_)\\d+(?=cm$)') |> 
  as.numeric() |> 
  sort()

stopifnot(solus_depths == c(0, 5, 15, 30, 60, 100, 150),
          nlyr(soil_stack) == length(properties) * 7 + 1) # check for now, would work with other depths

# what depths want to interpolate to (as in Alices/Daniels workflow)
layer_breaks <- c(0, 3, 10, 20, 30, 40, 60, 80, 100, 150, 201)  

# ---- per-pixel AWC function ---
compute_awc <- function(x) {
  sand_pts    <- x[1:7]/100 # percent to proportion
  clay_pts    <- x[8:14]/100
  fragvol_pts <- x[15:21]/100
  resdept     <- x[22]
  
  if (is.na(sand_pts[1]) || is.na(clay_pts[1]) || is.na(resdept)) {
    return(NA_real_)
  }
  
  # trapezoidal: 7-point SOLUS -> 10-layer Alice scheme
  # extrapolation = "closest" handles the 150 -> 201 extension
  sand_lyr    <- trapezoidSoilLayers(solus_depths, sand_pts,    layer_breaks, "closest")
  clay_lyr    <- trapezoidSoilLayers(solus_depths, clay_pts,    layer_breaks, "closest")
  fragvol_lyr <- trapezoidSoilLayers(solus_depths, fragvol_pts, layer_breaks, "closest")
  
  # truncate layer thicknesses at resdept
  layer_top    <- layer_breaks[-length(layer_breaks)]
  layer_bottom <- layer_breaks[-1]
  thickness    <- pmin(layer_bottom, resdept) - layer_top
  thickness[thickness <= 0] <- NA
  
  keep <- !is.na(thickness) &
    !is.na(sand_lyr) & !is.na(clay_lyr) & !is.na(fragvol_lyr)
  if (!any(keep)) return(NA_real_)
  
  sand_f    <- sand_lyr[keep]   
  clay_f    <- clay_lyr[keep]   
  fragvol_f <- fragvol_lyr[keep] 
  
  p <- rSOILWAT2::ptf_estimate(
    sand = sand_f, clay = clay_f, fcoarse = fragvol_f,
    swrc_name = "Campbell1974", ptf_name = "Cosby1984"
  )
  vwc <- rSOILWAT2::swrc_swp_to_vwc(
    c(-1.5, -0.033),
    fcoarse = fragvol_f,
    swrc = list(name = "Campbell1974", swrcp = p)
  )
  
  sum(thickness[keep] * as.vector(diff(vwc)))
}

# ---- per-pixel soil covariates ------------------------------------------

#' Calculate soil covariates for one pixel
#'
#' @param x Numeric vector of stacked SOLUS values: 7 depths each of sand,
#'   clay, coarse fragments and organic carbon (percent), then depth to
#'   restriction (cm).
#' @param solus_depths Point depths (cm) of the SOLUS predictions.
#' @param layer_breaks Target layer boundaries (cm).
#' @return Named numeric vector of length 7. Percentages stay as percentages;
#'   soilDepth and AWHC are in cm.
compute_soil_vars <- function(x, solus_depths, layer_breaks) {
  
  # every return path must have length 7, or app() errors with
  # "'fun' returns a list"
  na_out <- c(soilDepth = NA_real_, clay_surface = NA_real_, clay = NA_real_,
              sand = NA_real_, coarse = NA_real_, carbon = NA_real_,
              AWHC = NA_real_)
  
  sand_pts   <- x[1:7]
  clay_pts   <- x[8:14]
  coarse_pts <- x[15:21]
  carbon_pts <- x[22:28]
  resdept    <- x[29]
  
  if (is.na(sand_pts[1]) || is.na(clay_pts[1]) || is.na(resdept)) {
    return(na_out)
  }
  
  # trapezoidal: 7-point SOLUS -> layer averages
  # extrapolation = "closest" handles the 150 -> 201 extension
  to_layers <- function(v) {
    trapezoidSoilLayers(solus_depths, v, layer_breaks, "closest")
  }
  sand_lyr   <- to_layers(sand_pts)
  clay_lyr   <- to_layers(clay_pts)
  coarse_lyr <- to_layers(coarse_pts)
  carbon_lyr <- if (all(is.na(carbon_pts))) rep(NA_real_, length(layer_breaks) - 1) else to_layers(carbon_pts)
  
  # truncate layer thicknesses at resdept
  layer_top    <- layer_breaks[-length(layer_breaks)]
  layer_bottom <- layer_breaks[-1]
  thickness    <- pmin(layer_bottom, resdept) - layer_top
  thickness[thickness <= 0] <- NA
  
  keep <- !is.na(thickness) &
    !is.na(sand_lyr) & !is.na(clay_lyr) & !is.na(coarse_lyr)
  if (!any(keep)) return(na_out)
  
  thick      <- thickness[keep]
  profile_thickness <- sum(thick)
  # profile_thickness and soilDepth, generally should be the same
  soilDepth = pmin(resdept, max(layer_breaks)) 
  # thickness-weighted profile means
  wmean <- function(v) sum(v[keep] * thick) / profile_thickness
  
  # AWHC; the PTF takes proportions, not percentages
  fcoarse <- coarse_lyr[keep] / 100
  p <- rSOILWAT2::ptf_estimate(
    sand = sand_lyr[keep] / 100,
    clay = clay_lyr[keep] / 100,
    fcoarse = fcoarse,
    swrc_name = "Campbell1974", ptf_name = "Cosby1984"
  )
  vwc <- rSOILWAT2::swrc_swp_to_vwc(
    c(-1.5, -0.033),
    fcoarse = fcoarse,
    swrc = list(name = "Campbell1974", swrcp = p)
  )
  
  c(soilDepth    = soilDepth,
    clay_surface = unname(clay_lyr[1]),    # 0-3 cm
    clay         = wmean(clay_lyr),
    sand         = wmean(sand_lyr),
    coarse       = wmean(coarse_lyr),
    carbon       = unname(carbon_lyr[1]),  # 0-3 cm
    AWHC         = sum(thick * as.vector(diff(vwc))))
}


# ---- calculate and write ------------------------------------------------

soil_vars <- app(soil_stack, fun = compute_soil_vars,
                 solus_depths = solus_depths, layer_breaks = layer_breaks)
names(soil_vars) <- c("soilDepth", "clay_surface", "clay", "sand", "coarse",
                      "carbon", "AWHC")

writeRaster(soil_vars,
            file.path(soil_dir, "soil_covariates_solus100_1000m.tif"),
            overwrite = TRUE)

# testing -----------------------------------------------------------------

if (FALSE) {
  # layer order in x:
  # [1:7]   sand at SOLUS depths (0, 5, 15, 30, 60, 100, 150)
  # [8:14]  clay
  # [15:21] coarse fragments
  # [22:28] organic carbon
  # [29]    resdept
  
  # deep loamy soil, no bedrock restriction
  x_deep <- c(
    40, 42, 45, 48, 50, 52, 53,      # sand %
    20, 22, 23, 24, 25, 26, 27,      # clay %
    5, 5, 7, 8, 10, 12, 15,          # coarse %
    3, 2.5, 1.5, 1, 0.6, 0.4, 0.3,   # carbon %
    201                              # resdept (cm)
  )
  
  x_shallow     <- replace(x_deep, 29, 35)  # bedrock mid-layer
  x_veryshallow <- replace(x_deep, 29, 15)
  
  # sandy soil, deep (expect lower AWHC)
  x_sandy <- c(
    85, 86, 87, 88, 89, 90, 90,
    5, 5, 5, 6, 6, 7, 7,
    2, 3, 5, 8, 10, 12, 15,
    1, 0.8, 0.5, 0.3, 0.2, 0.2, 0.1,
    201
  )
  
  compute_soil_vars(x_deep, solus_depths, layer_breaks)
  compute_soil_vars(x_shallow, solus_depths, layer_breaks)
  compute_soil_vars(x_veryshallow, solus_depths, layer_breaks)
  compute_soil_vars(x_sandy, solus_depths, layer_breaks)
}
