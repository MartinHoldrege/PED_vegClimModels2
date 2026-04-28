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

# files created in 00_soils_aggregate.R
clay   <- rast(file.path(soil_dir, "claytotal_SOLUS100_1000m.tif"))
sand   <- rast(file.path(soil_dir, "sandtotal_SOLUS100_1000m.tif"))
coarse <- rast(file.path(soil_dir, "fragvol_SOLUS100_1000m.tif"))
resdept <- rast(file.path(soil_dir, "resdept_SOLUS100_1000m.tif"))
# stack: 7 clay bands, 7 sand bands, 7 coarse bands, 1 depth to restriction band
soil_stack <- c(clay, sand, coarse, resdept)

# layer depths ------------

solus_depths <- stringr::str_extract(names(clay), '(?<=_)\\d+(?=cm$)') |> 
  as.numeric() |> 
  sort()

stopifnot(solus_depths == c(0, 5, 15, 30, 60, 100, 150)) # check for now, would work with other depths

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

# project AWC to the Daymet grid
awc_rast <- app(soil_stack, fun = compute_awc)
names(awc_rast) <- "awc_cm"
writeRaster(awc_rast, file.path(paths$large, 
                               "./Data_processed/soils/", 
                               "awc_SOLUS100_1000m.tif"), 
                               overwrite = TRUE)


# testing
if(FALSE) {
  # layer order in x:
  # [1:7]   sand at SOLUS depths (0, 5, 15, 30, 60, 100, 150)
  # [8:14]  clay
  # [15:21] fragvol
  # [22]    resdept
  
  # Case 1: deep loamy soil, no bedrock issue
  x_deep <- c(
    # sand %
    40, 42, 45, 48, 50, 52, 53,
    # clay %
    20, 22, 23, 24, 25, 26, 27,
    # fragvol %
    5, 5, 7, 8, 10, 12, 15,
    # resdept (cm)
    201
  )
  
  # Case 2: moderately shallow soil, bedrock at 35 cm (mid-layer)
  x_shallow <- x_deep
  x_shallow[22] <- 35
  
  # Case 3: very shallow soil, bedrock at 15 cm
  x_veryshallow <- x_deep
  x_veryshallow[22] <- 15
  
  # Case 4: sandy soil, deep
  x_sandy <- c(
    85, 86, 87, 88, 89, 90, 90,
    5, 5, 5, 6, 6, 7, 7,
    2, 3, 5, 8, 10, 12, 15,
    201
  )
  
  # run them
  compute_awc(x_deep)        
  compute_awc(x_shallow)     
  compute_awc(x_veryshallow) 
  compute_awc(x_sandy)       # (sandy = lower AWC)
}
