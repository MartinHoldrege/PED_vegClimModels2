# purpose: aggregate 100m solus data to 1km daymet grid (specifically
# the layers needed for calculating available water content)
# files being loaded are listed here: https://storage.googleapis.com/solus100pub/index.html

# started: April 28, 2026

library(terra)
source('Functions/init.R')
source_functions()

# ---- config ----
out_dir <- file.path(paths$large, "data_processed/soils")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

properties <- c("sandtotal", "claytotal", "fragvol")
depths     <- c(0, 5, 15, 30, 60, 100, 150)

base_url <- "https://storage.googleapis.com/solus100pub/"


#' Stream + project a single SOLUS100 layer to Daymet 1 km grid
#'
#' Uses area-weighted mean ('average') so each Daymet 1 km cell is the mean
#' of the ~100 underlying SOLUS100 100 m cells that overlap it.
aggregate_solus_layer <- function(property, depth, target, unit_name = 'cm') {
  url <- paste0(base_url, property, "_", depth, "_cm_p.tif")
  vsi <- paste0("/vsicurl/", url)
  
  message("Streaming: ", property, " ", depth, " cm")
  r <- rast(vsi)
  
  out <- project(r, target, method = "mean")
  names(out) <- paste0(property, "_", depth, unit_name)
  out
}

# test <- aggregate_solus_layer(property = "sandtotal", depth = 0, target = daymet_grid)
# plot(test)

build_property_stack <- function(property, depths, target, out_dir, ...) {
  out_path <- file.path(out_dir, paste0(property, "_solus100_1000m.tif"))
  
  if (file.exists(out_path)) {
    message("Skipping (exists): ", basename(out_path))
    return(rast(out_path))
  }
  
  layers <- purrr::map(depths, function(d) aggregate_solus_layer(property, d, target, ...))
  
  stk <- rast(layers)
  writeRaster(stk, out_path, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  message("Wrote: ", basename(out_path))
  stk
}

options(timeout = 3600)

stacks <- purrr::map(
  properties,
  ~ build_property_stack(.x, depths, daymet_grid, out_dir)
) |>
  set_names(properties)

# download resdept
build_property_stack('resdept', 'all', daymet_grid, out_dir, unit_name = '')
