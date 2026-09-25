# 03_predict_raster.R
#
# Wall-to-wall CONUS predictions from one fitted cover model, chosen by
# --cover_type, --cover_model, --vc and --vmc (params.R). Uses
# predict_raster(), which dispatches on the model's class, so the same script
# serves every model type that has a predict_raster method.
#
# Predicts on every CONUS cell (snap mask); masking to the available area is
# left to the maps.
#
# Inputs:
#   fitted model - 02_fit_classification.R (read_cover_model())
#   DaymetClimateData_1991-2020_CLIM.tif, soils - read_climate_raster()
#
# Output (one file per climate scenario: current, BNU-ESM, IPSL-CM5A-MR):
#   Data_processed/CoverData/Predictions/
#     <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# read --------------------------------------------------------------------

fit <- read_cover_model(opt$cover_type, opt$cover_model, opt$vc, opt$vmc)
if(is.null(fit$config$spec$engine)) fit$config$spec$engine <- 'glmnet' # for legacy reasons
mask_r <- read_mask()

# climate scenarios to predict for (see read_climate_raster())
scenarios <- c("current", "BNU-ESM", "IPSL-CM5A-MR")

out_dir <- file.path(paths$large, "Data_processed", "CoverData", "Predictions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# predict and write, one scenario at a time ---------------------------------

walk(scenarios, \(scenario) {
  
  # normals and soils in original units, on the snap grid
  rast <- read_climate_raster(scenario) |>
    align_raster(mask_r) |>
    terra::mask(mask_r)
  
  pred <- predict_raster(fit, rast)
  
  p_out <- file.path(out_dir,
                     paste0(opt$cover_type, "_", opt$cover_model, "_",
                            opt$vc, "-", opt$vmc, "_", scenario, ".tif"))
  # float for both layers: a GeoTIFF has one data type for all bands, and
  # the logical class layer otherwise makes terra write integers, which
  # truncates every probability to 0
  terra::writeRaster(pred, p_out, overwrite = TRUE, datatype = "FLT4S")
})
