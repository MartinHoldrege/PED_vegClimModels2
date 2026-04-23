# Purpose:
# Create wall-to-wall CONUS prediction rasters from fitted cwexp models.
# For each model (herb/woody), produces a multi-layer GeoTIFF with:
#   - total predicted biomass
#   - per-group cover-weighted biomass
#   - per-group potential biomass (cover = 100%)
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

vd <- opt$vd
vm <- opt$vm
vp <- opt$vp
model_type <- opt$model_type  # "herb" or "woody"

stopifnot(model_type %in% c("herb", "woody"))

# load fitted model -------------------------------------------------------

suffix <- paste0(model_type,'_', vd, '-', vp, '-', vm)

p_fit <- file.path(
  paths$large,
  "Data_processed/BiomassQuantityData/Fit",
  paste0("fitted_model_", suffix, ".rds")
)
stopifnot(file.exists(p_fit))

m <- readRDS(p_fit)

fit <- m$fit
config <- m$config
scale_df <- m$scale
cover_cols <- config$model$cover_cols
pred_vars <- config$model$pred_vars

cat("  Model type:", model_type, "\n")
cat("  Formula:", deparse(config$model$formula), "\n")
cat("  Cover cols:", paste(cover_cols, collapse = ", "), "\n")
cat("  Pred vars:", paste(pred_vars, collapse = ", "), "\n\n")

# load CONUS rasters -------------------------------------------------------

rasters <- load_conus_rasters()
scale_df <- rasters$scale_df

# standardize climate
r_climate_std <- standardize_raster(
  rasters$climate[[pred_vars]],
  scale_df,
  vars = pred_vars
)

# build prediction stack: standardized climate + cover
r_pred_stack <- c(r_climate_std, rasters$cover[[cover_cols]])

# predict ------------------------------------------------------------------


r_total <- predict_raster(fit, rast = r_pred_stack, type = "total")

r_by_group <- predict_raster(fit, r_pred_stack, type = "by_group")

r_potential <- predict_raster(fit, r_pred_stack, type = "potential")

# rename potential layers to distinguish from weighted
names(r_potential) <- paste0(names(r_potential), "_potential")

# stack all layers into one raster
r_all <- c(r_total, r_by_group, r_potential)

cat("  Output layers:", paste(names(r_all), collapse = ", "), "\n")

# save ---------------------------------------------------------------------

out_dir <- file.path(
  paths$large,
  "Data_processed/BiomassQuantityData/Predictions"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_out <- file.path(
  out_dir,
  paste0("predicted_biomass_", suffix, ".tif")
)

cat("Writing prediction raster:", p_out, "\n")
terra::writeRaster(r_all, p_out, overwrite = TRUE)

cat("Done.\n")
