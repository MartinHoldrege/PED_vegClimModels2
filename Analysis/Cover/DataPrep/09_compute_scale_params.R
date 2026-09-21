# 09_compute_scale_params.R
#
# Scaling parameters (mean, sd) for normalizing predictors, so every model
# (cover and biomass) uses the same constants. Read with read_scale_params();
# applied by read_cover_training(normalize = TRUE).
#
# Climate normals and soils: mean and sd across all CONUS cells (snap mask) of
# the 1991-2020 Daymet normals raster and the soils raster. Not LCMAP- or
# fire-masked. Written once and not overwritten unless rerun = TRUE, because
# changing them changes the scale of every fitted model.
#
# Anomalies (_3yrAnom): sd across all training pixel-years (all sources). mean
# is set to 0, so normalizing only divides by the sd and 0 stays "normal". The
# observed mean is stored as mean_obs, for reference.
#
# Inputs:
#   DaymetClimateData_1991-2020_CLIM.tif - 01b_summarise_daymet_climate_data.R
#                                          (read_climate_raster())
#   soil_covariates_solus100_1000m.tif   - 02_soils_calculate_variables.R
#   daymet_conus_snap_1000m.tif          - 00_create_snap_raster.R
#                                          (read_mask())
#   cover_clim_soils_<vc>.csv            - 08_combine_cover_and_covariates.R
#
# Outputs (Data_processed/scale_params/):
#   scale_params_climate_soils.csv
#   scale_params_anomalies_<vc>.csv
#
# September, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

vc <- opt$vc  # cover data version

# params ------------------------------------------------------------------

rerun <- FALSE  # TRUE recomputes the climate/soils constants

out_dir <- file.path(paths$large, "Data_processed/scale_params")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
p_global <- file.path(out_dir, "scale_params_climate_soils.csv")
p_anom   <- file.path(out_dir, paste0("scale_params_anomalies_", vc, ".csv"))

# climate normals and soils (CONUS) ---------------------------------------

if (!file.exists(p_global) || rerun) {
  mask_r <- read_mask()
  # short names (climate_name_lookup()), soils included
  r <- read_climate_raster("current") |>
    align_raster(mask_r) |>
    terra::mask(mask_r)
  
  stats <- terra::global(r, c("mean", "sd", "notNA"), na.rm = TRUE)
  
  global <- tibble(variable = names(r),
                   mean = stats$mean,
                   sd = stats$sd,
                   n = stats$notNA,
                   source = "CONUS raster")
  
  stopifnot(!anyDuplicated(global$variable), all(global$sd > 0))
  write_csv(global, p_global)
}

# anomalies (training pixel-years) ----------------------------------------

anom <- read_cover_training(vc) |>
  select(ends_with("_3yrAnom"))

anom_params <- tibble(variable = names(anom),
                      mean = 0,
                      sd = map_dbl(anom, \(z) sd(z, na.rm = TRUE)),
                      n = colSums(!is.na(anom)),
                      source = paste("training pixel-years", vc),
                      mean_obs = map_dbl(anom, \(z) mean(z, na.rm = TRUE)))

stopifnot(ncol(anom) > 0, all(anom_params$sd > 0))
write_csv(anom_params, p_anom)


