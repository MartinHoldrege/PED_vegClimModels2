# 09_compute_scale_params.R
#
# Scaling parameters (mean, sd) for normalizing predictors, so every model
# (cover and biomass) uses the same constants. Read with read_scale_params();
# applied by read_cover_training(normalize = TRUE).
#
# Climate normals and soils: mean and sd across all CONUS cells (snap mask) of
# the 1991-2020 Daymet normals raster and the soils raster, plus the same for
# the log1p version of variables that could reasonably be logged. Not LCMAP-
# or fire-masked. Written once and not overwritten unless rerun = TRUE, because
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

rerun <- TRUE # FALSE  # recomputes the climate/soils constants

# a log1p version this correlated with the original is a near-duplicate, so
# no parameters are stored for it
log_cor_max <- 0.99

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
  
  # log1p versions, for models that log a predictor. Logging happens in
  # original units, before scaling, so each needs its own mean and sd.
  # log1p (log(x + 1)) is used rather than log, so a variable that is 0
  # somewhere (driest-month precipitation, wet degree days) still works.
  # Candidates are non-negative with a meaningful zero, including all the
  # soil variables; signed variables (water deficit, PrecipTempCorr) are
  # never logged. log1p rather than log so nothing fails on a zero in any
  # dataset the equations are later applied to (future climate, gridMET). Computing a few that no
  # model uses is harmless.
  log_vars <- c("MAP", "P_wettestMonth", "P_driestMonth",
                      "P_seasonality", "WDD_mean", "WDD_p05",
                      "frost_free_days", "frost_free_days_p05",
                      "VPD_mean", "VPD_max", "VPD_max_p95",
                      "soilDepth", "clay_surface", "clay", "sand", "coarse",
                      "carbon", "awc")
  
  stopifnot(log_vars %in% names(r))

  # these are non-negative by construction; log1p needs values above -1
  mins <- terra::global(r[[log_vars]], "min", na.rm = TRUE)$min
  stopifnot(all(mins >= 0))
  
  # A log1p version that is nearly a straight line rescaling of the original
  # adds nothing to a model, so drop it. Correlation is measured on a sample
  # of CONUS cells, i.e. over the range the variable actually takes.
  set.seed(1)
  samp <- terra::spatSample(r[[log_vars]], size = 1e5, method = "random",
                            na.rm = TRUE)
  
  log_cor <- map_dbl(log_vars, \(v) cor(samp[[v]], log1p(samp[[v]])))
  names(log_cor) <- log_vars
  
  message("correlation of each variable with its log1p version:")
  print(round(sort(log_cor), 3))
  
  dropped <- log_vars[log_cor > log_cor_max]
  log_vars <- log_vars[log_cor <= log_cor_max]
  if (length(dropped) > 0) {
    message("not logged (r > ", log_cor_max, " with the original): ",
            paste(dropped, collapse = ", "))
  }
  
  r_log <- log1p(r[[log_vars]])
  names(r_log) <- paste0("log1p_", log_vars)
  
  stats_log <- terra::global(r_log, c("mean", "sd", "notNA"), na.rm = TRUE)
  
  global_log <- tibble(variable = names(r_log),
                       mean = stats_log$mean,
                       sd = stats_log$sd,
                       n = stats_log$notNA,
                       source = "CONUS raster (log1p)")
  
  global <- bind_rows(global, global_log)
  
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