# 03_predict_raster.R
#
# Wall-to-wall predictions from one fitted cover model, chosen by
# --cover_type, --cover_model, --vc and --vmc (params.R). Uses
# predict_raster(), which dispatches on the model's class, so the same script
# serves every model type that has a predict_raster method.
#
# Predicts on every cell with climate and soils; masking to the available
# area is left to the maps. The climate rasters are normals, so the anomaly
# predictors are set to 0 (predict_raster()).
#
# --climate_set (params.R; "minimal" unless passed) chooses the climates:
#   minimal  1 km Daymet grid: current (Daymet 1991-2020 normals) and the
#            BNU-ESM and IPSL-CM5A-MR RCP 8.5, 2070-2099 futures
#            (MACAv2-METDATA interpolated onto the Daymet grid;
#            read_maca_daymet(), cached). Before November 2026 the futures
#            came from ForecastedClimateData_*, so predictions differ, and
#            1 km cells next to an NA 4 km cell are now NA.
#   full     current (Daymet, 1 km), as above; gridMET 1991-2020 normals
#            (4 km); and the MACAv2-METDATA futures of every GCM, for each
#            RCP and period in `rcps` and `periods` below (4 km, gridMET
#            grid). gridMET is the historical for the futures (MACA is
#            bias-corrected to gridMET); the GCMs' own historical runs aren't
#            used. read_climate_raster() reads them, with the soils averaged
#            from 1 km onto the gridMET grid (cached; read_soil_gridmet()).
#
# GCM summaries (full set): for each RCP and period, the per-cell min, median
# and max across GCMs, and each GCM's change from gridMET (delta: CONUS mean
# of GCM prediction - gridMET prediction, area-weighted, over the cells where
# both have a value), with the GCMs that are lowest, median and highest by
# delta. Not made for the components of combined models (the tree
# models and the proportion models): the combine scripts make them from the
# combined predictions.
#
# Inputs:
#   fitted model - 02_fit_classification.R or 02_fit_cover.R
#                  (read_cover_model())
#   Data_processed/WallToWallClimateData/ (read_climate_raster()):
#     DaymetClimateData_1991-2020_CLIM.tif (and ForecastedClimateData_*)
#     gridmet/gridmet_1991-2020_CLIM.tif  - 01_summarise_gridmet_climate_data.R
#     maca/maca_<gcm>_<rcp>_<period>_CLIM.tif
#                                    - 01_summarise_projected_climate_data.R
#   soil_covariates_solus100_1000m.tif - read_soil_raster()
#
# Outputs (Data_processed/CoverData/Predictions/; <stem> =
# <cover_type>_<cover_model>_<vc>-<vmc>; .tif names from cover_pred_path()):
#   <stem>_current.tif  layers prob or cover, and class for classification
#                       models
#   minimal set:
#     <stem>_BNU-ESM.tif, <stem>_IPSL-CM5A-MR.tif  as current
#   full set:
#     <stem>_gridmet.tif                as current
#     <stem>_<rcp>_<period>.tif         one layer per GCM (named by GCM): the
#                                       prediction (prob or cover)
#     <stem>_<rcp>_<period>_summary.tif one per RCP and period: layers min,
#                                       median, max across GCMs, per cell
#                                       (not for combined-model components)
#     <stem>_gcm_deltas.csv             per RCP, period and GCM: CONUS means
#                                       of the GCM and gridMET predictions,
#                                       delta, rank (1 = lowest delta) and
#                                       role (lowest, median = rank n/2,
#                                       i.e. 10th of 20, highest) (not for
#                                       combined-model components)
#
# September, 2026 (gridMET and MACA GCMs added November 2026)


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------

climate_set <- opt$climate_set %||% "minimal"
stopifnot(climate_set %in% c("minimal", "full"))

# full set: MACA futures
rcps <- c("rcp45", "rcp85")
periods <- "2070-2099"  # end of century (also available: "2040-2069")
n_gcms <- 20            # GCMs expected for each RCP and period

# components of combined models: GCM ranking and per-cell summaries are made
# by the combine scripts, from the combined predictions
gcm_summaries <- !(opt$cover_type %in% c("classification", "proportion") ||
                     opt$cover_model %in% c("tree_forest", "tree_nonforest"))

clim_dir <- file.path(paths$large, "Data_processed/WallToWallClimateData")
out_dir <- file.path(paths$large, "Data_processed", "CoverData", "Predictions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- paste0(opt$cover_type, "_", opt$cover_model, "_", opt$vc, "-",
               opt$vmc)


if (FALSE) {
  # parameters for testing/development (run these lines by hand, then step
  # through the script)
  opt <- list(cover_type = "cover",
              cover_model = "shrub",  # gets the GCM summaries
              vc = "c08",
              vmc = "m05.51",
              climate_set = "full")   # or "minimal"
  # a combined-model component (no GCM summaries), e.g.:
  # opt <- list(cover_type = "classification", cover_model = "forest",
  #             vc = "c07", vmc = "m01.51", climate_set = "full")
  climate_set <- opt$climate_set
}
# read --------------------------------------------------------------------

fit <- read_cover_model(opt$cover_type, opt$cover_model, opt$vc, opt$vmc)
if(is.null(fit$config$spec$engine)) fit$config$spec$engine <- 'glmnet' # for legacy reasons
mask_r <- read_mask()

# functions ---------------------------------------------------------------

#' Write a prediction raster to cover_pred_path()
#' @param r SpatRaster.
#' @param scenario Scenario part of the file name, e.g. "current".
write_pred <- function(r, scenario) {
  p_out <- cover_pred_path(opt$cover_type, opt$cover_model, opt$vc, opt$vmc,
                           scenario)
  # float for every layer: a GeoTIFF has one data type for all bands, and
  # the logical class layer otherwise makes terra write integers, which
  # truncates every probability to 0
  terra::writeRaster(r, p_out, overwrite = TRUE, datatype = "FLT4S")
}

#' Area-weighted mean of a raster over its non-NA cells
#' @param r Single-layer SpatRaster.
#' @param area SpatRaster of cell areas, on the same grid.
#' @return Numeric.
area_mean <- function(r, area) {
  area <- terra::mask(area, r)
  sum_r_area <- terra::global(r * area, "sum", na.rm = TRUE)[[1]]
  sum_area <- terra::global(area, "sum", na.rm = TRUE)[[1]]
  sum_r_area / sum_area
}

# 1 km Daymet grid ----------------------------------------------------------

scenarios_1km <- "current"
if (climate_set == "minimal") {
  scenarios_1km <- c("current", "BNU-ESM", "IPSL-CM5A-MR")
}

walk(scenarios_1km, \(scenario) {
  # normals and soils in original units, on the snap grid
  rast <- read_climate_raster(scenario) |>
    align_raster(mask_r) |>
    terra::mask(mask_r)
  
  write_pred(predict_raster(fit, rast), scenario)
})

# 4 km gridMET grid (full set) ----------------------------------------------

if (climate_set == "full") {

  # gridMET normals, and soils averaged onto the gridMET grid (cached)
  pred_gridmet <- predict_raster(fit, read_climate_raster("gridmet"))
  write_pred(pred_gridmet, "gridmet")
  score_gridmet <- pred_gridmet[[score_name(fit)]]

  # cell areas, for CONUS means (cells of a lon/lat grid differ in area)
  area <- terra::cellSize(score_gridmet, unit = "km")

  rcp_periods <- expand_grid(rcp = rcps, period = periods)

  deltas <- pmap(rcp_periods, \(rcp, period) {
    # the GCMs with a MACA file for this RCP and period
    files <- list.files(file.path(clim_dir, "maca"),
                        pattern = paste0("^maca_.+_", rcp, "_", period,
                                         "_CLIM\\.tif$"))
    gcms <- str_match(files, paste0("^maca_(.+)_", rcp))[, 2]
    stopifnot(length(gcms) == n_gcms)
    message(rcp, " ", period, ": ", length(gcms), " GCMs")

    # one layer per GCM: the prediction (prob or cover), on the gridMET grid
    preds <- map(gcms, \(g) {
      rast <- read_climate_raster("maca", gcm = g, rcp = rcp, period = period)
      predict_raster(fit, rast)[[score_name(fit)]]
    }) |>
      terra::rast()
    names(preds) <- gcms

    scenario <- paste0(rcp, "_", period)
    write_pred(preds, scenario)

    if (!gcm_summaries) return(NULL)

    # across GCMs, per cell
    summary_r <- c(min(preds), median(preds), max(preds))
    names(summary_r) <- c("min", "median", "max")
    write_pred(summary_r, paste0(scenario, "_summary"))

    # CONUS means over the cells where both the GCM and gridMET predictions
    # have a value
    map(gcms, \(g) {
      delta_r <- preds[[g]] - score_gridmet
      tibble(rcp = rcp, period = period, gcm = g,
             mean_gcm = area_mean(terra::mask(preds[[g]], delta_r), area),
             delta = area_mean(delta_r, area))
    }) |>
      bind_rows() |>
      mutate(rank = rank(delta, ties.method = "first"),
             role = case_when(rank == 1 ~ "lowest",
                              rank == ceiling(n() / 2) ~ "median",
                              rank == n() ~ "highest"))
  })

  if (gcm_summaries) {
    deltas <- bind_rows(deltas)
    print(filter(deltas, !is.na(role)))
    write_csv(deltas, file.path(out_dir, paste0(stem, "_gcm_deltas.csv")))
  }
}