# 04_combine_herb.R
#
# Herbaceous cover by functional type (C3 grass, C4 grass, forb), for current
# and future climate, and at the pixel-years of the share models' data
# version. Combines four models:
#   herbaceous cover (cover/herbaceous)
#   C3, C4 and forb shares of herbaceous cover (proportion/c3, c4, forb)
# The combination (which version of each model) is
# cover_specs$combined$herbaceous[[--vmc]].
#
# The three shares are fit separately, so they don't sum to 100%; they are
# rescaled to do so, and cover of each type is its rescaled share of the
# predicted herbaceous cover (split_cover(), used for both the rasters and the
# pixel-years). The share_sum layer keeps their sum before rescaling, to show
# how much rescaling did.
#
# --climate_set (params.R; "minimal" unless passed) chooses the climates, as
# in 03_predict_raster.R (whose outputs for the four models are read here):
#   minimal  current, BNU-ESM and IPSL-CM5A-MR (1 km)
#   full     current (1 km), gridMET (4 km), and every MACA GCM for each RCP
#            and period in `rcps` and `periods` below (4 km)
# For the GCMs only the cover of each type is written (the shares follow from
# the covers), with what 03_predict_raster.R makes for single models: the
# per-cell min, median and max across GCMs, and each GCM's change from gridMET
# (gcm_deltas()), with the GCMs ranked separately for each type.
#
# Pixel-years: every row of the share models' data version (all three must
# share one, e.g. c04, the only data with observed cover by type) that has
# every predictor, for comparing with the observations
# (read_cover_predictions()). As in 04_combine_tree.R, each component's
# prediction at a row is the one saved with the fitted model where the row was
# in that model's training data (pred_in: in-sample, out-of-bag for ranger;
# pred_oof: out-of-fold), and elsewhere the full model's prediction, for both,
# so there pred_oof isn't out-of-fold (e.g. the herbaceous model, fit to c08,
# at c04 rows not in c08). fallback marks the rows where a share model's
# prediction isn't out-of-fold; the herbaceous model is ignored (fit to other
# data, so at most rows its prediction is a full-model one). The pixel-years
# use their observed climate anomalies; the rasters are normals (anomalies 0).
#
#
# Inputs:
#   Data_processed/CoverData/Predictions/ (from 03_predict_raster.R):
#     <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif for the four models
#   fitted models (saved predictions) - 02_fit_cover.R
#   cover_clim_soils_<vc>.csv - 08_combine_cover_and_covariates.R
#
# Outputs (paths from combined_path()):
#   Data_processed/CoverData/Predictions/combined_herbaceous_<vmc>_<scenario>.tif
#     scenario current, BNU-ESM and IPSL-CM5A-MR (minimal set), or current and
#     gridmet (full set): layers share_c3, share_c4, share_forb (% of
#     herbaceous cover, sum to 100), cover_c3, cover_c4, cover_forb (% cover),
#     share_sum (%)
#   full set, for each RCP and period (same folder and stem):
#     ..._<rcp>_<period>.tif          layers cover_<type>_<gcm> (% cover), e.g.
#                                     cover_c3_CCSM4
#     ..._<rcp>_<period>_summary.tif  layers cover_<type>_min, _median, _max:
#                                     across GCMs, per cell
#     combined_herbaceous_<vmc>_gcm_deltas.csv
#       per RCP, period, layer (cover_<type>) and GCM: CONUS mean of the GCM
#       prediction, delta (from gridMET), rank (1 = lowest delta, within the
#       layer) and role (lowest, median, highest), as in 03_predict_raster.R
#   Data_processed/CoverData/Fit/combined_herbaceous_<vmc>_components.csv
#     the component models
#   Data_processed/CoverData/Fit/combined_herbaceous_<vmc>_<vc>_pixel_predictions.csv
#     cell, year, layer (as the raster layers), pred_in, pred_oof, fallback
#
# October, 2026 (GCMs added November 2026)


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

# spec --------------------------------------------------------------------
if(FALSE) {
  # forest testing
  opt <- list(cover_type = 'combined',
              cover_model = 'herbaceous',
              vmc = 'm21.51',
              climate_set = 'full')
}
stopifnot(opt$cover_type == "combined", opt$cover_model == "herbaceous")
spec <- cover_specs$combined$herbaceous[[opt$vmc]]
if (is.null(spec)) stop("no spec for cover_specs$combined$herbaceous$", opt$vmc)

comp <- spec$components  # one row per component model
print(comp)

components_file <- combined_path("herbaceous", opt$vmc, "components")
dir.create(dirname(components_file), recursive = TRUE, showWarnings = FALSE)
write_csv(comp, components_file)

# params ------------------------------------------------------------------

climate_set <- opt$climate_set %||% "minimal"
stopifnot(climate_set %in% c("minimal", "full"))

# climates with every layer (shares, covers, share_sum)
if (climate_set == "minimal") {
scenarios <- c("current", "BNU-ESM", "IPSL-CM5A-MR")
} else {
  scenarios <- c("current", "gridmet")
}

# full set: MACA futures, as in 03_predict_raster.R
rcps <- c("rcp45", "rcp85")
periods <- "2070-2099"
n_gcms <- 20

types <- c("c3", "c4", "forb")
cover_lyrs <- paste0("cover_", types)  # the layers written for the GCMs

# functions ---------------------------------------------------------------

# prediction of one component for one scenario: % cover, or % of
# herbaceous cover for the shares
read_component <- function(component, scenario) {
  read_component_pred(comp[comp$component == component, ], scenario)
}

# predictions of one component for one RCP and period (scenario
# "<rcp>_<period>"): one layer per GCM, named by GCM
read_component_gcms <- function(component, scenario) {
  m <- comp[comp$component == component, ]
  terra::rast(cover_pred_path(m$cover_type, m$cover_model, m$vc, m$vmc,
                              scenario))
}

write_combined <- function(r, scenario) {
  terra::writeRaster(r, combined_path("herbaceous", opt$vmc, "prediction",
                                      scenario),
                     overwrite = TRUE, datatype = "FLT4S")
  message("Wrote ", scenario)
}

# combine and write, one scenario at a time --------------------------------

walk(scenarios, \(scenario) {
  shares <- map(types, \(type) read_component(type, scenario)) |>
    terra::rast()
  names(shares) <- types
  herb <- read_component("herbaceous", scenario)

  write_combined(split_cover(shares, herb), scenario)
})

# GCMs (full set) -------------------------------------------------------------

if (climate_set == "full") {

  # the historical, for the changes
  gridmet <- terra::rast(combined_path("herbaceous", opt$vmc, "prediction",
                                       "gridmet"))

  rcp_periods <- expand_grid(rcp = rcps, period = periods)

  deltas <- pmap(rcp_periods, \(rcp, period) {
    scenario <- paste0(rcp, "_", period)

    herb <- read_component_gcms("herbaceous", scenario)
    shares <- map(set_names(types), \(type) {
      read_component_gcms(type, scenario)
    })
    gcms <- names(herb)
    stopifnot(length(gcms) == n_gcms,
              map_lgl(shares, \(r) identical(names(r), gcms)))
    message(rcp, " ", period, ": ", length(gcms), " GCMs")

    # cover of each type, one raster per GCM
    cover <- map(set_names(gcms), \(g) {
      shares_g <- terra::rast(map(shares, \(r) r[[g]]))
      names(shares_g) <- types
      split_cover(shares_g, herb[[g]])[[cover_lyrs]]
    })

    # one raster per type, one layer per GCM
    by_type <- map(set_names(cover_lyrs), \(lyr) {
      r <- terra::rast(map(cover, \(x) x[[lyr]]))
      names(r) <- gcms
      r
    })

    # every type in one file, layers cover_<type>_<gcm>
    gcm_r <- terra::rast(unname(by_type))
    names(gcm_r) <- paste0(rep(cover_lyrs, each = length(gcms)), "_", gcms)
    write_combined(gcm_r, scenario)

    # across GCMs, per cell: layers cover_<type>_min, _median, _max
    summary_r <- map(cover_lyrs, \(lyr) {
      r <- by_type[[lyr]]
      s <- c(min(r), terra::median(r), max(r))
      names(s) <- paste0(lyr, c("_min", "_median", "_max"))
      s
    }) |>
      terra::rast()
    write_combined(summary_r, paste0(scenario, "_summary"))

    # change from gridMET, with the GCMs ranked separately for each type
    # (Functions/spatial.R)
    map(cover_lyrs, \(lyr) {
      gcm_deltas(by_type[[lyr]], gridmet[[lyr]]) |>
        mutate(layer = lyr, .before = 1)
    }) |>
      bind_rows() |>
      mutate(rcp = rcp, period = period, .before = 1)
  })

  deltas <- bind_rows(deltas)
  print(filter(deltas, !is.na(role)), n = Inf)
  write_csv(deltas, combined_path("herbaceous", opt$vmc, "gcm_deltas"))
}

# predictions at the pixel-years --------------------------------------------

# fitted component models
fits <- map(set_names(comp$component), \(k) {
  m <- comp[comp$component == k, ]
  read_cover_model(m$cover_type, m$cover_model, m$vc, m$vmc)
})

# the share models' data version (see header)
pred_vc <- unique(comp$vc[comp$component %in% types])
stopifnot(length(pred_vc) == 1)

# rows with every predictor of every component
dat <- read_cover_training(pred_vc)
source_vars <- map(fits, \(f) f$config$spec$pred_vars) |>
  unlist() |>
  str_remove("^log1p_") |>
  unique()
dat <- dat[complete.cases(dat[source_vars]), ]
stopifnot(!anyDuplicated(dat[c("cell", "year")]))

# component_predictions(): Functions/models/classification_helpers.R
preds <- map(fits, \(f) component_predictions(f, dat))

# split_cover() of the in-sample and of the out-of-fold predictions, one row
# per pixel-year and layer
pixel_pred <- map(c("pred_in", "pred_oof"), \(col) {
  shares <- as_tibble(map(preds[types], col))  # columns c3, c4, forb
  split_cover(shares, preds$herbaceous[[col]]) |>
    mutate(cell = dat$cell, year = dat$year) |>
    pivot_longer(-c(cell, year), names_to = "layer", values_to = col)
}) |>
  reduce(left_join, by = c("cell", "year", "layer"))

# pixel-years where a share model's prediction is a full-model (not
# out-of-fold) prediction. The herbaceous model is ignored: it's fit to other
# data (c08), so at most rows it only has a full-model prediction.
# TRUE where any of the share models' predictions isn't out-of-fold
fb <- tibble(cell = dat$cell, year = dat$year,
             fallback = !preds$c3$saved | !preds$c4$saved | !preds$forb$saved)

pixel_pred <- left_join(pixel_pred, fb, by = c("cell", "year"))

message(nrow(dat), " pixel-years of ", pred_vc, "; fallback at ",
        round(100 * mean(fb$fallback), 1), "%")

write_csv(pixel_pred, combined_path("herbaceous", opt$vmc,
                                    "pixel_predictions", vc = pred_vc))