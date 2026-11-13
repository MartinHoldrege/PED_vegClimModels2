# 04_combine_tree.R
#
# Combines the tree cover predictions of four models into one tree cover
# prediction, for current and future climate, and at the pixel-years of the
# component models' data version:
#   forest classification (classification/forest): where to use which model
#   tree cover in forest (cover/tree_forest)
#   tree cover in non-forest (cover/tree_nonforest)
#   zero-tree classification (classification/zero_tree)
# The combination (which version of each model, and the blending width) is
# cover_specs$combined$tree[[--vmc]].
#
# Rule, per cell (combine_tree(), used for both the rasters and the
# pixel-years):
#   1. tree = w * forest prediction + (1 - w) * non-forest prediction, where
#      the weight w comes from the forest model's continuous score (probability
#      of forest, or predicted cover for continuous forest models): 0 below
#      score_low, 1 above score_high, and linear in between through 0.5 at the
#      forest threshold (as Alice's 05_FinalSyntheticModelPreds.Rmd did, with
#      approxfun()). So away from the boundary each model is used as is, and
#      the two are blended near it.
#   score_low and score_high are set from the current-climate predictions: the
#   scores below and above the threshold that take in spec$blend_prop of the
#   CONUS cells on each side of it (e.g. 5% of cells on each side). They are
#   saved and reused for the future scenarios and the pixel-years.
#   2. tree = 0 where the cell is classified non-forest and zero tree (no
#      "tree drizzle" in treeless areas). Applied after blending, as in
#      05_FinalSyntheticModelPreds.Rmd, so zero-tree cells just below the
#      forest threshold are 0 while their neighbours just above it are blended
#      values. A cell is in a class where its score is above the model's
#      threshold, as in the class layers from 03_predict_raster.R.
#
# Pixel-years: every row of the components' data version (all four must share
# one, e.g. c07) that has every predictor, for comparing with observed tree
# cover (read_cover_predictions()). Each component's prediction at a row is
# the one saved with the fitted model (pred_in: in-sample, out-of-bag for
# ranger; pred_oof: out-of-fold) where the row was in that model's training
# data. Elsewhere (e.g. tree_forest, fit only where tree cover is > 10%, at the
# rows <= 10%) it is the full model's prediction, for both pred_in and
# pred_oof, because the fold models aren't saved; so there pred_oof isn't
# out-of-fold. fallback marks the rows where such a prediction gets weight in
# the combined pred_oof. Each model makes its own CV folds, so the
# out-of-fold predictions of different components at a row come from
# different folds, but each comes from a fit that didn't use that row. The
# pixel-years use their observed climate anomalies; the rasters are normals
# (anomalies 0).
# The same is done at the rows of the field data (c04), all of them outside
# the components' data, for tree cover by leaf type (05_combine_tree_all.R).
#
# Inputs:
#   Data_processed/CoverData/Predictions/ (from 03_predict_raster.R):
#     <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif for the four models
#   fitted models (thresholds, saved predictions) - 02_fit_classification.R,
#     02_fit_cover.R
#   cover_clim_soils_<vc>.csv - 08_combine_cover_and_covariates.R
#
# Outputs:
#   Data_processed/CoverData/Predictions/combined_tree_<vmc>_<scenario>.tif
#     layers: cover (combined tree cover, %), weight (w, weight on the forest
#     model)
#   Data_processed/CoverData/Fit/combined_tree_<vmc>_blend.csv
#     the component models and the blending parameters
#   Data_processed/CoverData/Fit/combined_tree_<vmc>_<vc>_pixel_predictions.csv
#     one per data version (e.g. c07 and c04): cell, year, pred_in, pred_oof
#     (combined tree cover, %), fallback
#   (paths from combined_path())
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

# spec --------------------------------------------------------------------

if(FALSE) {
  # parameters for testing/development
  opt <- list(cover_type = 'combined',
              cover_model = 'tree',
              vmc = 'm24.3')
}

stopifnot(opt$cover_type == "combined", opt$cover_model == "tree")
spec <- cover_specs$combined$tree[[opt$vmc]]
if (is.null(spec)) stop("no spec for cover_specs$combined$tree$", opt$vmc)

comp <- spec$components  # one row per component model
print(comp)

scenarios <- c("current", "BNU-ESM", "IPSL-CM5A-MR")

# read --------------------------------------------------------------------

#' Prediction raster of one component model for one scenario
read_component <- function(component, scenario) {
  m <- comp[comp$component == component, ]
  terra::rast(cover_pred_path(m$cover_type, m$cover_model, m$vc, m$vmc,
                              scenario))
}

# fitted component models (thresholds, and the predictions at pixel-years)
fits <- map(set_names(comp$component), \(k) {
  m <- comp[comp$component == k, ]
  read_cover_model(m$cover_type, m$cover_model, m$vc, m$vmc)
})
threshold <- fits$forest$threshold
zero_threshold <- fits$zero_tree$threshold
stopifnot(!is.null(threshold), !is.null(zero_threshold))
# "prob", or "cover" for continuous models
score_lyr <- score_name(fits$forest)
zero_score_lyr <- score_name(fits$zero_tree)

# blending parameters, from current climate ------------------------------------

score_cur <- read_component("forest", "current")[[score_lyr]]
v <- terra::values(score_cur, mat = FALSE, na.rm = TRUE)

# share of cells at or below the threshold (non-forest), and the scores that
# take in blend_prop more cells on each side
q_threshold <- mean(v <= threshold) # quantile of the forest break
probs <- q_threshold + c(-1, 1) * spec$blend_prop
stopifnot(probs[1] > 0, probs[2] < 1)
score_cut <- quantile(v, probs, names = FALSE)
rm(v)

blend <- comp |>
  mutate(forest_threshold = threshold,
         forest_score = score_lyr,
         blend_prop = spec$blend_prop,
         q_threshold = q_threshold,
         score_low = score_cut[1],
         score_high = score_cut[2])
print(select(blend, forest_threshold:score_high) |> distinct())
stopifnot(blend$score_low[1] < threshold, threshold < blend$score_high[1])

blend_file <- combined_tree_path(opt$vmc, "blend")
dir.create(dirname(blend_file), recursive = TRUE, showWarnings = FALSE)
write_csv(blend, blend_file)

#' Weight on the forest model: 0 below score_low, 1 above score_high, linear
#' through 0.5 at the threshold in between
blend_weight <- function(score) {
  stats::approx(x = c(score_cut[1], threshold, score_cut[2]),
                y = c(0, 0.5, 1), xout = score, rule = 2)$y
}

#' Combined tree cover from the four component predictions
#'
#' The rule in the header, for vectors (raster values or pixel-years). Uses
#' threshold, zero_threshold and blend_weight() from above.
#'
#' @param forest_score,zero_score Scores of the forest and zero-tree models
#'   (probability, or % cover for continuous models).
#' @param tree_forest,tree_nonforest Tree cover (%) predicted by the forest
#'   and non-forest cover models.
#' @return List: cover (combined tree cover, %), weight (w, weight on the
#'   forest model) and zeroed (TRUE where step 2 set cover to 0).
combine_tree <- function(forest_score, zero_score, tree_forest,
                         tree_nonforest) {
  stopifnot(length(forest_score) == length(zero_score),
            length(forest_score) == length(tree_forest),
            length(forest_score) == length(tree_nonforest))

  # 1. blend
  w <- blend_weight(forest_score)
  cover <- w * tree_forest + (1 - w) * tree_nonforest

  # 2. no trees in non-forest, zero-tree cells
  zeroed <- forest_score <= threshold & zero_score > zero_threshold
  cover[which(zeroed)] <- 0

  list(cover = cover, weight = w, zeroed = zeroed)
}

# combine and write, one scenario at a time --------------------------------

walk(scenarios, \(scenario) {

  forest_r <- read_component("forest", scenario)
  zero_r   <- read_component("zero_tree", scenario)
  tree_f   <- read_component("tree_forest", scenario)[["cover"]]
  tree_nf  <- read_component("tree_nonforest", scenario)[["cover"]]
  stopifnot(terra::compareGeom(forest_r, zero_r, tree_f, tree_nf))

  res <- combine_tree(
    forest_score = terra::values(forest_r[[score_lyr]], mat = FALSE),
    zero_score = terra::values(zero_r[[zero_score_lyr]], mat = FALSE),
    tree_forest = terra::values(tree_f, mat = FALSE),
    tree_nonforest = terra::values(tree_nf, mat = FALSE))

  out <- terra::rast(tree_f, nlyrs = 2)
  terra::values(out) <- cbind(res$cover, res$weight)
  names(out) <- c("cover", "weight")
  terra::writeRaster(out, combined_tree_path(opt$vmc, "prediction", scenario),
                     overwrite = TRUE, datatype = "FLT4S")
  message("Wrote ", scenario)
})

# predictions at the pixel-years --------------------------------------------

# the components' data version (see header), and the field data (c04), for
# tree cover by leaf type (05_combine_tree_all.R)
stopifnot(n_distinct(comp$vc) == 1)
pixel_vcs <- unique(c(comp$vc, "c04"))

source_vars <- map(fits, \(f) f$config$spec$pred_vars) |>
  unlist() |>
  str_remove("^log1p_") |>
  unique()

walk(pixel_vcs, \(pred_vc) {

  # rows with every predictor of every component
  dat <- read_cover_training(pred_vc)
dat <- dat[complete.cases(dat[source_vars]), ]
stopifnot(!anyDuplicated(dat[c("cell", "year")]))

  # component_predictions(): Functions/models/classification_helpers.R
preds <- map(fits, \(f) component_predictions(f, dat))

  # combined tree cover from the in-sample and from the out-of-fold
  # predictions
comb <- map(c(pred_in = "pred_in", pred_oof = "pred_oof"), \(col) {
  combine_tree(forest_score = preds$forest[[col]],
               zero_score = preds$zero_tree[[col]],
               tree_forest = preds$tree_forest[[col]],
               tree_nonforest = preds$tree_nonforest[[col]])
})

  # rows where a full-model (not out-of-fold) prediction gets weight in
  # pred_oof (the cover models get none where step 2 set cover to 0)
w <- comb$pred_oof$weight
used <- !comb$pred_oof$zeroed
fallback <- !preds$forest$saved | !preds$zero_tree$saved |
  (!preds$tree_forest$saved & w > 0 & used) |
  (!preds$tree_nonforest$saved & w < 1 & used)

pixel_pred <- tibble(cell = dat$cell, year = dat$year,
                     pred_in = comb$pred_in$cover,
                     pred_oof = comb$pred_oof$cover,
                     fallback = fallback)

  message(nrow(pixel_pred), " pixel-years of ", pred_vc, "; fallback ",
          "(pred_oof partly from a full-model prediction) at ",
        round(100 * mean(fallback), 1), "%")

write_csv(pixel_pred, combined_path("tree", opt$vmc, "pixel_predictions",
                                    vc = pred_vc))
})