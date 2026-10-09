# 04_combine_tree.R
#
# Combines the tree cover predictions of four models into one tree cover
# prediction, for current and future climate:
#   forest classification (classification/forest): where to use which model
#   tree cover in forest (cover/tree_forest)
#   tree cover in non-forest (cover/tree_nonforest)
#   zero-tree classification (classification/zero_tree)
# The combination (which version of each model, and the blending width) is
# cover_specs$combined$tree[[--vmc]].
#
# Rule, per cell:
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
#   saved and reused for the future scenarios.
#   2. tree = 0 where the cell is classified non-forest and zero tree (no
#      "tree drizzle" in treeless areas). Applied after blending, as in
#      05_FinalSyntheticModelPreds.Rmd, so zero-tree cells just below the
#      forest threshold are 0 while their neighbours just above it are blended
#      values.
#
# Inputs (Data_processed/CoverData/Predictions/, from 03_predict_raster.R):
#   <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif for the four models
#   fitted forest model (threshold) - 02_fit_classification.R
#
# Outputs:
#   Data_processed/CoverData/Predictions/combined_tree_<vmc>_<scenario>.tif
#     layers: cover (combined tree cover, %), weight (w, weight on the forest
#     model)
#   Data_processed/CoverData/Fit/combined_tree_<vmc>_blend.csv
#     the component models and the blending parameters
#   (paths from combined_tree_path())
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

# spec --------------------------------------------------------------------

if(FALSE) {
  # parameters for testing/development
  opt <- list(cover_type == 'combined', 
              cover_model == 'tree',
              vmc = 'm01.0')
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

f <- comp[comp$component == "forest", ]
forest_fit <- read_cover_model(f$cover_type, f$cover_model, f$vc, f$vmc)
threshold <- forest_fit$threshold
score_lyr <- score_name(forest_fit)  # "prob", or "cover" for continuous models

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

# combine and write, one scenario at a time --------------------------------

walk(scenarios, \(scenario) {

  forest_r <- read_component("forest", scenario)
  zero_r   <- read_component("zero_tree", scenario)
  tree_f   <- read_component("tree_forest", scenario)[["cover"]]
  tree_nf  <- read_component("tree_nonforest", scenario)[["cover"]]
  stopifnot(terra::compareGeom(forest_r, zero_r, tree_f, tree_nf))

  # 1. blend
  w <- terra::app(forest_r[[score_lyr]], blend_weight)
  tree <- w * tree_f + (1 - w) * tree_nf

  # 2. no trees in non-forest, zero-tree cells
  tree <- terra::ifel(forest_r$class == 0 & zero_r$class == 1, 0, tree)

  out <- c(tree, w)
  names(out) <- c("cover", "weight")
  terra::writeRaster(out, combined_tree_path(opt$vmc, "prediction", scenario),
                     overwrite = TRUE, datatype = "FLT4S")
  message("Wrote ", scenario)
})
