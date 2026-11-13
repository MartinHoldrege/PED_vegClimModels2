# 05_combine_tree_all.R
#
# Needleleaf and broadleaf tree cover at the pixel-years of the needleleaf
# share model's data version (e.g. c04, the field data, the only data with
# observed cover by leaf type), for comparing with the observations
# (read_cover_predictions()). Combines the two components of
# cover_specs$combined$tree_all[[--vmc]]:
#   combined tree cover: its predictions at these pixel-years, from
#     04_combine_tree.R (run that first)
#   needleleaf share of tree cover (proportion/needle)
# Tree cover is split with split_cover(): needleleaf is the share, broadleaf
# the rest, as tree_leaf_cover() does for the rasters. The rasters aren't
# written; tree_leaf_cover() computes them on the fly.
#
# The needleleaf share at a row is the one saved with the fitted model where
# the row was in its training data, and elsewhere the full model's prediction
# (component_predictions()). fallback marks the rows where the needleleaf
# share isn't out-of-fold; the combined tree cover is ignored (its components
# are fit to other data, e.g. c07, so at these rows its predictions are
# full-model ones).
#
# Inputs:
#   combined_tree_<tree vmc>_<vc>_pixel_predictions.csv - 04_combine_tree.R
#   fitted needleleaf share model - 02_fit_cover.R
#   cover_clim_soils_<vc>.csv - 08_combine_cover_and_covariates.R
#
# Output (path from combined_path()):
#   Data_processed/CoverData/Fit/combined_tree_all_<vmc>_<vc>_pixel_predictions.csv
#     cell, year, layer (share_needle, share_broad: % of tree cover;
#     cover_needle, cover_broad: % cover), pred_in, pred_oof, fallback
#
# November, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# spec --------------------------------------------------------------------
if(FALSE) {
  opt$cover_type = "combined"
  opt$cover_model = "tree_all"
  opt$vmc = 'm21.5'
}
stopifnot(opt$cover_type == "combined", opt$cover_model == "tree_all")
spec <- cover_specs$combined$tree_all[[opt$vmc]]
if (is.null(spec)) stop("no spec for cover_specs$combined$tree_all$", opt$vmc)

comp <- spec$components  # one row per component model
print(comp)

tree_m <- comp[comp$component == "tree", ]
needle_m <- comp[comp$component == "needle", ]
pred_vc <- needle_m$vc  # the needleleaf share model's data version

# read --------------------------------------------------------------------

needle_fit <- read_cover_model(needle_m$cover_type, needle_m$cover_model,
                               needle_m$vc, needle_m$vmc)

# combined tree cover at the pixel-years of pred_vc (04_combine_tree.R)
tree_pred <- read_cover_predictions("combined", "tree", pred_vc, tree_m$vmc)

# rows with combined tree cover and every predictor of the needleleaf model
dat <- read_cover_training(pred_vc) |>
  semi_join(tree_pred, by = c("cell", "year"))
source_vars <- unique(str_remove(needle_fit$config$spec$pred_vars, "^log1p_"))
dat <- dat[complete.cases(dat[source_vars]), ]
stopifnot(!anyDuplicated(dat[c("cell", "year")]))

# predictions at the pixel-years -------------------------------------------

needle <- component_predictions(needle_fit, dat)
tree <- left_join(select(dat, cell, year), tree_pred, by = c("cell", "year"))

# split_cover() of the in-sample and of the out-of-fold predictions, one row
# per pixel-year and layer
pixel_pred <- map(c("pred_in", "pred_oof"), \(col) {
  share_needle <- needle[[col]]
  shares <- tibble(needle = share_needle, broad = 100 - share_needle)
  split_cover(shares, tree[[col]]) |>
    select(-share_sum) |>  # the two shares sum to 100 already
    mutate(cell = dat$cell, year = dat$year) |>
    pivot_longer(-c(cell, year), names_to = "layer", values_to = col)
}) |>
  reduce(left_join, by = c("cell", "year", "layer"))

# pixel-years where the needleleaf share isn't out-of-fold (see header)
fb <- tibble(cell = dat$cell, year = dat$year, fallback = !needle$saved)
pixel_pred <- left_join(pixel_pred, fb, by = c("cell", "year"))

message(nrow(dat), " pixel-years of ", pred_vc, "; fallback at ",
        round(100 * mean(fb$fallback), 1), "%")

write_csv(pixel_pred, combined_path("tree_all", opt$vmc, "pixel_predictions",
                                    vc = pred_vc))
