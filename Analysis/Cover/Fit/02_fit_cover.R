# 03_fit_cover.R
#
# Fits one continuous cover model (tree, shrub or herbaceous cover) and saves
# it. The counterpart of 02_fit_classification.R, sharing its folds and
# fitting engines (fit_engine.R), but the response is % cover and there is no
# class or threshold.
#
# The model is chosen by --cover_type (here "cover"), --cover_model (e.g.
# tree_forest) and --vmc, which index cover_specs; --vc is the cover data
# version. Defaults are in params.R; main_cover.R can run other combinations.
#
# spec$rows picks the pixel-years to fit to: tree cover is fit separately
# above and at or below 10% (cover_specs$cover$tree_forest and
# tree_nonforest), to be blended at 10% after prediction.
#
# spec$family: "quasibinomial" (glmnet or gam): cover / 100 as a proportion,
# predicting mean cover (%). "gaussian": cover after spec$response_transform;
# with engine "ranger", a regression forest.
#
# Predictors: as in 02_fit_classification.R (prepare_predictors()).
#
# Inputs:
#   cover_clim_soils_<vc>.csv  - 08_combine_cover_and_covariates.R, or
#                                05_rap_training_sample.R (c02)
#   scale_params_*.csv         - 09_compute_scale_params.R
#   cover_specs                - Functions/models/cover_specs.R
#
# Output:
#   Data_processed/CoverData/Fit/<cover_type>_<cover_model>_<vc>-<vmc>.rds
#   (path from cover_model_path()), class "cover_continuous"
#
# October, 2026


# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")

library(glmnet)

# spec --------------------------------------------------------------------
test_run <- FALSE
spec <- cover_specs[[opt$cover_type]][[opt$cover_model]][[opt$vmc]]

if (is.null(spec)) {
  stop("no spec for cover_specs$", opt$cover_type, "$", opt$cover_model, "$",
       opt$vmc)
}
stopifnot(opt$cover_type == "cover")

out_file <- cover_model_path(opt$cover_type, opt$cover_model, opt$vc,
                             opt$vmc)

# data --------------------------------------------------------------------

# original units; prepare_predictors() does the transforming and scaling
dat <- read_cover_training(opt$vc)

# which pixel-years to fit to (e.g. tree cover above 10%)
dat <- select_rows(dat, spec$rows)

# rows with the response and every source column; log1p_MAP comes from MAP
source_vars <- unique(c(str_remove(spec$pred_vars, "^log1p_"),
                        spec$cv$cluster_vars))
dat <- dat[complete.cases(dat[c(spec$response, source_vars)]), ]

if (test_run) {
  dat <- sample_n(dat, 1000)
}

cover <- dat[[spec$response]]
stopifnot(all(cover >= 0 & cover <= 100))

# what the model is fit to
y_fit <- switch(spec$family,
                quasibinomial = cover / 100,
                gaussian = switch(spec$response_transform,
                                  identity = cover,
                                  log1p = log1p(cover)),
                stop("family not set up for cover models: ", spec$family))

scale_df <- read_scale_params()  # climate and soils; no anomalies here

x <- prepare_predictors(dat,
                        pred_vars = spec$pred_vars,
                        scale_df = scale_df,
                        squares = spec$squares,
                        interactions = spec$interactions,
                        interact_log1p = spec$interact_log1p)
stopifnot(nrow(x) == length(y_fit))

# folds -------------------------------------------------------------------

# environmental blocking (see make_cover_folds())
folds <- make_cover_folds(dat, spec)
foldid <- folds$foldid
clusters <- folds$clusters
rm(folds)

# fit ---------------------------------------------------------------------

# engine and family from the spec (see fit_engine()). Predictions are % cover.
# pred_in: in-sample (out-of-bag for ranger); pred_oof: out-of-fold
res <- fit_engine(spec, x = x, y_fit = y_fit, foldid = foldid)

# save --------------------------------------------------------------------

out <- list(
  fit = res$fit,
  lambda = res$lambda,
  config = list(cover_type = opt$cover_type,
                cover_model = opt$cover_model,
                vc = opt$vc, vmc = opt$vmc, spec = spec,
                x_colnames = colnames(x)),
  scale_df = scale_df,
  clustering = clusters,  # centers and scaling, for assign_to_clusters()
  predictions = tibble(cell = dat$cell, year = dat$year,
                       pred_in = res$pred_in, pred_oof = res$pred_oof)
)

if (!test_run) {
  class(out) <- "cover_continuous"
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, out_file)
}
