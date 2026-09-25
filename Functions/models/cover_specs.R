# Specifications for the cover models.
#
# Nested: family -> response -> model version. The family decides which
# fitting script runs, the response names what is modelled, and the version
# (m01, m02, ...) is one set of decisions for that response.
#
#   classification  binomial lasso           02_fit_classification.R
#   cover           continuous cover (beta)  03_fit_cover.R
#   proportion      shares of a total (beta) 04_fit_proportion.R
#
# Every spec has:
#   response      column in the cover training data
#   rows          which pixel-years to fit to (see the fitting script)
#   engine        "glmnet" (penalized regression) or "ranger" (random forest)
#   pred_vars     predictors (short climate/soil names), with the log1p_
#                 versions of log1p_vars appended by the constructor
#   log1p_vars    predictors to add as log1p (e.g. MAP means log1p_MAP is used)
#   squares       add squared terms
#   interactions  add pairwise interactions
#   alpha         glmnet penalty (1 = lasso); ignored by ranger
#   ranger        list of arguments passed on to ranger::ranger(); NULL for
#                 glmnet
#   cv            list: cluster_vars, k_clusters, select_rule
#
# classification specs add:
#   cover_threshold   % cover dividing the two classes
#   threshold_method  rule for turning probabilities into a class
#
# Specs are built with a constructor per model (e.g. .defaults_class_forest()),
# so each version lists only what differs from the defaults.


.default_cluster_vars <- c("MAT", "MAP", "PrecipTempCorr", "awc")

#' Spec for the forest / non-forest classification
#'
#' Arguments override the defaults. Appends the `log1p_` names of
#' `log1p_vars` to `pred_vars`, and nests the CV settings under `cv`.
#'
#' @param response Column in the cover training data.
#' @param rows Which pixel-years to fit to.
#' @param engine "glmnet" or "ranger".
#' @param cover_threshold Percent cover dividing the two classes, as stored.
#' @param threshold_method See `?PresenceAbsence::optimal.thresholds`.
#' @param pred_vars Predictors, in original units.
#' @param log1p_vars Predictors to also add as `log1p_<var>`.
#' @param squares,interactions Add squared terms / pairwise interactions.
#' @param alpha glmnet penalty mixing (1 = lasso).
#' @param ranger List of arguments for `ranger::ranger()`; required when
#'   `engine = "ranger"`.
#' @param cluster_vars,k_clusters Variables and number of environmental
#'   clusters for the CV folds.
#' @param select_rule "min" or "1se".
#' @return A spec list.
.defaults_class_forest <- function(response = "cov_tree",
                                   rows = "all",
                                   engine = "glmnet",
                                   cover_threshold = 10,
                                   threshold_method = "PredPrev=Obs",
                                   pred_vars = c("MAT", "MAP", "PrecipTempCorr",
                                                 "isothermality", "WD_p95",
                                                 "clay_surface", "awc"),
                                   log1p_vars = c("MAP", "clay_surface", "awc"),
                                   squares = TRUE,
                                   interactions = TRUE,
                                   alpha = 1,
                                   ranger = NULL,
                                   cluster_vars = .default_cluster_vars,
                                   k_clusters = 10,
                                   select_rule = "min") {
  stopifnot(engine %in% c("glmnet", "ranger"),
            select_rule %in% c("min", "1se"),
            engine != "ranger" || is.list(ranger))
  
  list(
    response = response,
    rows = rows,
    engine = engine,
    cover_threshold = cover_threshold,
    threshold_method = threshold_method,
    # recycle0: with no log1p_vars, adds nothing rather than "log1p_"
    pred_vars = c(pred_vars, paste0("log1p_", log1p_vars, recycle0 = TRUE)),
    log1p_vars = log1p_vars,
    squares = squares,
    interactions = interactions,
    alpha = alpha,
    ranger = ranger,
    cv = list(cluster_vars = cluster_vars,
              k_clusters = k_clusters,
              select_rule = select_rule)
  )
}


cover_specs <- list(
  
  classification = list(
    
    # forest / non-forest: tree cover above or below 10%
    forest = list(
      
      m01 = .defaults_class_forest(),
      
      # elastic net
      m02 = .defaults_class_forest(alpha = 0.5),
      
      # less climate extrapolation
      m03 = .defaults_class_forest(k_clusters = 50),
      
      # random forest, as a benchmark for how well the same predictors can do
      # without the constraint of an equation. No log1p, squares or
      # interactions: tree splits are invariant to monotone transforms, and
      # the forest finds interactions itself. Not a candidate for prediction
      # under future climate, which falls outside the training range.
      # min.node.size is the smallest node that can be split, not the smallest
      # leaf. min.bucket (smallest leaf) was tried and roughly tripled fit time
      m04 = .defaults_class_forest(engine = "ranger",
                                   log1p_vars = character(0),
                                   squares = FALSE,
                                   interactions = FALSE,
                                   ranger = list(num.trees = 300,
                                                 min.node.size = 100))
    )
    
    # zero_tree: trees vs no trees in non-forest, trained on a binarized
    # RAP raster. Needs its own data-reading function; not specified yet.
  ),
  
  # continuous cover: tree cover in forest and in non-forest, herbaceous,
  # shrub. Not specified yet.
  cover = list(),
  
  # shares: needleleaf (of tree), forb / C3 / C4 (of herbaceous, rescaled to
  # sum to 1 at prediction time). Not specified yet.
  proportion = list()
)