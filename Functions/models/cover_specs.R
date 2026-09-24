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
#   pred_vars     predictors (short climate/soil names). 
#   log1p_vars    predictors to add as log1p (e.g. MAP means log1_MAP get used)
#   squares       add squared terms
#   interactions  add pairwise interactions
#   alpha         glmnet penalty (1 = lasso); ignored by ranger
#   cv            list: cluster_vars, k_clusters, select_rule
#
# classification specs add:
#   cover_threshold   % cover dividing the two classes
#   threshold_method  rule for turning probabilities into a class
#
# ranger specs add:
#   ranger        list of arguments passed on to ranger::ranger()


cover_specs <- list(
  
  classification = list(
    
    # forest / non-forest: tree cover above or below 10%
    forest = list(
      m01 = list(
        response = "cov_tree",
        rows = "all",
        engine = "glmnet",
        cover_threshold = 10,  # percent, as stored
        threshold_method = "PredPrev=Obs",  # see ?PresenceAbsence::optimal.thresholds
        pred_vars = c("MAT", "MAP", "PrecipTempCorr", "isothermality", 'WD_p95',
                      "clay_surface", "awc"), 
        log1p_vars = c('MAP', 'clay_surface', 'awc'),
        squares = TRUE,
        interactions = TRUE,
        alpha = 1,
        cv = list(cluster_vars = c("MAT", "MAP", "PrecipTempCorr",
                                   'awc'),
                  k_clusters = 10,
                  select_rule = "min")  # "1se" or "min"
      ),
      m02 = NULL, # adjusted below
      m03 = NULL, # adjusted below
      m04 = NULL  # adjusted below
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

# mo2 elastic net
cover_specs$classification$forest$m02 <- cover_specs$classification$forest$m01
cover_specs$classification$forest$m02$alpha = 0.5 # elastic net version

# m03 less climate extrapolation
cover_specs$classification$forest$m03 <- cover_specs$classification$forest$m01
cover_specs$classification$forest$m03$cv$k_clusters <- 50

# m04 random forest, as a benchmark for how well the same predictors can do
# without the constraint of an equation. No log1p, squares or interactions:
# tree splits are invariant to monotone transforms, and the forest finds
# interactions itself. Not a candidate for prediction under future climate,
# which falls outside the training range.
cover_specs$classification$forest$m04 <- cover_specs$classification$forest$m01
cover_specs$classification$forest$m04$engine <- "ranger"
cover_specs$classification$forest$m04$log1p_vars <- character(0)
cover_specs$classification$forest$m04$squares <- FALSE
cover_specs$classification$forest$m04$interactions <- FALSE
cover_specs$classification$forest$m04$ranger <- list(num.trees = 300,
                                                     min.node.size = 100)

#' Append log1p-transformed predictor names to every spec's pred_vars
#'
#' @param x A spec list or nested list of specs.
#' @return The same structure with `pred_vars` expanded.
expand_log1p_preds <- function(x) {
  if (!is.list(x)) return(x)
  if (!is.null(x$pred_vars)) {
    x$pred_vars <- c(x$pred_vars, paste0("log1p_", x$log1p_vars, recycle0 = TRUE))
    return(x)
  }
  lapply(x, expand_log1p_preds)
}

cover_specs <- expand_log1p_preds(cover_specs)