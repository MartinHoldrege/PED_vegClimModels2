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
#   pred_vars     predictors (short climate/soil names). A name starting
#                 "log1p_" (e.g. "log1p_MAP") is made from the original
#                 variable before standardizing; parameters for it must exist
#                 in 09_compute_scale_params.R output
#   squares       add squared terms
#   interactions  add pairwise interactions
#   alpha         glmnet penalty (1 = lasso)
#   cv            list: cluster_vars, k_clusters, select_rule
#
# classification specs add:
#   cover_threshold   % cover dividing the two classes
#   threshold_method  rule for turning probabilities into a class


cover_specs <- list(
  
  classification = list(
    
    # forest / non-forest: tree cover above or below 10%
    forest = list(
      m01 = list(
        response = "cov_tree",
        rows = "all",
        cover_threshold = 10,  # percent, as stored
        threshold_method = "PredPrev=Obs",  # see ?PresenceAbsence::optimal.thresholds
        pred_vars = c("MAT", "log1p_MAP", "PrecipTempCorr", "isothermality",
                      "clay_surface", "awc"),  # placeholder
        squares = TRUE,
        interactions = TRUE,
        alpha = 1,
        cv = list(cluster_vars = c("MAT", "MAP", "PrecipTempCorr",
                                   'awc'),
                  k_clusters = 10,
                  select_rule = "1se")  # "1se" or "min"
      )
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
