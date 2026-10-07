# Specifications for the cover models.
#
# Nested: family -> response -> model version. The family decides which
# fitting script runs, the response names what is modelled, and the version
# (m01, m02, ...) is one set of decisions for that response.
#
#   classification  binomial lasso           02_fit_classification.R
#   cover           continuous cover         03_fit_cover.R
#                   (quasibinomial, or a ranger regression forest; beta
#                   later)
#   proportion      shares of a total (beta) 04_fit_proportion.R
#
# Every spec has:
#   response      column in the cover training data
#   rows          which pixel-years to fit to: "all", or a quoted condition on
#                 the training data in original units, e.g. quote(MAP < 700)
#                 (see select_rows())
#   engine        "glmnet" (penalized regression), "ranger" (random forest),
#                 "grpreg" (hierarchical group lasso: x^2 only enters with
#                 x, and x:z with x and z; see hier_groups()) or "gam" (one
#                 smooth per predictor; see gam_helpers.R)
#   pred_vars     predictors (short climate/soil names), with the log1p_
#                 versions of log1p_vars appended by the constructor
#   log1p_vars    predictors to add as log1p (e.g. MAP means log1p_MAP is used)
#   squares       add squared terms
#   interactions  add pairwise interactions
#   interact_log1p  let log1p_x interact even when x is also a predictor (never
#                 x:log1p_x); see make_design_matrix()
#   alpha         glmnet / grpreg penalty (1 = lasso); ignored by ranger
#   ranger        list of arguments passed on to ranger::ranger(); NULL
#                 otherwise
#   grpreg        list of arguments passed on to grpreg::grpreg() (e.g.
#                 nlambda); NULL otherwise
#   gam           list: k (knots per smooth), bs (basis), log10_mult
#                 (multipliers on REML's smoothing parameters to try, as
#                 log10), metric (CV score), interactions (moved here from
#                 the interactions argument: a smooth of x * z for every
#                 pair); NULL otherwise. For gam, interactions, squares and
#                 interact_log1p are then FALSE, so x has main effects only
#   cv            list: cluster_vars, n_folds, env_clusters, select_rule
#
# classification specs add (cover specs, from .defaults_cover(), don't have
# cover_threshold or threshold_method):
#   cover_threshold     % cover dividing the two classes
#   threshold_method    rule for turning predictions into a class
#   family              "binomial": fit to the class (cover above or below
#                       cover_threshold). "gaussian" or "poisson": fit to
#                       cover itself, then cut the predicted cover where the
#                       predicted fraction of forest matches the observed one
#                       (gaussian: glmnet, grpreg or ranger, a regression
#                       forest; poisson: grpreg only).
#                       "quasibinomial": fit to cover as a
#                       proportion with a binomial likelihood (logit link;
#                       glmnet only), predicting mean cover (0-100%), cut the
#                       same way as the continuous families
#   response_transform  applied to cover before fitting: "identity" or
#                       "log1p"; predictions are back-transformed to % cover
#
# Specs are built with a constructor per model (e.g. .defaults_class_forest()),
# so each version lists only what differs from the defaults.

source('Functions/models/predictors.R')
.default_cluster_vars <- c("MAT", "MAP", "PrecipTempCorr", "awc")

# complex model (for comparisons)
.pred_vars_complex1 <- c("MAT", "P_wettestMonth", "PrecipTempCorr", "isothermality",
                         "WDD_mean", "soilDepth", "clay_surface", "sand", "coarse",
                         "carbon")

.pred_vars_complex2	<-  c("MAT", "MAP", "P_driestMonth", "PrecipTempCorr", 
                          "isothermality",  "WD_mean", "VPD_max", "clay_surface", 
                          "sand", "coarse", "awc")

# zero tree model defaults
.thresh_zt <- 90 # zt stands for 'zero tree', % of pixel threshold 

#' Spec for the forest / non-forest classification
#'
#' Arguments override the defaults. Appends the `log1p_` names of
#' `log1p_vars` to `pred_vars`, and nests the CV settings under `cv`.
#'
#' @param response Column in the cover training data.
#' @param rows Which pixel-years to fit to: "all", or a quoted condition such
#'   as `quote(MAP < 700)`; see `select_rows()`.
#' @param engine "glmnet", "ranger", "grpreg" or "gam".
#' @param family "binomial" (fit to the class), "gaussian" / "poisson"
#'   (fit to cover; gaussian with glmnet, grpreg or ranger, poisson with
#'   grpreg only), or "quasibinomial" (fit to cover / 100 with
#'   a binomial likelihood; glmnet only).
#' @param response_transform "identity" or "log1p", applied to cover before
#'   fitting (gaussian / poisson only).
#' @param cover_threshold Percent cover dividing the two classes, as stored.
#' @param threshold_method See `?PresenceAbsence::optimal.thresholds`; only
#'   "PredPrev=Obs" for continuous families.
#' @param pred_vars Predictors, in original units.
#' @param log1p_vars Predictors to also add as `log1p_<var>`.
#' @param squares,interactions Add squared terms / pairwise interactions.
#' @param interact_log1p Let `log1p_x` interact even when `x` is also a
#'   predictor.
#' @param alpha glmnet / grpreg penalty mixing (1 = lasso).
#' @param ranger List of arguments for `ranger::ranger()`; required when
#'   `engine = "ranger"`.
#' @param grpreg List of extra arguments for `grpreg::grpreg()`. For grpreg
#'   specs, `max.iter` defaults to 1e5: grpreg's own default (1e4) caps the
#'   iterations for the whole lambda path, which can end a binomial path
#'   early.
#' @param cluster_vars Variables for the environmental clusters.
#' @param n_folds Number of CV folds.
#' @param env_clusters Number of environmental clusters (k-means), randomly
#'   grouped into the `n_folds` folds. NULL (the default) uses `n_folds`:
#'   one cluster per fold.
#' @param gam List of settings for `engine = "gam"`, overriding the
#'   defaults (see function body for defaults) `interactions = TRUE` adds a smooth of
#'   the product of each pair (see gam_formula()); it is stored as
#'   `gam$interactions`, and `interactions`, `squares` and `interact_log1p`
#'   are set to FALSE, so the design matrix holds the main effects only.
#' @param select_rule "min" or "1se" ("min" only for grpreg).
#' @return A spec list.
.defaults_class_forest <- function(response = "cov_tree",
                                   rows = "all",
                                   engine = "glmnet",
                                   family = "binomial",
                                   response_transform = "identity",
                                   cover_threshold = 10,
                                   threshold_method = "PredPrev=Obs",
                                   pred_vars = c("MAT", "MAP", "PrecipTempCorr",
                                                 "isothermality", "WD_p95",
                                                 "clay_surface", "awc"),
                                   log1p_vars = .possible_log_vars(pred_vars),
                                   squares = TRUE,
                                   interactions = TRUE,
                                   interact_log1p = FALSE,
                                   alpha = 1,
                                   ranger = NULL,
                                   grpreg = list(nlambda = 25),
                                   gam = NULL,
                                   cluster_vars = .default_cluster_vars,
                                   n_folds = 10,
                                   env_clusters = NULL,
                                   select_rule = "min") {
  stopifnot(identical(rows, "all") || is.language(rows),
            engine %in% c("glmnet", "ranger", "grpreg", "gam"),
            family %in% c("binomial", "quasibinomial", "gaussian", "poisson"),
            response_transform %in% c("identity", "log1p"),
            select_rule %in% c("min", "1se"),
            is.null(env_clusters) || env_clusters >= n_folds,
            engine != "ranger" || is.list(ranger),
            # ranger: a probability forest on the class, or (gaussian) a
            # regression forest on cover; poisson is grpreg only
            engine != "ranger" || family %in% c("binomial", "gaussian"),
            engine == "grpreg" || family != "poisson",
            # grpreg's binomial needs a 0/1 response, and ranger isn't set up
            # for proportions
            family != "quasibinomial" || engine %in% c("glmnet", "gam"),
            # gam smooths the predictors, so no log1p transforms
            engine != "gam" || length(log1p_vars) == 0,
            !family %in% c("binomial", "quasibinomial") ||
              response_transform == "identity",
            # predicted cover isn't a probability, so only prevalence matching
            family == "binomial" || threshold_method == "PredPrev=Obs",
            # grpreg models use cv.grpreg()'s lambda.min
            engine != "grpreg" || select_rule == "min")
  
  if (engine == "grpreg") {
    grpreg <- utils::modifyList(list(max.iter = 1e5), as.list(grpreg))
  }
  if (engine == "gam") {
    gam_default <- list(
      k = 5,
      bs = "cs",
      log10_mult = seq(-1, 4, by = 0.5),
      # as glmnet's CV: deviance, or mean squared error for gaussian
      metric = if (family == "gaussian") "mse" else "deviance_binomial",
      # smooths of products, built by gam_formula() rather than as columns
      # of the design matrix
      interactions = interactions
    )
    gam <- utils::modifyList(gam_default, as.list(gam))
    # x: main effects only (a smooth already covers x^2)
    interactions <- FALSE
    squares <- FALSE
    interact_log1p <- FALSE
  }
  
  list(
    response = response,
    rows = rows,
    engine = engine,
    family = family,
    response_transform = response_transform,
    cover_threshold = cover_threshold,
    threshold_method = threshold_method,
    # recycle0: with no log1p_vars, adds nothing rather than "log1p_"
    pred_vars = c(pred_vars, paste0("log1p_", log1p_vars, recycle0 = TRUE)),
    log1p_vars = log1p_vars,
    squares = squares,
    interactions = interactions,
    interact_log1p = interact_log1p,
    alpha = alpha,
    ranger = ranger,
    grpreg = grpreg,
    gam = gam,
    cv = list(cluster_vars = cluster_vars,
              n_folds = n_folds,
              env_clusters = if (is.null(env_clusters)) n_folds else env_clusters,
              select_rule = select_rule)
  )
}

.defaults_zero_tree = function(response = "pct_zero_tree",
                               cover_threshold = .thresh_zt, ...) {
  spec <- .defaults_class_forest(response = response, 
                                 cover_threshold = cover_threshold, ...)
  spec
}

#' Spec for a continuous cover model (03_fit_cover.R)
#'
#' As `.defaults_class_forest()`, which it calls (all its arguments can be
#' given), but quasibinomial by default (fit to cover / 100, predicting mean %
#' cover), and without the classification fields (`cover_threshold`,
#' `threshold_method`). For a ranger regression forest use `engine = "ranger"`
#' with `family = "gaussian"`.
#'
#' @param response Cover column in the training data, e.g. "cov_tree".
#' @param family "quasibinomial" or "gaussian" (not "binomial").
#' @param ... Other arguments to `.defaults_class_forest()`.
#' @return A spec list.
.defaults_cover <- function(response = "cov_tree", family = "quasibinomial",
                            ...) {
  stopifnot(family != "binomial")
  spec <- .defaults_class_forest(response = response, family = family, ...)
  spec$cover_threshold <- NULL
  spec$threshold_method <- NULL
  spec
}


cover_specs <- list(
  
  classification = list(
    
    # forest / non-forest: tree cover above or below 10%
    forest = list(
      # binomial GLM
      m01 = .defaults_class_forest(),
      m01.0 = .defaults_class_forest(interact_log1p = TRUE),
      # as m05, but fit to the class: the hierarchy alone, for comparison
      # with m01 (same response, no hierarchy) and m05 (same hierarchy,
      # continuous response)
      m01.1 = .defaults_class_forest(engine = "grpreg",
                                     interact_log1p = TRUE),
      # complex covarariate comparison
      m01.2 = .defaults_class_forest(
        pred_vars = .pred_vars_complex2,
        log1p_vars = .possible_log_vars(.pred_vars_complex2),
        interact_log1p = TRUE
        ),
      # less climate extrapolation
      m01.3 = .defaults_class_forest(n_folds = 10, env_clusters = 20),

      # elastic net
      m02 = .defaults_class_forest(alpha = 0.5),
      
      # less climate extrapolation
      m03 = .defaults_class_forest(n_folds = 50),
      
      # random forest, as a benchmark for how well the same predictors can do
      # without the constraint of an equation.  Not a candidate for prediction.
      # min.node.size is the smallest node that can be split, not the smallest
      # leaf. min.bucket (smallest leaf) was tried and roughly tripled fit time
      m04 = .defaults_class_forest(engine = "ranger",
                                   log1p_vars = character(0),
                                   squares = FALSE,
                                   interactions = FALSE,
                                   ranger = list(num.trees = 300,
                                                 min.node.size = 100,
                                                 importance = "permutation")),
      m04.1 = .defaults_class_forest(engine = "ranger",
                                   pred_vars = .pred_vars_complex2,
                                   log1p_vars = character(0),
                                   squares = FALSE,
                                   interactions = FALSE,
                                   ranger = list(num.trees = 300,
                                                 min.node.size = 100,
                                                 importance = "permutation")),
      
      # hierarchical group lasso on continuous tree cover. x^2 only enters
      # with x, and an interaction only with both of its terms. Squares don't
      # interact. log1p_x can enter without x and interacts on its own
      # (log1p_MAP:MAT needs log1p_MAP, not MAP).
      # Fit to log1p(cover): handles zeros, and spreads out the low covers
      # around the 10% threshold while compressing high ones. Predicted cover
      # is then cut where the predicted fraction of forest matches the
      # observed one, as the probabilities are for m01
      m05 = .defaults_class_forest(engine = "grpreg",
                                   family = "gaussian",
                                   response_transform = "log1p",
                                   interact_log1p = TRUE),
      
      # as m01 (glmnet lasso, same predictors), but fit to cover / 100 with a
      # binomial likelihood (what quasibinomial estimates; the dispersion
      # doesn't matter for the lasso). Predicts mean cover, which stays in
      # 0-100% unlike m05, and is cut where the predicted fraction of forest
      # matches the observed one, as m05 is
      m06 = .defaults_class_forest(family = "quasibinomial"),
      
      # models fit to a filtered subset of the data
      m07.0 = .defaults_class_forest(rows = quote(MAP < 700)),
      m07.1 = .defaults_class_forest(rows = quote(MAP >= 700)),
      m07.2 =  .defaults_class_forest(engine = "ranger",
                                     rows = quote(MAP < 700),
                                    log1p_vars = character(0),
                                    squares = FALSE,
                                    interactions = FALSE,
                                    ranger = list(num.trees = 300,
                                                  min.node.size = 100)),
      m07.3 =  .defaults_class_forest(engine = "ranger",
                                     rows = quote(MAP >= 700),
                                    log1p_vars = character(0),
                                    squares = FALSE,
                                    interactions = FALSE,
                                    ranger = list(num.trees = 300,
                                                  min.node.size = 100)),
      m07.4 = .defaults_class_forest(rows = quote(MAP < 800)),
      m07.5 = .defaults_class_forest(rows = quote(MAP >= 800)),
      m07.6 =  .defaults_class_forest(engine = "ranger",
                                      rows = quote(MAP < 800),
                                      log1p_vars = character(0),
                                      squares = FALSE,
                                      interactions = FALSE,
                                      ranger = list(num.trees = 300,
                                                    min.node.size = 100)),
      m07.7 =  .defaults_class_forest(engine = "ranger",
                                      rows = quote(MAP >= 800),
                                      log1p_vars = character(0),
                                      squares = FALSE,
                                      interactions = FALSE,
                                      ranger = list(num.trees = 300,
                                                    min.node.size = 100)),
      # same idea as 7.0 and 7.1, but with continuous model
      m07.8 = .defaults_class_forest(rows = quote(MAP < 700),
                                   family = "gaussian",
                                   response_transform = "log1p"),
      m07.9 = .defaults_class_forest(rows = quote(MAP >= 700),
                                     family = "gaussian",
                                     response_transform = "log1p"),
      
      # GAM: one smooth per predictor in place of log1p and squares; main
      # effects only here (interactions = TRUE adds a smooth of each pair's
      # product). REML sets each term's smoothing, and CV on the
      # environmental folds scales it all up or down (see gam_helpers.R).
      # Same predictors as m01, fit to the class
      m08.0 = .defaults_class_forest(engine = "gam",
                                   log1p_vars = character(0),
                                   gam = list(
                                     log10_mult = seq(-1, 4, by = 0.5)
                                     ),
                                   interactions = TRUE),
      # as m08, but fit to cover / 100, as m06
      m08.1 = .defaults_class_forest(engine = "gam",
                                     family = "quasibinomial",
                                     log1p_vars = character(0),
                                     gam = list(
                                       log10_mult = seq(-1, 4, by = 0.5)
                                       ),
                                     interactions = TRUE),
      # 8.0 but complex set of predictors
      m08.2 = .defaults_class_forest(engine = "gam",
                                     pred_vars = .pred_vars_complex1,
                                     log1p_vars = character(0),
                                     gam = list(
                                       log10_mult = seq(-1, 4, by = 0.5)
                                     ),
                                     interactions = TRUE),
      # 8.0 but w/ less climat extrapolation
      m08.3 = .defaults_class_forest(engine = "gam",
                                     log1p_vars = character(0),
                                     gam = list(
                                       log10_mult = seq(-1, 4, by = 0.5)
                                     ),
                                     interactions = TRUE,
                                     n_folds = 10,
                                     env_clusters = 20)
    ),
    
    # zero tree: % of a cell's natural land with < 3% RAP tree cover; above
    # 90% = zero tree. Data: c02 only (05_rap_training_sample.R)
    zero_tree = list(
      # binomial GLM
      m01.0 = .defaults_zero_tree(interact_log1p = TRUE),
      m01.2 = .defaults_zero_tree(
        pred_vars = .pred_vars_complex2,
        log1p_vars = .possible_log_vars(.pred_vars_complex2),
        interact_log1p = TRUE
      ),
      # random forest
      m04.0 = .defaults_zero_tree(
        engine = "ranger",
       log1p_vars = character(0),
       squares = FALSE,
       interactions = FALSE,
       ranger = list(num.trees = 300,
                     min.node.size = 100,
                     importance = "permutation")),
      m04.1 = .defaults_zero_tree(
        engine = "ranger",
         pred_vars = .pred_vars_complex2,
         log1p_vars = character(0),
         squares = FALSE,
         interactions = FALSE,
         ranger = list(num.trees = 300,
                       min.node.size = 100,
                       importance = "permutation"))
    )

  ),
  
    # continuous cover. Tree cover is fit separately above and at or below 10%
  # (the forest threshold), so neither is pulled toward the other's mean; the
  # two are blended at 10%, and zero-tree areas set to 0 afterwards. Zeros are
  # kept. In c02, shrub and herbaceous cover are NA wherever tree cover is
  # >= 10% (RAP understorey), so those models are non-forest only there.
  # m01: glmnet lasso, quasibinomial. m02: ranger regression forest, as a
  # benchmark (as forest m04)
  cover = list(
    tree_forest = list(
      m01.0 = .defaults_cover(rows = quote(cov_tree > 10), interact_log1p = TRUE),
      m04.0 = .defaults_cover(rows = quote(cov_tree > 10),
                            engine = "ranger", family = "gaussian",
                            log1p_vars = character(0), squares = FALSE,
                            interactions = FALSE,
                            ranger = list(num.trees = 300,
                                          min.node.size = 100))
    ),
    tree_nonforest = list(
      m01.0 = .defaults_cover(rows = quote(cov_tree <= 10), interact_log1p = TRUE),
      m04.0 = .defaults_cover(rows = quote(cov_tree <= 10),
                            engine = "ranger", family = "gaussian",
                            log1p_vars = character(0), squares = FALSE,
                            interactions = FALSE,
                            ranger = list(num.trees = 300,
                                          min.node.size = 100))
    )
  ),
  # shares: needleleaf (of tree), forb / C3 / C4 (of herbaceous, rescaled to
  # sum to 1 at prediction time). Not specified yet.
  proportion = list()
)
