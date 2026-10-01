# Helpers for the "gam" engine: one smooth per predictor, using a cubic
# regression spline basis ("cs"), which can be written out as knots and knot
# values.
#
# Smoothing is set at two levels. REML picks a smoothing parameter (sp) for
# each term, which decides which terms are wiggly and which shrink out
# ("cs" lets a term shrink to zero). Cross-validation then picks one
# multiplier applied to all of them, which sets the overall amount of
# smoothing (> 1 is smoother than REML). REML is rerun on each fold's
# training rows, so held-out rows never inform it.
#
# Not specific to classification: takes a data frame, a family object and
# fold ids, so the cover and proportion models can use it too.


#' Formula with one smooth per predictor
#'
#' @param response Name of the response column.
#' @param pred_vars Predictor columns.
#' @param k Basis dimension (number of knots) of each smooth.
#' @param bs Basis: "cs" (cubic regression spline that can shrink to zero)
#'   or "cr" (same basis, can only shrink to a straight line).
#' @param interactions Logical; also add a smooth of the product of each pair
#'   of predictors, e.g. s(I(MAT * MAP)). mgcv computes the product from the
#'   predictor columns, when fitting and when predicting.
#' @return A formula.
#' @examples
#' response <- ".y"
#' pred_vars <- c("MAT", "MAP", "awc")
#' k <- 5
#' bs <- "cs"
#' interactions <- TRUE
#' gam_formula(response = response, pred_vars = pred_vars, k = k, bs = bs,
#'             interactions = interactions)
gam_formula <- function(response, pred_vars, k, bs = "cs",
                        interactions = FALSE) {
  stopifnot(bs %in% c("cs", "cr"))
  vars <- pred_vars
  if (interactions) {
    pairs <- utils::combn(pred_vars, 2)  # 2 rows, one column per pair
    vars <- c(vars, paste0("I(", pairs[1, ], " * ", pairs[2, ], ")"))
  }
  terms <- paste0("s(", vars, ", k = ", k, ", bs = \"", bs, "\")")
  stats::reformulate(terms, response = response)
}


#' Fit a GAM, by REML or with REML's smoothing parameters scaled
#'
#' `bam()` with discretized covariates, which is much faster on large data.
#' With `sp_reml = NULL` the smoothing parameters are estimated by REML (read
#' them from `$sp` of the result); otherwise they are fixed at
#' `sp_reml * mult`.
#'
#' @param formula From `gam_formula()`.
#' @param data Data frame with the response and predictor columns.
#' @param family Family object, e.g. `binomial()`.
#' @param sp_reml Per-term smoothing parameters from a REML fit, or NULL.
#' @param mult Multiplier on `sp_reml`; ignored when `sp_reml` is NULL.
#' @return A `bam` object.
#' @examples
#' set.seed(1)
#' data <- data.frame(a = rnorm(2000), b = rnorm(2000))
#' data$.y <- rbinom(2000, 1, plogis(sin(data$a)))
#' formula <- gam_formula(".y", c("a", "b"), k = 5)
#' family <- binomial()
#' m0 <- fit_gam(formula, data, family)
#' sp_reml <- m0$sp
#' mult <- 15
#' m <- fit_gam(formula = formula, data = data, family = family,
#'              sp_reml = sp_reml, mult = mult)
#'  par(mfrow = c(2, 2))
#' plot(m0)
#' plot(m)
fit_gam <- function(formula, data, family, sp_reml = NULL, mult = 1) {
  sp <- if (is.null(sp_reml)) NULL else sp_reml * mult
  mgcv::bam(formula, data = data, family = family, method = "fREML",
            discrete = TRUE, sp = sp)
}


#' Choose the smoothing multiplier by cross-validation
#'
#' For each fold, REML on the training rows gives the per-term smoothing
#' parameters; the model is then refit at each multiplier and predicted on the
#' held-out rows. Knots are placed on the training rows, so held-out values
#' beyond them are extrapolated (linearly), as future climates would be.
#'
#' Fold scores are averaged weighted by fold size, as `cv.glmnet()` does, and
#' the multiplier chosen with `select_lambda()`. Those functions expect the
#' tuning value in a column named `lambda`, so the multiplier is stored there.
#'
#' @param formula,data,family As in `fit_gam()`.
#' @param foldid Fold of each row of `data`.
#' @param mults Multipliers to try.
#' @param metric Name of a metric in `get_metric_fun()`, computed on the
#'   scale the model is fit on (e.g. probability, proportion, log1p cover).
#' @param rule Selection rule passed to `select_lambda()`: "min" or "1se".
#' @return List: `mult` (the selected multiplier), `summary` (one row per
#'   multiplier: weighted mean `score` and its `score_se`), `scores` (one row
#'   per fold and multiplier) and `oof` (out-of-fold predictions at the
#'   selected multiplier, on the response scale).
#' @examples
#' set.seed(1)
#' data <- data.frame(a = rnorm(5000), b = rnorm(5000))
#' data$.y <- rbinom(5000, 1, plogis(sin(data$a)))
#' formula <- gam_formula(".y", c("a", "b"), k = 5)
#' family <- binomial()
#' foldid <- kmeans(data[c("a", "b")], centers = 5)$cluster
#' mults <- 10^seq(-3, 3, by = 0.5)
#' metric <- "deviance_binomial"
#' rule <- "min"
#' cv <- cv_gam_mult(formula = formula, data = data, family = family,
#'                   foldid = foldid, mults = mults, metric = metric,
#'                   rule = rule)
#' cv$summary
cv_gam_mult <- function(formula, data, family, foldid, mults, metric, rule) {
  stopifnot(length(foldid) == nrow(data), !anyNA(foldid))
  y <- data[[all.vars(formula)[1]]]
  metric_fun <- get_metric_fun(metric)
  
  # out-of-fold predictions: one column per multiplier
  oof <- matrix(NA_real_, nrow = nrow(data), ncol = length(mults))
  for (f in unique(foldid)) {
    test <- foldid == f
    sp_reml <- fit_gam(formula, data[!test, ], family)$sp
    for (j in seq_along(mults)) {
      m <- fit_gam(formula, data[!test, ], family, sp_reml = sp_reml, 
                   mult = mults[j])
      oof[test, j] <- predict(m, data[test, ], type = "response",
                              discrete = FALSE)
    }
  }
  
  scores <- purrr::map(unique(foldid), \(f) {
    test <- foldid == f
    tibble::tibble(
      fold_id = f,
      n_test = sum(test),
      lambda = mults,
      score = apply(oof[test, , drop = FALSE], 2,
                    \(p) metric_fun(y[test], p))
    )
  }) |>
    dplyr::bind_rows()
  
  summary <- summarize_scores(scores, "score", weight_col = "n_test")
  mult <- select_lambda(summary, metric = "score", rule = rule)$lambda
  if (mult %in% range(mults)) {
    warning("selected multiplier (", signif(mult, 3), ") is at the edge ",
            "of the grid; the best value may lie outside it")
  }
  
  list(mult = mult,
       summary = summary,
       scores = scores,
       oof = oof[, match(mult, mults)])
}


#' Drop the per-row parts of a fitted GAM
#'
#' Cuts the size of a `bam` object by ~95% for saving; `predict()` on new data
#' is unaffected.
#'
#' @param m A `bam` object.
#' @return `m` without its model frame, response, fitted values, residuals
#'   and weights.
strip_gam <- function(m) {
  m[c("model", "y", "fitted.values", "linear.predictors", "residuals",
      "prior.weights", "weights", "offset", "wt")] <- NULL
  m
}


#' Effective degrees of freedom of each smooth
#'
#' About 1 is close to a straight line; about 0 means the term has dropped
#' out. With k knots the maximum is k - 1.
#'
#' @param m A fitted `gam` or `bam` object.
#' @return Tibble with `term` and `edf`.
gam_edf <- function(m) {
  s_table <- summary(m)$s.table
  tibble::tibble(term = rownames(s_table), edf = s_table[, "edf"])
}
