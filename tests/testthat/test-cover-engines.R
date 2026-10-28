# Tests that the cover model engines (glmnet, grpreg, grplasso), run through
# fit_engine() as the cover pipeline runs them, fit the model they should.
#
# Simulated 0/1 data with a large n and strong effects, so the penalty CV
# picks is small and every engine should land close to an unpenalized glm().
# The engines pick their penalties separately, so close, not identical.

library(glmnet)  # fit_engine() calls cv.glmnet() without glmnet::
library(purrr)
library(stringr)
library(tibble)
library(tidyr)
root <- here::here() # needed b/ testthat often treats wd as tests/testthat
source(file.path(root, "Functions/data/prepare_data.R"))
source(file.path(root, "Functions/models/fit_engine.R"))
source(file.path(root, "Functions/models/classification_helpers.R"))
# cover_specs.R sources predictors.R by a path relative to the project root
withr::with_dir(root, source("Functions/models/cover_specs.R"))

engines <- c("glmnet", "grpreg", "grplasso")

# simulated predictors already have mean 0 and sd 1, so the scaling in
# prepare_predictors() leaves them unchanged and glm() can use the same matrix
set.seed(1)
n <- 1000
dat <- tibble(x = rnorm(n), z = rnorm(n), w = rnorm(n))
scale_df <- tibble(variable = c("x", "z", "w"), mean = 0, sd = 1)
foldid <- sample(rep(1:10, length.out = n))

#' Fit one engine through fit_engine(), and return its final coefficients
#'
#' @param engine "glmnet", "grpreg" or "grplasso".
#' @param pred_vars Predictors (columns of dat).
#' @param interactions Logical; add pairwise interactions.
#' @param y 0/1 response.
#' @return List: `coef` (one per design-matrix column, intercept first, as in
#'   the final equation), `x` (the design matrix), `pred_in` (fit_engine()'s
#'   in-sample predictions, from predict_score()).
#' @examples
#' engine <- "grplasso"
#' pred_vars <- c("x", "z")
#' interactions <- TRUE
#' y <- rbinom(n, 1, plogis(dat$x))
#' fit_one(engine = engine, pred_vars = pred_vars,
#'         interactions = interactions, y = y)
fit_one <- function(engine, pred_vars, interactions, y) {
  spec <- .defaults_class_forest(response = "y", engine = engine,
                                 pred_vars = pred_vars,
                                 log1p_vars = character(0), squares = FALSE,
                                 interactions = interactions)
  x <- prepare_predictors(dat, pred_vars = pred_vars, scale_df = scale_df,
                          squares = FALSE, interactions = interactions)

  # with effects this strong, CV can pick the smallest penalty on the path,
  # which fit_engine() warns about; expected here
  res <- suppressWarnings(fit_engine(spec, x = x, y_fit = y, foldid = foldid))

  b <- if (engine == "glmnet") {
    b <- coef(res$fit, s = res$lambda)
    setNames(as.numeric(b), rownames(b))
  } else {
    res$fit$coef  # already summed over each column's repeats in the groups
  }
  list(coef = b[c("(Intercept)", colnames(x))], x = x, pred_in = res$pred_in)
}

test_that("main effects only: every engine is close to glm()", {
  # w has no effect
  y <- rbinom(n, 1, plogis(-0.5 + 1 * dat$x - 0.5 * dat$z))
  x <- prepare_predictors(dat, pred_vars = c("x", "z", "w"),
                          scale_df = scale_df, squares = FALSE,
                          interactions = FALSE)
  glm_coef <- unname(coef(glm(y ~ x, family = binomial())))

  for (engine in engines) {
    fit <- fit_one(engine, c("x", "z", "w"), interactions = FALSE, y = y)
    expect_equal(unname(fit$coef), glm_coef, tolerance = 0.15,
                 label = engine)
  }
})

test_that("x, z and x:z all real: every engine is close to glm()", {
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * dat$x - 0.6 * dat$z +
                             0.5 * dat$x * dat$z))
  x <- prepare_predictors(dat, pred_vars = c("x", "z"), scale_df = scale_df,
                          squares = FALSE, interactions = TRUE)
  glm_coef <- unname(coef(glm(y ~ x, family = binomial())))

  for (engine in engines) {
    fit <- fit_one(engine, c("x", "z"), interactions = TRUE, y = y)
    expect_equal(unname(fit$coef), glm_coef, tolerance = 0.15,
                 label = engine)
  }
})

test_that("hierarchical engines keep x and z whenever x:z is in", {
  # only the interaction has an effect, so without a hierarchy x and z could
  # be left out
  y <- rbinom(n, 1, plogis(-0.5 + 0.5 * dat$x * dat$z))

  for (engine in c("grpreg", "grplasso")) {
    b <- fit_one(engine, c("x", "z"), interactions = TRUE, y = y)$coef
    expect_true(b[["x:z"]] != 0, label = engine)
    expect_true(b[["x"]] != 0 && b[["z"]] != 0, label = engine)
  }
})

test_that("predictions are made from the reported coefficients", {
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * dat$x - 0.6 * dat$z +
                             0.5 * dat$x * dat$z))

  for (engine in engines) {
    fit <- fit_one(engine, c("x", "z"), interactions = TRUE, y = y)
    by_hand <- plogis(fit$coef[[1]] + drop(fit$x %*% fit$coef[-1]))
    # unname(): the design matrix's row names carry over to by_hand
    expect_equal(unname(fit$pred_in), unname(by_hand), label = engine)
  }
})
