# Building model predictors from the cover training data (or from rasters
# converted to a data frame).
#
# The order is always:
#   1. log1p versions, made from the original units (names starting "log1p_")
#   2. standardize with the fixed global parameters (read_scale_params())
#   3. squares and pairwise interactions, via the design matrix
# so the coefficients mean the same thing wherever the model is applied.
#
# Each @examples block assigns every argument, so can run it and then step
# through the function body line by line.


#' Design matrix of main effects, squares and pairwise interactions
#'
#' Builds the formula from the pieces requested, then `model.matrix()`. The
#' intercept column is dropped, because glmnet fits its own intercept.
#'
#' `log1p_x` is never squared: log1p is already its transformation. If both
#' `x` and `log1p_x` are predictors, `log1p_x` is a main effect only and
#' interactions use `x`. If only `log1p_x` is present, it stands in for `x`
#' in interactions.
#'
#' @param dat Data frame containing `pred_vars`, with no missing values in
#'   them (model.matrix would silently drop those rows).
#' @param pred_vars Predictor columns.
#' @param squares Logical; add squared terms.
#' @param interactions Logical; add pairwise interactions.
#' @return Numeric matrix, one row per row of `dat`.
#' @examples
#' dat <- mtcars
#' dat$log1p_hp <- log1p(dat$hp)
#' pred_vars <- c("wt", "hp", "log1p_hp")
#' squares <- TRUE
#' interactions <- TRUE
#' x <- make_design_matrix(dat = dat, pred_vars = pred_vars,
#'                         squares = squares, interactions = interactions)
#' colnames(x)  # log1p_hp: main effect only
make_design_matrix <- function(dat, pred_vars, squares, interactions) {
  stopifnot(all(pred_vars %in% names(dat)),
            !anyNA(dat[pred_vars]))
  
  terms_main <- pred_vars
  
  is_log1p <- str_detect(pred_vars, "^log1p_")
  
  # squares: never for log1p_ variables
  square_vars <- pred_vars[!is_log1p]
  
  # interactions: log1p_x only if x itself isn't a predictor
  source_names <- str_remove(pred_vars, "^log1p_")
  inter_vars <- pred_vars[!(is_log1p & source_names %in% pred_vars)]
  
  terms_squares <- if (squares) paste0("I(", square_vars, "^2)") else NULL
  
  terms_inter <- if (interactions && length(inter_vars) >= 2) {
    pairs <- utils::combn(inter_vars, 2)  # 2 rows, one column per pair
    paste(pairs[1, ], pairs[2, ], sep = ":")
  } else {
    NULL
  }
  
  form <- stats::reformulate(c(terms_main, terms_squares, terms_inter))
  
  x <- stats::model.matrix(form, data = dat)
  x[, colnames(x) != "(Intercept)", drop = FALSE]
}


#' Transform, standardize and expand predictors into a design matrix
#'
#' Predictors named "log1p_<var>" are made from `<var>` in its original units.
#' All predictors are then standardized with the fixed global parameters and
#' expanded into squares and interactions.
#'
#' @param dat Data frame in original units (e.g.
#'   `read_cover_training(vc, normalize = FALSE)`), with no missing values in
#'   the predictors' source columns.
#' @param pred_vars Predictors, e.g. `c("MAT", "log1p_MAP")`.
#' @param scale_df Scaling parameters (`read_scale_params()`), which must
#'   include every predictor, "log1p_" ones included.
#' @param squares Logical; add squared terms.
#' @param interactions Logical; add pairwise interactions.
#' @return Numeric matrix, one row per row of `dat`.
#' @examples
#' dat <- mtcars
#' pred_vars <- c("wt", "log1p_hp")
#' scale_df <- data.frame(variable = c("wt", "log1p_hp"),
#'                        mean = c(mean(mtcars$wt), mean(log1p(mtcars$hp))),
#'                        sd = c(sd(mtcars$wt), sd(log1p(mtcars$hp))))
#' squares <- TRUE
#' interactions <- TRUE
#' x <- prepare_predictors(dat = dat, pred_vars = pred_vars,
#'                         scale_df = scale_df, squares = squares,
#'                         interactions = interactions)
#' head(x)
prepare_predictors <- function(dat, pred_vars, scale_df, squares,
                               interactions) {
  stopifnot(all(pred_vars %in% scale_df$variable))
  
  # 1. log1p versions from the original units
  log1p_vars <- str_subset(pred_vars, "^log1p_")
  source_vars <- str_remove(log1p_vars, "^log1p_")
  stopifnot(all(source_vars %in% names(dat)))
  dat[log1p_vars] <- map(dat[source_vars], log1p)
  
  # 2. standardize with the fixed global parameters
  dat <- standardize(dat, vars = pred_vars, scale_df = scale_df)$data
  
  # 3. squares and interactions
  make_design_matrix(dat, pred_vars = pred_vars, squares = squares,
                     interactions = interactions)
}
