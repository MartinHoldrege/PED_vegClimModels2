# main functions for cross validation

#' Run path fitting, scoring, and lambda selection for one fold
#'
#' This function runs the full train/test workflow for one fold:
#' it computes `lambda_max_l1` on training data, builds a lambda path,
#' fits the path on training data, scores the path on test data, and
#' selects a best lambda.
#'
#' @param data Data frame used for modeling.
#' @param fold Fold object containing row indices, such as from
#'   `make_cluster_folds()`.
#' @param en_alpha Elastic-net mixing parameter.
#' @param lambda_max_fun Function that computes `lambda_max_l1` on training data.
#' @param lambda_max_args Named list of additional arguments passed to
#'   `lambda_max_fun`.
#' @param lambda_path_fun Function that creates a lambda path.
#' @param lambda_path_args Named list of additional arguments passed to
#'   `lambda_path_fun`.
#' @param fit_path_fun Function that fits models across a lambda path.
#' @param fit_path_args Named list of additional arguments passed to
#'   `fit_path_fun`.
#' @param score_fun Function that scores a fitted path on test data.
#' @param score_args Named list of additional arguments passed to `score_fun`.
#' @param select_fun Function that selects the best lambda from scored results.
#' @param select_args Named list of additional arguments passed to `select_fun`.
#' @param keep_best_fit Logical. If `TRUE`, return the full fit object for the
#'   selected lambda.
#' @param keep_path_fit Logical. If `TRUE`, return the full path-fit object.
#'
#' @return A list containing fold metadata, scored results, selected lambda,
#'   lightweight parameter objects for all lambdas, and optionally the best
#'   full fit and/or the full path-fit object.
#' @examples
#'
#' # simulate small dataset
#' set.seed(1)
#' dummy <- cwexp_make_dummy_data(n = 1000)
#' dat <- dummy$data
#' cover_cols <- dummy$spec$cover_cols
#'
#' # environmental clusters
#' env_obj <- make_env_clusters(
#'   data = dat,
#'   vars = c("tmean", "ppt", "vpd", "sand"),
#'   k = 6,
#'   seed = 1
#' )
#'
#' # create cluster-based folds
#' folds <- make_cluster_folds(
#'   env_cluster = env_obj$env_cluster,
#'   n_folds = 2,
#'   seed = 1
#' )
#'
#' # compile model
#' dll <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
#'
#' # run workflow for one fold
#' res <- run_fold_path(
#'   data = dat,
#'   fold = folds[[1]],
#'   en_alpha = 0.5,
#'
#'   lambda_max_fun = cwexp_lambda_max_l1_tmb,
#'   lambda_max_args = list(
#'     formula = totalBio ~ tmean + ppt + vpd + sand,
#'     cover_cols = cover_cols,
#'     dll_en = dll
#'   ),
#'
#'   lambda_path_args = list(
#'     n_lambda = 5,
#'     include_zero = TRUE
#'   ),
#'
#'   fit_path_fun = cwexp_fit_lambda_path_tmb,
#'   fit_path_args = list(
#'     formula = totalBio ~ tmean + ppt + vpd + sand,
#'     cover_cols = cover_cols,
#'     dll = dll,
#'     include_report = FALSE
#'   )
#' )
#'
#' # inspect scored lambda path
#' res$scores
#'
#' # selected lambda
#' res$selected
#'
#' # parameters for all lambdas
#' res$path_par
#'
#' # best fitted model
#' res$best_fit
run_fold_path <- function(data,
                          fold,
                          en_alpha,
                          lambda_max_fun,
                          lambda_max_args = list(),
                          lambda_path_fun = cwexp_lambda_path,
                          lambda_path_args = list(),
                          fit_path_fun,
                          fit_path_args = list(),
                          score_fun = score_path,
                          score_args = list(),
                          select_fun = select_lambda,
                          select_args = list(
                            metric = "mae_log",
                            rule = "min",
                            tol = 0
                          ),
                          keep_best_fit = TRUE,
                          keep_path_fit = FALSE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.list(fold) || is.null(fold$train_rows) || is.null(fold$test_rows)) {
    stop("fold must be a list with elements 'train_rows' and 'test_rows'.")
  }
  if (!is.function(lambda_max_fun)) {
    stop("lambda_max_fun must be a function.")
  }
  if (!is.function(lambda_path_fun)) {
    stop("lambda_path_fun must be a function.")
  }
  if (!is.function(fit_path_fun)) {
    stop("fit_path_fun must be a function.")
  }
  if (!is.function(score_fun)) {
    stop("score_fun must be a function.")
  }
  if (!is.function(select_fun)) {
    stop("select_fun must be a function.")
  }
  if (!is.finite(en_alpha) || en_alpha <= 0 || en_alpha > 1) {
    stop("en_alpha must be in (0, 1].")
  }
  
  split_dat <- subset_fold_data(data = data, fold = fold)
  train_dat <- split_dat$train
  test_dat <- split_dat$test
  
  lambda_max_obj <- do.call(
    lambda_max_fun,
    c(list(data = train_dat), lambda_max_args)
  )
  
  if (is.null(lambda_max_obj$lambda_max_l1)) {
    stop("lambda_max_fun output must contain 'lambda_max_l1'.")
  }
  
  lambda_seq <- do.call(
    lambda_path_fun,
    c(
      list(
        lambda_max_l1 = lambda_max_obj$lambda_max_l1,
        en_alpha = en_alpha
      ),
      lambda_path_args
    )
  )
  
  path_fit <- do.call(
    fit_path_fun,
    c(
      list(
        data = train_dat,
        en_alpha = en_alpha,
        lambda_seq = lambda_seq
      ),
      fit_path_args
    )
  )
  
  score_df <- do.call(
    score_fun,
    c(
      list(
        path_fit = path_fit,
        newdata = test_dat
      ),
      score_args
    )
  )
  
  selected <- do.call(
    select_fun,
    c(
      list(score_df = score_df),
      select_args
    )
  )
  
  # lightweight parameter objects for every lambda
  path_par <- extract_path_par_objects(path_fit)
  
  # identify selected model
  best_idx <- match(selected$lambda[[1]], path_fit$summary$lambda)
  
  spec = path_fit$fits[[1]]$spec
  spec$lambda <- NA # generic specification info
  out <- list(
    fold_id = fold$fold_id,
    train_rows = fold$train_rows,
    test_rows = fold$test_rows,
    train_clusters = fold$train_clusters,
    test_clusters = fold$test_clusters,
    n_train = nrow(train_dat),
    n_test = nrow(test_dat),
    en_alpha = en_alpha,
    lambda_max_l1 = lambda_max_obj$lambda_max_l1,
    lambda_seq = lambda_seq,
    scores = score_df,
    selected = selected,
    path_par = path_par,
    spec = spec
  )
  
  if (keep_best_fit) {
    out$best_fit <- path_fit$fits[[best_idx]]
  }
  
  if (keep_path_fit) {
    out$path_fit <- path_fit
  }
  
  class(out) <- "fold_path_result"
  out
}
