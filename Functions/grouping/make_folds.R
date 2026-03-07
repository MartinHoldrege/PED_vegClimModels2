
#' Make cross-validation folds from environmental clusters
#'
#' Creates train/test folds by assigning whole clusters to test folds.
#' Each cluster appears in the test set exactly once.Clusters can
#' be created w/ make_env_clusters()
#'
#' @param env_cluster Integer or character vector of cluster assignments, one
#'   value per row of the data.
#' @param n_folds Integer number of folds.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list of folds. Each fold is a list with elements:
#' \describe{
#'   \item{fold_id}{Fold number.}
#'   \item{test_clusters}{Vector of cluster IDs in the test set.}
#'   \item{train_clusters}{Vector of cluster IDs in the training set.}
#'   \item{test_rows}{Integer row indices in the test set.}
#'   \item{train_rows}{Integer row indices in the training set.}
#' }
#'
#' @examples
#' env_cluster <- sample(1:20, size = 100, replace = TRUE)
#'
#' folds <- make_cluster_folds(env_cluster = env_cluster, n_folds = 3, seed = 1)
#'
#' folds[[1]]
make_cluster_folds <- function(env_cluster,
                               n_folds,
                               seed = NULL) {
  if (length(env_cluster) < 1) {
    stop("env_cluster must have length >= 1.")
  }
  if (!is.numeric(n_folds) || length(n_folds) != 1 ||
      n_folds < 2 || n_folds != as.integer(n_folds)) {
    stop("n_folds must be an integer >= 2.")
  }
  if (anyNA(env_cluster)) {
    stop("env_cluster must not contain missing values.")
  }
  
  cluster_ids <- unique(env_cluster)
  n_cluster <- length(cluster_ids)
  
  if (n_folds > n_cluster) {
    stop("n_folds cannot be greater than the number of unique clusters.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Randomize cluster order
  cluster_ids <- sample(cluster_ids, size = n_cluster, replace = FALSE)
  
  # Assign clusters to folds, split as evenly as possible across folds
  fold_id_by_cluster <- rep(seq_len(n_folds), length.out = n_cluster)
  
  folds <- vector("list", n_folds)
  
  for (i in seq_len(n_folds)) {
    test_clusters <- cluster_ids[fold_id_by_cluster == i]
    train_clusters <- setdiff(cluster_ids, test_clusters)
    
    test_rows <- which(env_cluster %in% test_clusters)
    train_rows <- which(env_cluster %in% train_clusters)
    
    folds[[i]] <- list(
      fold_id = i,
      test_clusters = test_clusters,
      train_clusters = train_clusters,
      test_rows = test_rows,
      train_rows = train_rows,
      n_test = length(test_rows),
      n_train = length(train_rows)
    )
  }
  
  class(folds) <- "cluster_folds"
  folds
}

#' Subset data into training and test sets for one fold
#'
#' @param data Data frame used for modeling.
#' @param fold Fold object containing `train_rows` and `test_rows`.
#'
#' @return A list with elements `train` and `test`.
subset_fold_data <- function(data, fold) {
   stopifnot(is.data.frame(data),
             is.list(fold),
             !is.null(fold$train_rows),
            !is.null(fold$test_rows),
             length(fold$train_rows) + length(fold$test_rows) == nrow(data))
  
  list(
    train = data[fold$train_rows, , drop = FALSE],
    test = data[fold$test_rows, , drop = FALSE]
  )
}


#' Check a set of cluster-based folds
#'
#' Performs basic validation checks on folds created from whole-cluster
#' assignments.
#'
#' @param folds A fold object such as returned by `make_cluster_folds()`.
#' @param env_cluster Vector of cluster assignments, one per row of the data.
#'
#' @return Invisibly returns `TRUE` if all checks pass. Otherwise throws an
#'   error describing the problem.
#'
#' @examples
#' set.seed(1)
#' env_cluster <- sample(1:10, size = 100, replace = TRUE)
#' folds <- make_cluster_folds(env_cluster = env_cluster, n_folds = 3, seed = 1)
#' check_cluster_folds(folds, env_cluster)
check_cluster_folds <- function(folds, env_cluster) {
  if (!is.list(folds) || length(folds) < 1) {
    stop("folds must be a non-empty list.")
  }
  if (length(env_cluster) < 1) {
    stop("env_cluster must have length >= 1.")
  }
  if (anyNA(env_cluster)) {
    stop("env_cluster must not contain missing values.")
  }
  
  n_row <- length(env_cluster)
  all_test_rows <- integer(0)
  all_test_clusters <- c()
  
  for (i in seq_along(folds)) {
    fold_i <- folds[[i]]
    
    needed_names <- c("fold_id", "test_clusters", "train_clusters", "test_rows", "train_rows")
    if (!all(needed_names %in% names(fold_i))) {
      stop("Fold ", i, " is missing required elements.")
    }
    
    # No overlap in train/test rows within fold
    if (length(intersect(fold_i$test_rows, fold_i$train_rows)) > 0) {
      stop("Fold ", i, " has overlapping test_rows and train_rows.")
    }
    
    # Rows should be in range
    if (length(fold_i$test_rows) > 0 &&
        (min(fold_i$test_rows) < 1 || max(fold_i$test_rows) > n_row)) {
      stop("Fold ", i, " has test_rows outside valid row range.")
    }
    if (length(fold_i$train_rows) > 0 &&
        (min(fold_i$train_rows) < 1 || max(fold_i$train_rows) > n_row)) {
      stop("Fold ", i, " has train_rows outside valid row range.")
    }
    
    if((length(fold_i$train_rows) + length(fold_i$test_rows)) != n_row) {
      stop('incorrect total rows in fold ', i)
    }
    
    # test_rows should match test_clusters
    expected_test_rows <- which(env_cluster %in% fold_i$test_clusters)
    if (!setequal(fold_i$test_rows, expected_test_rows)) {
      stop("Fold ", i, ": test_rows do not match test_clusters.")
    }
    
    # train_rows should match train_clusters
    expected_train_rows <- which(env_cluster %in% fold_i$train_clusters)
    if (!setequal(fold_i$train_rows, expected_train_rows)) {
      stop("Fold ", i, ": train_rows do not match train_clusters.")
    }
    
    # Train/test clusters should not overlap
    if (length(intersect(fold_i$test_clusters, fold_i$train_clusters)) > 0) {
      stop("Fold ", i, " has overlapping test_clusters and train_clusters.")
    }
    
    all_test_rows <- c(all_test_rows, fold_i$test_rows)
    all_test_clusters <- c(all_test_clusters, fold_i$test_clusters)
  }
  
  # Every row should appear in a test set exactly once
  if (!setequal(all_test_rows, seq_len(n_row))) {
    stop("Test rows across folds do not cover all rows exactly once.")
  }
  if (length(all_test_rows) != n_row) {
    stop("Some rows appear in test sets more than once.")
  }
  
  # Every cluster should appear in a test set exactly once
  unique_clusters <- unique(env_cluster)
  if (!setequal(all_test_clusters, unique_clusters)) {
    stop("Test clusters across folds do not cover all clusters exactly once.")
  }
  if (length(all_test_clusters) != length(unique_clusters)) {
    stop("Some clusters appear in test sets more than once.")
  }
  
  invisible(TRUE)
}



