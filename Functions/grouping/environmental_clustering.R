# functions for making clusters based on environmental variables



#' Make environmental clusters with k-means
#'
#' Creates row-level environmental cluster assignments using k-means on
#' standardized predictor variables.
#'
#' @param data Data frame containing environmental variables.
#' @param vars Character vector of column names used for clustering.
#' @param k Integer number of clusters.
#' @param nstart Integer passed to `stats::kmeans()`.
#' @param iter.max Integer passed to `stats::kmeans()`.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list with elements:
#' \describe{
#'   \item{env_cluster}{Integer vector of cluster assignments, one per row of `data`.}
#'   \item{vars}{Character vector of variables used for clustering.}
#'   \item{k}{Number of clusters.}
#'   \item{centers}{Cluster centers on the standardized scale.}
#'   \item{scale_center}{Named numeric vector of means used for standardization.}
#'   \item{scale_scale}{Named numeric vector of standard deviations used for standardization.}
#' }
#'
#' @examples
#' set.seed(1)
#' n = 100
#' dat <- data.frame(
#'   tmean = rnorm(n),
#'   ppt = rnorm(n, 100, 20),
#'   vpd = rnorm(n, 2, 0.5)
#' )
#'
#' clust <- make_env_clusters(
#'   data = dat,
#'   vars = c("tmean", "ppt", "vpd"),
#'   k = 5,
#'   seed = 1
#' )
#'
#' head(clust$env_cluster)
make_env_clusters <- function(data,
                              vars,
                              k,
                              nstart = 20,
                              iter.max = 100,
                              seed = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  if (!is.character(vars) || length(vars) < 1) {
    stop("vars must be a character vector of length >= 1.")
  }
  if (!all(vars %in% names(data))) {
    stop("all vars must be column names in data.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != as.integer(k)) {
    stop("k must be a positive integer.")
  }
  if (!is.numeric(nstart) || length(nstart) != 1 || nstart < 1) {
    stop("nstart must be a positive integer.")
  }
  if (!is.numeric(iter.max) || length(iter.max) != 1 || iter.max < 1) {
    stop("iter.max must be a positive integer.")
  }

  x <- data[, vars, drop = FALSE]

  # Ensure all clustering variables are numeric
  is_num <- vapply(x, is.numeric, logical(1))
  if (!all(is_num)) {
    stop("all clustering variables must be numeric.")
  }

  # Do not allow missing values for now
  if (anyNA(x)) {
    stop("clustering variables contain missing values.")
  }

  # Check for zero-variance columns
  sds <- vapply(x, stats::sd, numeric(1))
  if (any(sds == 0)) {
    bad_vars <- names(sds)[sds == 0]
    stop("clustering variables must have non-zero variance. Problem variables: ",
         paste(bad_vars, collapse = ", "))
  }

  # Standardize variables
  x_scaled <- scale(x)
  scale_center <- attr(x_scaled, "scaled:center")
  scale_scale <- attr(x_scaled, "scaled:scale")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  km <- stats::kmeans(
    x = x_scaled,
    centers = k,
    nstart = nstart,
    iter.max = iter.max
  )

  out <- list(
    env_cluster = km$cluster,
    vars = vars,
    k = k,
    centers = km$centers,
    scale_center = scale_center,
    scale_scale = scale_scale
  )

  class(out) <- "env_clusters"
  out
}

#' Assign new data to existing environmental clusters
#'
#' Uses the saved k-means centers and scaling parameters from
#' `make_env_clusters()` to assign each row of new data to its
#' nearest cluster.
#'
#' @param data Data frame containing the clustering variables.
#' @param clustering Output of `make_env_clusters()`, with elements
#'   `centers`, `scale_center`, `scale_scale`, `vars`.
#'
#' @return Integer vector of cluster assignments, one per row of `data`.
#' @export
assign_to_clusters <- function(data, clustering) {
  stopifnot(all(clustering$vars %in% names(data)))
  
  x <- scale(
    as.matrix(data[, clustering$vars, drop = FALSE]),
    center = clustering$scale_center,
    scale = clustering$scale_scale
  )
  
  centers <- clustering$centers
  
  # for each row, find the index of the nearest cluster center
  apply(x, 1, function(row) {
    which.min(rowSums(sweep(centers, 2, row)^2))
  })
}

#' Map cluster assignments to CV fold test sets
#'
#' Given cluster assignments and fold definitions, determines which
#' fold each observation belongs to as a test observation.
#'
#' @param clusters Integer vector of cluster assignments.
#' @param folds Fold object from `make_cluster_folds()`.
#'
#' @return Integer vector of fold assignments (which fold's test set
#'   each observation falls in). Length equals length of `clusters`.
#' @export
clusters_to_fold_id <- function(clusters, folds) {
  fold_id <- rep(NA_integer_, length(clusters))
  
  for (fold in folds) {
    fold_id[clusters %in% fold$test_clusters] <- fold$fold_id
  }
  
  if (any(is.na(fold_id))) {
    warning(sum(is.na(fold_id)), " observations not assigned to any fold ",
            "(cluster not in any test set).")
  }
  
  fold_id
}


#' Create k-means region labels from selected variables
#'
#'
#' for making artificial regions--for 
# sampling purer cover pixels will occur within regions
# so it's stratified
# Note--this is initially for making regions in simulated data,
# but could be useful elsewhere
# I'll be using the name 'region' for these, whether
# defined by clustering or by externally defined (e.g. EPA) 
#'
#'
#'
#'
#' Assigns cluster-based region labels using k-means on scaled
#' versions of specified variables. Rows with missing values in
#' clustering variables receive `NA` for the region label.
#'
#' @param dat A data.frame or tibble containing clustering variables.
#' @param vars Character vector of column names used for clustering.
#' @param k Integer; number of clusters.
#' @param seed Integer; random seed for k-means initialization.
#' @param nstart Integer; number of random starts for k-means.
#' @param max_iter Integer; maximum number of k-means iterations.
make_region_kmeans <- function(dat,
                               vars,
                               k = 10,
                               seed = 1,
                               nstart = 10,
                               max_iter = 1000) {
  stopifnot(all(vars %in% names(dat)))
  
  x <- dat %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.data.frame()
  
  # drop rows with NA in clustering vars (keep mapping back)
  ok <- stats::complete.cases(x)
  x_ok <- x[ok, , drop = FALSE]
  
  # z-score within this call (so it can later be done within CV training folds)
  x_scaled <- scale(x_ok)
  centers <- attr(x_scaled, "scaled:center")
  scales  <- attr(x_scaled, "scaled:scale")
  
  set.seed(seed)
  km <- stats::kmeans(x_scaled, centers = k, nstart = nstart, iter.max = max_iter)
  
  region <- rep(NA_integer_, nrow(dat))
  region[ok] <- km$cluster
  

  
  list(
    data = region,
    info = list(
      vars = vars,
      k = k,
      centers = centers,
      scales = scales,
      kmeans = km
    )
  )
}



