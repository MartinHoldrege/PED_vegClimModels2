# Functions for making artificial regions--for 
# sampling purer cover pixels will occur within regions
# so it's stratified
# Note--this is initially for making regions in simulated data,
# but could be useful elsewhere
# I'll be using the name 'region' for these, whether
# defined by clustering or by externally defined (e.g. EPA) 



#' Create k-means region labels from selected variables
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



# Optional quick sanity check
# table(dat4$region, useNA = "ifany")