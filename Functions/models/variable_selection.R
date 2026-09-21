# Functions for choosing a small set of predictors: PCA summaries,
# correlations, and subset quality from the subselect package (coverage = RM^2,
# the share of total variance explained by a set; comparable to the cumulative
# variance of the same number of PCs).


# PCA summaries -----------------------------------------------------------

#' Variance explained by each principal component
#'
#' @param pca A `prcomp` object.
#' @return Tibble with `pc_num`, `PC`, `variance`, `prop_var`, `cum_var`.
pca_variance <- function(pca) {
  v <- pca$sdev^2
  tibble::tibble(pc_num = seq_along(v),
                 PC = factor(paste0("PC", pc_num), levels = paste0("PC", pc_num)),
                 variance = v,
                 prop_var = v / sum(v),
                 cum_var = cumsum(prop_var))
}

#' PCA loadings in long format, variables ordered by clustering
#'
#' Variables are ordered by hierarchical clustering on their loadings, so
#' variables with similar loading patterns sit next to each other.
#'
#' @param pca A `prcomp` object.
#' @param n_pc Number of PCs to keep (default all).
#' @return Tibble with `variable` (factor), `PC` (factor), `loading`.
pca_loadings_long <- function(pca, n_pc = ncol(pca$rotation)) {
  L <- pca$rotation[, seq_len(n_pc), drop = FALSE]
  var_order <- rownames(L)[stats::hclust(stats::dist(L))$order]
  
  tibble::as_tibble(L, rownames = "variable") |>
    tidyr::pivot_longer(-variable, names_to = "PC", values_to = "loading") |>
    dplyr::mutate(PC = factor(PC, levels = colnames(L)),
                  variable = factor(variable, levels = var_order))
}

# correlation -------------------------------------------------------------

#' Correlation matrix in long format, variables ordered by clustering
#'
#' @param cor_mat Correlation matrix.
#' @return Tibble with `var1`, `var2` (factors ordered by hierarchical
#'   clustering on 1 - |r|) and `r`.
cor_long <- function(cor_mat) {
  ord <- rownames(cor_mat)[stats::hclust(stats::as.dist(1 - abs(cor_mat)))$order]
  
  tibble::as_tibble(cor_mat, rownames = "var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r") |>
    dplyr::mutate(var1 = factor(var1, levels = ord),
                  var2 = factor(var2, levels = ord))
}

# subsets (subselect package) --------------------------------------------
# Coverage of a set = RM^2 from subselect::rm.coef(): the share of the total
# variance of all variables explained by regressing them on the set (McCabe's
# second criterion for "principal variables"). See ?subselect::rm.coef.

#' Best k-variable subsets by RM, as a table
#'
#' Wraps `subselect::improve()` (criterion "RM"), which searches for the
#' k-variable subset with the highest RM for each k. It's a local search, so
#' not guaranteed optimal. See ?subselect::improve.
#'
#' `improve()` needs a non-singular correlation matrix, so variables that are
#' (near-)exact linear combinations of others must be passed in `exclude`
#' (see `subselect::trim.matrix()`).
#'
#' @param cor_mat Correlation matrix of all variables.
#' @param kmax Largest subset size.
#' @param include Variables forced into every subset (optional).
#' @param exclude Variables removed before the search (optional).
#' @return Tibble with `k`, `coverage` (RM^2) and `variables`.
best_subsets <- function(cor_mat, kmax, include = NULL, exclude = NULL) {
  stopifnot(all(exclude %in% colnames(cor_mat)))
  keep <- setdiff(colnames(cor_mat), exclude)
  cor_mat <- cor_mat[keep, keep]
  
  include_idx <- if (is.null(include)) NULL else match(include, colnames(cor_mat))
  stopifnot(!anyNA(include_idx))
  
  res <- subselect::improve(cor_mat, kmin = length(include) + 1,
                            kmax = kmax, include = include_idx,
                            criterion = "RM", setseed = TRUE)
  sets <- res$bestsets
  tibble::tibble(
    k = unname(rowSums(sets > 0)),
    coverage = unname(res$bestvalues)^2,
    variables = apply(sets, 1, \(i) paste(colnames(cor_mat)[i[i > 0]],
                                          collapse = ", "))
  )
}

#' Coverage, largest |r| within the set, and VIF for a set of variables
#'
#' @param dat Data frame of numeric variables (coverage uses all columns).
#' @param vars Variables in the set, in the order to add them.
#' @return Tibble with `variable`, `cum_coverage` (RM^2 of the set up to and
#'   including this variable), `max_abs_r` (with the rest of the set) and
#'   `vif` (1 / (1 - R^2) from `lm()` of the variable on the rest of the set).
summarise_set <- function(dat, vars) {
  stopifnot(all(vars %in% names(dat)))
  cor_mat <- stats::cor(dat)
  idx <- match(vars, colnames(cor_mat))
  
  r <- abs(cor_mat[vars, vars, drop = FALSE])
  diag(r) <- NA
  
  vif <- purrr::map_dbl(vars, \(v) {
    if (length(vars) == 1) return(1)
    fit <- stats::lm(stats::reformulate(setdiff(vars, v), v), data = dat)
    1 / (1 - summary(fit)$r.squared)
  })
  
  tibble::tibble(
    variable = vars,
    cum_coverage = purrr::map_dbl(seq_along(idx), \(i) {
      subselect::rm.coef(cor_mat, idx[1:i])^2
    }),
                 max_abs_r = apply(r, 1, max, na.rm = TRUE),
    vif = vif
  )
}
