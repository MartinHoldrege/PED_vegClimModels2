# Functions for choosing a small set of predictors.
#
# PCA summaries, correlations, and how well a subset of variables stands in for
# all of them. "Coverage" of a subset = RM^2 from the subselect package: the
# share of the total variance of all variables that is explained by regressing
# them on the subset (McCabe's criterion for "principal variables"). It's on
# the same scale as the cumulative variance explained by the same number of
# PCs, which is the best any k-dimensional summary can do.
# See ?subselect::rm.coef and ?subselect::improve.
#
# Each @examples block assigns every argument, so you can run it and then step
# through the function body line by line.


# PCA summaries -----------------------------------------------------------

#' Variance explained by each principal component
#'
#' @param pca A `prcomp` object.
#' @return Tibble with one row per PC: `pc_num`, `PC` (factor), `variance`
#'   (eigenvalue), `prop_var` (share of total variance) and `cum_var`
#'   (cumulative share).
#' @examples
#' pca <- prcomp(mtcars, scale. = TRUE)
#' pca_variance(pca = pca)
pca_variance <- function(pca) {
  variance <- pca$sdev^2
  pc_names <- colnames(pca$rotation)  # "PC1", "PC2", ...
  
  tibble(
    pc_num = seq_along(variance),
    PC = factor(pc_names, levels = pc_names),
    variance = variance,
    prop_var = variance / sum(variance),
    cum_var = cumsum(variance) / sum(variance)
  )
}

#' Number of PCs needed to reach a share of the total variance
#'
#' @param pca A `prcomp` object.
#' @param prop Share of total variance (0-1).
#' @return Integer.
#' @examples
#' pca <- prcomp(mtcars, scale. = TRUE)
#' prop <- 0.95
#' n_pcs_needed(pca = pca, prop = prop)
n_pcs_needed <- function(pca, prop) {
  cum_var <- pca_variance(pca)$cum_var
  min(which(cum_var >= prop))
}

#' PCA loadings in long format
#'
#' Variables are ordered by hierarchical clustering on their loadings, so
#' variables with similar loadings sit next to each other in plots.
#'
#' @param pca A `prcomp` object.
#' @param n_pc Number of PCs to keep.
#' @return Tibble with `variable` (factor), `PC` (factor) and `loading`.
#' @examples
#' pca <- prcomp(mtcars, scale. = TRUE)
#' n_pc <- 3
#' pca_loadings_long(pca = pca, n_pc = n_pc)
pca_loadings_long <- function(pca, n_pc) {
  loadings <- pca$rotation[, 1:n_pc, drop = FALSE]  # variables x PCs
  
  clusters <- hclust(dist(loadings))
  var_order <- rownames(loadings)[clusters$order]
  
  as_tibble(loadings, rownames = "variable") |>
    pivot_longer(-variable, names_to = "PC", values_to = "loading") |>
    mutate(PC = factor(PC, levels = colnames(loadings)),
           variable = factor(variable, levels = var_order))
}

#' Variables with the largest absolute loadings on each PC
#'
#' @param pca A `prcomp` object.
#' @param n_pc Number of PCs to show.
#' @param n_top Number of variables to show per PC.
#' @return Tibble with a `rank` column and one column per PC; entries are
#'   "variable (loading)".
#' @examples
#' pca <- prcomp(mtcars, scale. = TRUE)
#' n_pc <- 3
#' n_top <- 4
#' top_loadings(pca = pca, n_pc = n_pc, n_top = n_top)
top_loadings <- function(pca, n_pc, n_top) {
  pca_loadings_long(pca, n_pc) |>
    group_by(PC) |>
    slice_max(abs(loading), n = n_top, with_ties = FALSE) |>
    mutate(rank = row_number(),
           entry = paste0(variable, " (", sprintf("%+.2f", loading), ")")) |>
    ungroup() |>
    select(rank, PC, entry) |>
    pivot_wider(names_from = PC, values_from = entry)
}


# correlation -------------------------------------------------------------

#' Correlation matrix in long format, for a heatmap
#'
#' Variables are ordered by hierarchical clustering on 1 - |r|, so strongly
#' correlated variables (either sign) sit next to each other.
#'
#' @param cor_mat Correlation matrix.
#' @return Tibble with `var1`, `var2` (factors) and `r`.
#' @examples
#' cor_mat <- cor(mtcars)
#' cor_long(cor_mat = cor_mat)
cor_long <- function(cor_mat) {
  clusters <- hclust(as.dist(1 - abs(cor_mat)))
  var_order <- rownames(cor_mat)[clusters$order]
  
  as_tibble(cor_mat, rownames = "var1") |>
    pivot_longer(-var1, names_to = "var2", values_to = "r") |>
    mutate(var1 = factor(var1, levels = var_order),
           var2 = factor(var2, levels = var_order))
}

#' Variables that are (near-)exact linear combinations of others
#'
#' `subselect::improve()` stops if any variable is a (near-)exact linear
#' combination of others (e.g. MAT = mean of T_min and T_max). This names one
#' variable per such dependency, using `subselect::trim.matrix()` with the
#' tolerance `improve()` checks against. Empty if there are none.
#'
#' @param cor_mat Correlation matrix.
#' @return Character vector of variable names.
#' @examples
#' dat <- mtcars
#' dat$wt_hp <- dat$wt + dat$hp  # exact linear combination
#' cor_mat <- cor(dat)
#' linear_dependencies(cor_mat = cor_mat)
linear_dependencies <- function(cor_mat) {
  trimmed <- subselect::trim.matrix(cor_mat, tolval = 1000 * .Machine$double.eps)
  trimmed$names.discarded
}


# subsets -----------------------------------------------------------------

#' Best subsets of given sizes, by coverage
#'
#' For each subset size k, searches for the k variables with the highest
#' coverage (RM^2), using `subselect::improve()` with criterion "RM". It's a
#' local search from random starting subsets, so not guaranteed to find the
#' best subset; more starts (`n_starts`, improve()'s `nsol`) make a miss less
#' likely. The best subset of one size doesn't have to contain the best subset
#' of a smaller size.
#'
#' There are two ways to keep a variable out of the subsets. `exclude` removes
#' it from the correlation matrix, so it is neither selectable nor counted in
#' coverage (use this for linear dependencies). `exclude_improve` passes it to
#' `improve()`'s own `exclude`, so it can't be selected but coverage still
#' measures how well the subset explains it.
#'
#' @param cor_mat Correlation matrix of all variables.
#' @param sizes Subset sizes (k) to return, e.g. `4:8`. Each must be larger
#'   than the number of `include` variables.
#' @param include Variables that must be in every subset, or NULL.
#' @param exclude Variables removed from the correlation matrix, or NULL. Not
#'   selectable and not counted in coverage. Must cover any linear
#'   dependencies (see `linear_dependencies()`).
#' @param exclude_improve Variables not selectable but still counted in
#'   coverage, or NULL.
#' @param n_starts Number of random starting subsets per size.
#' @return Tibble with `k`, `coverage` and `variables` (comma-separated).
#' @examples
#' cor_mat <- cor(mtcars)
#' sizes <- 2:4
#' include <- "wt"
#' exclude <- NULL
#' exclude_improve <- "qsec"
#' n_starts <- 10
#' best_subsets(cor_mat = cor_mat, sizes = sizes, include = include,
#'              exclude = exclude, exclude_improve = exclude_improve,
#'              n_starts = n_starts)
best_subsets <- function(cor_mat, sizes, include = NULL, exclude = NULL,
                         exclude_improve = NULL,
                         n_starts = 10) {
  stopifnot(all(c(include, exclude, exclude_improve) %in% colnames(cor_mat)),
            !any(include %in% c(exclude, exclude_improve)),
            !any(exclude_improve %in% exclude))
  
  # remove excluded variables: neither selectable nor counted in coverage
  keep <- setdiff(colnames(cor_mat), exclude)
  cor_mat <- cor_mat[keep, keep]
  
  # each size must be bigger than the required set, and small enough that
  # enough selectable variables are left
  if (any(sizes <= length(include))) {
    stop("every size must be larger than the number of `include` variables (",
         length(include), ")")
  }
  if (any(sizes > ncol(cor_mat) - length(exclude_improve))) {
    stop("sizes exceed the number of selectable variables (",
         ncol(cor_mat) - length(exclude_improve), ")")
  }
  
  # improve() would stop on these; say which ones to exclude instead
  dependent <- linear_dependencies(cor_mat)
  if (length(dependent) > 0) {
    stop("Linearly dependent variables; add to `exclude`: ",
         paste(dependent, collapse = ", "))
  }
  
  # improve() takes column positions, not names
  include_idx <- if (is.null(include)) NULL else match(include, colnames(cor_mat))
  exclude_improve_idx <- if (is.null(exclude_improve)) NULL else {
    match(exclude_improve, colnames(cor_mat))
  }
  
  # searches every size from min(sizes) to max(sizes)
  result <- subselect::improve(cor_mat,
                               kmin = min(sizes), kmax = max(sizes),
                               include = include_idx, nsol = n_starts,
                               exclude = exclude_improve_idx,
                               criterion = "RM", setseed = TRUE)
  
  # result$bestsets: best subset of each size (one row per size), as column
  #   positions padded with 0
  # result$bestvalues: RM of each best subset (coverage is RM^2)
  all_sizes <- min(sizes):max(sizes)
  variables <- map_chr(seq_along(all_sizes), \(i) {
    idx <- result$bestsets[i, ]
    idx <- idx[idx > 0]
    paste(colnames(cor_mat)[idx], collapse = ", ")
  })
  
  out <- tibble(k = all_sizes,
                coverage = unname(result$bestvalues)^2,
                variables = variables)
  
  # keep only the requested sizes
  out[out$k %in% sizes, ]
}

#' Coverage, largest |r| and VIF for a chosen set of variables
#'
#' @param dat Data frame of numeric variables. Coverage is measured against
#'   all of its columns, so pass every candidate variable, not just the set.
#' @param vars Variables in the set, in the order to add them (at least 2).
#' @return Tibble with one row per variable in `vars`:
#'   `cum_coverage` (coverage of the set up to and including this variable),
#'   `max_abs_r` (largest |r| with another variable in the set) and
#'   `vif` (1 / (1 - R^2) from regressing this variable on the rest of the
#'   set).
#' @examples
#' dat <- mtcars
#' vars <- c("wt", "hp", "qsec")
#' summarise_set(dat = dat, vars = vars)
summarise_set <- function(dat, vars) {
  stopifnot(all(vars %in% names(dat)), length(vars) >= 2)
  cor_mat <- cor(dat)  # correlations among all variables in dat
  
  # coverage of the first i variables in the set
  # (rm.coef takes column positions, not names; it returns RM, so square it)
  cum_coverage <- map_dbl(seq_along(vars), \(i) {
    first_i <- match(vars[1:i], colnames(cor_mat))
    subselect::rm.coef(cor_mat, first_i)^2
  })
  
  # largest |r| between each variable and the rest of the set
  max_abs_r <- map_dbl(vars, \(v) {
    others <- setdiff(vars, v)
    max(abs(cor_mat[v, others]))
  })
  
  # VIF: 1 / (1 - R^2) from regressing each variable on the rest of the set
  vif <- map_dbl(vars, \(v) {
    others <- setdiff(vars, v)
    fit <- lm(reformulate(others, response = v), data = dat)
    1 / (1 - summary(fit)$r.squared)
  })
  
  tibble(variable = vars, cum_coverage, max_abs_r, vif)
}

#' Best subsets for several combinations of sizes and required variables
#'
#' Runs `best_subsets()` once per spec and stacks the results.
#'
#' @param cor_mat Correlation matrix of all variables.
#' @param specs List of specs, each a list with `sizes` (subset sizes),
#'   `include` (NULL or variables required in every subset) and, optionally,
#'   `exclude_improve` (variables that spec can't select but that coverage
#'   still measures).
#' @param exclude Variables removed from the correlation matrix for every
#'   search, or NULL.
#' @return Tibble with `required`, `excluded` (the spec's `exclude_improve`),
#'   `k`, `coverage` and `variables`.
#' @examples
#' cor_mat <- cor(mtcars)
#' specs <- list(
#'   list(sizes = 2:3, include = NULL),
#'   list(sizes = 3, include = c("wt", "hp"), exclude_improve = "qsec")
#' )
#' exclude <- NULL
#' best_subsets_by_spec(cor_mat = cor_mat, specs = specs, exclude = exclude)
best_subsets_by_spec <- function(cor_mat, specs, exclude = NULL) {
  map(specs, \(spec) {
    stopifnot(!is.null(spec$sizes))  # catches a misnamed element
    
    subsets <- best_subsets(cor_mat, sizes = spec$sizes,
                            include = spec$include, exclude = exclude,
                            exclude_improve = spec$exclude_improve)
    
    required <- if (is.null(spec$include)) "none" else paste(spec$include, collapse = ", ")
    excluded <- if (is.null(spec$exclude_improve)) "none" else {
      paste(spec$exclude_improve, collapse = ", ")
    }
    mutate(subsets, required = required, excluded = excluded, .before = 1)
  }) |>
    bind_rows()
}