# selected 'purer' pixels for model fitting




# Select 'purer' pixels for model fitting


#' Select near-pure pixels by region for model fitting
#'
#' For each region and PFT, identifies pixels where that PFT dominates
#' (top quantile of scaled cover). Returns a de-duplicated set of pixels
#' selected as "pure" for at least one PFT. PFTs that never exceed
#' `min_raw_cover` in a region are skipped for that region.
#'
#' @param dat A data.frame/tibble.
#' @param cover_vars Character vector of cover column names (must contain "Cov").
#' @param region_col Character; column name defining regions.
#' @param q Numeric in (0, 1); quantile cutoff for scaled cover (default 0.9).
#' @param min_raw_cover Numeric; minimum raw cover to be eligible (default 0.05).
#' @param max_n_per_region Integer or NULL; if set, randomly downsample regions
#'   with more than this many selected pixels. Helps prevent data-rich regions
#'   from dominating the training set.
#' @param id_col Optional character; unique row ID column. If NULL, row numbers
#'   are used.
#' @param seed Integer or NULL; seed for reproducible downsampling.
#'
#' @return A list with:
#' \describe{
#'   \item{data}{De-duplicated selected rows with all original columns plus
#'     a `selected_for` column (comma-separated PFT names the pixel qualified
#'     for).}
#'   \item{thresholds}{Data frame of scaled-cover thresholds by region and PFT.}
#'   \item{selection_summary}{Data frame summarizing n selected per region and
#'     PFT, before and after any downsampling.}
#' }
#' @examples
#' set.seed(1)
#' n <- 200
#' dat <- dplyr::tibble(
#'   region = sample(1:3, n, replace = TRUE),
#'   treeCov = runif(n, 0, 0.6),
#'   shrubCov = runif(n, 0, 0.4),
#'   grassCov = runif(n, 0, 0.3),
#'   totalBio = rlnorm(n, 3, 0.5)
#' )
#'
#' res <- select_purer_by_region(
#'   dat = dat,
#'   cover_vars = c("treeCov", "shrubCov", "grassCov"),
#'   region_col = "region",
#'   q = 0.9,
#'   min_raw_cover = 0.05
#' )
#'
#' # selected pixels (de-duplicated)
#' nrow(res$data)
#'
#' # which PFTs each pixel qualified for
#' table(res$data$selected_for)
#'
#' # thresholds used per region x PFT
#' res$thresholds
#'
#' # with a per-region cap
#' res2 <- select_purer_by_region(
#'   dat = dat,
#'   cover_vars = c("treeCov", "shrubCov", "grassCov"),
#'   region_col = "region",
#'   max_n_per_region = 10,
#'   seed = 42
#' )
#' res2$data |> dplyr::count(region)
#' @export
select_purer_by_region <- function(dat,
                                   cover_vars,
                                   region_col = "region",
                                   q = 0.9,
                                   min_raw_cover = 0.05,
                                   max_n_per_region = NULL,
                                   id_col = NULL,
                                   seed = NULL) {
  
  # --- input checks --------------------------------------------------------
  stopifnot(
    is.data.frame(dat),
    length(cover_vars) >= 1,
    all(cover_vars %in% names(dat)),
    all(grepl("Cov", cover_vars)),
    region_col %in% names(dat),
    is.numeric(q), length(q) == 1, q >= 0, q < 1,
    is.numeric(min_raw_cover), length(min_raw_cover) == 1,
    min_raw_cover >= 0, min_raw_cover <= 1
  )
  
  if (!is.null(max_n_per_region)) {
    stopifnot(is.numeric(max_n_per_region), length(max_n_per_region) == 1,
              max_n_per_region > 0)
    max_n_per_region <- as.integer(max_n_per_region)
  }
  
  dat <- dplyr::as_tibble(dat)
  
  # i.e. no 'filtering'
  if(q == 0) {
    return(dat)
  }
  
  # --- row IDs -------------------------------------------------------------
  if (is.null(id_col)) {
    id_col <- ".row_id"
    dat <- dplyr::mutate(dat, .row_id = dplyr::row_number())
  } else {
    stopifnot(id_col %in% names(dat))
  }
  
  # --- scale covers to sum to 1 (across PFTs only, ignoring bare ground) ---
  cov_mat <- as.matrix(dat[, cover_vars])
  cov_sum <- rowSums(cov_mat, na.rm = TRUE)
  
  # avoid division by zero
  valid <- is.finite(cov_sum) & cov_sum > 0
  scaled_mat <- matrix(NA_real_, nrow = nrow(cov_mat), ncol = ncol(cov_mat))
  scaled_mat[valid, ] <- cov_mat[valid, , drop = FALSE] / cov_sum[valid]
  colnames(scaled_mat) <- paste0(cover_vars, "_scaled")
  
  # --- find thresholds and select pixels -----------------------------------
  # work in long form: one row per pixel x PFT
  regions <- dat[[region_col]]
  row_ids <- dat[[id_col]]
  
  # build long data frame efficiently
  n <- nrow(dat)
  n_pft <- length(cover_vars)
  
  long <- dplyr::tibble(
    .id = rep(row_ids, times = n_pft),
    region = rep(regions, times = n_pft),
    pft = rep(cover_vars, each = n),
    raw_cover = as.vector(cov_mat),
    scaled_cover = as.vector(scaled_mat)
  )
  names(long)[names(long) == "region"] <- region_col
  names(long)[names(long) == ".id"] <- id_col
  
  # filter to rows eligible for selection (meet min_raw_cover)
  eligible <- long |>
    dplyr::filter(!is.na(.data$scaled_cover),
                  .data$raw_cover >= min_raw_cover)
  
  # compute thresholds per region x PFT
  thresholds <- eligible |>
    dplyr::group_by(.data[[region_col]], .data$pft) |>
    dplyr::summarise(
      threshold_scaled = stats::quantile(.data$scaled_cover, probs = q,
                                         na.rm = TRUE),
      n_eligible = dplyr::n(),
      max_raw_cover = max(.data$raw_cover, na.rm = TRUE),
      .groups = "drop"
    )
  
  # select pixels above threshold
  selected_long <- eligible |>
    dplyr::inner_join(thresholds, by = c(region_col, "pft")) |>
    dplyr::filter(.data$scaled_cover >= .data$threshold_scaled)
  
  # --- selection summary (before downsampling) -----------------------------
  selection_summary <- selected_long |>
    dplyr::group_by(.data[[region_col]], .data$pft) |>
    dplyr::summarise(n_selected = dplyr::n(), .groups = "drop")
  
  # --- de-duplicate: one row per pixel, tag which PFTs it qualified for ----
  selected_tags <- selected_long |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      selected_for = paste(sort(unique(.data$pft)), collapse = ", "),
      .groups = "drop"
    )
  
  # --- optional: cap per region -------------------------------------------
  if (!is.null(max_n_per_region)) {
    # need region info for downsampling
    selected_tags <- selected_tags |>
      dplyr::left_join(
        dat |> dplyr::select(dplyr::all_of(c(id_col, region_col))),
        by = id_col
      )
    
    if (!is.null(seed)) set.seed(seed)
    
    selected_tags <- selected_tags |>
      dplyr::group_by(.data[[region_col]]) |>
      dplyr::slice_sample(n = max_n_per_region) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::all_of(region_col))
  }
  
  # --- join back to full data ----------------------------------------------
  selected_data <- selected_tags |>
    dplyr::left_join(dat, by = id_col)
  

  # clean up temporary row id if we created it
  if (id_col == ".row_id") {
    selected_data <- dplyr::select(selected_data, -dplyr::all_of(".row_id"))
  }
  
  list(
    data = selected_data,
    thresholds = thresholds,
    selection_summary = selection_summary
  )
}


#' Filter data for alpha pre-estimation
#'
#' Identifies pixels where specified cover variables are below a threshold
#' (e.g., tree-free pixels for estimating herbaceous alphas).
#'
#' @param data Data frame with cover columns.
#' @param exclude_cols Character vector of cover column names that should
#'   be near-zero in the output subset.
#' @param max_cover Numeric; either length 1 (same threshold for all
#'   exclude_cols) or same length as exclude_cols (per-column thresholds).
#'
#' @return A filtered data frame containing only rows where all excluded
#'   cover columns are at or below their threshold.
#' @examples
#' # pixels with essentially no trees
#' herb_pixels <- filter_for_alpha(
#'   data = dat_train,
#'   exclude_cols = c("needleLeavedTreeCov", "broadLeavedTreeCov"),
#'   max_cover = 0.01
#' )
#'
#' # different thresholds per column
#' herb_pixels <- filter_for_alpha(
#'   data = dat_train,
#'   exclude_cols = c("needleLeavedTreeCov", "broadLeavedTreeCov", "shrubCov"),
#'   max_cover = c(0.01, 0.01, 0.05)
#' )
#' @export
filter_for_alpha <- function(data, exclude_cols, max_cover,
                             verbose = FALSE) {
  stopifnot(is.data.frame(data),
            all(exclude_cols %in% names(data)))
  
  if (length(max_cover) == 1) {
    max_cover <- rep(max_cover, length(exclude_cols))
  }
  stopifnot(length(max_cover) == length(exclude_cols))
  
  mask <- rep(TRUE, nrow(data))
  for (j in seq_along(exclude_cols)) {
    mask <- mask & data[[exclude_cols[j]]] <= max_cover[j]
  }
  
  out <- data[mask, , drop = FALSE]
  
  if (nrow(out) == 0) {
    stop("No pixels remain after filtering. Consider relaxing max_cover.")
  }
  
  if(verbose) {
    cat("filter_for_alpha:\n")
    cat("  Filters:", paste(exclude_cols, "<=", max_cover, collapse = ", "), "\n")
    cat("  N input:", nrow(data), "\n")
    cat("  N output:", nrow(out), "\n")
  }

  
  out
}

#' Select training pixels — purer selection or random subsample (i.e. a wrapper)
#'
#' Dispatches to `select_purer_by_region()` when `q` is provided,
#' or performs a simple random subsample when `n_sample` is provided.
#'
#' @param dat Data frame of candidate pixels.
#' @param cover_vars Character vector of cover column names.
#' @param purer_spec List with either purer params (`q`, `min_raw_cover`)
#'   or subsample params (`n_sample`). If `q` is NULL, subsampling is used.
#' @param region_col Character; column name for regions (used by purer).
#' @param seed Integer; random seed.
#'
#' @return A list with `$data` (selected data frame) and metadata.
#' @export
select_training_pixels <- function(dat, cover_vars = NULL, purer_spec,
                                   region_col = "region",
                                   seed = 42) {
  
  q <- purer_spec$q
  
  if (is.null(q)) {
    # simple random subsample
    n_sample <- purer_spec$n_sample
    stopifnot(!is.null(n_sample), n_sample > 0)
    
    n <- min(n_sample, nrow(dat))
    set.seed(seed)
    idx <- sample.int(nrow(dat), n)
    
    list(
      data = dat[idx, , drop = FALSE],
      thresholds = NULL,
      selection_summary = NULL,
      method = "subsample",
      n_selected = n
    )
  } else {
    # purer selection
    purer <- select_purer_by_region(
      dat = dat,
      cover_vars = cover_vars,
      region_col = region_col,
      q = q,
      min_raw_cover = purer_spec$min_raw_cover,
      seed = seed
    )
    
    purer$method <- "purer"
    purer$n_selected <- nrow(purer$data)
    purer
  }
}
