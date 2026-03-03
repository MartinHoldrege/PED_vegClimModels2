# selected 'purer' pixels for model fitting




#' Select "near-pure" pixels by region and cover deciles
#'
#' Rescales `cover_vars` within each row to sum to 1 (ignoring bare ground and
#' other cover types). Within each region and for each cover variable, selects
#' rows in the top `q` quantile of the *scaled* cover, additionally requiring the
#' *raw* cover to be at least `min_raw_cover`. Rows may be selected for multiple
#' variables (returned as duplicated rows with a `selected_var` tag).
#'
#' @param dat A data.frame/tibble.
#' @param cover_vars Character vector of cover column names to use (e.g., 6 PFT covers).
#' @param region_col Character; region column name.
#' @param q Numeric in (0,1); quantile cutoff (default 0.9 = top decile).
#' @param min_raw_cover Numeric; minimum raw cover required for eligibility
#' Idea here is that lots of noise if really low cover (even if pure)
#' @param id_col Optional character; unique row ID column. If `NULL`, a temporary
#'   `.row_id` is created.
#' @param keep_cols Optional character vector of columns to retain (defaults to all).
#' @param scaled_suffix Character; suffix for scaled cover columns.
#' @param min_sum Numeric; rows with sum(cover_vars) < min_sum get NA scaled covers
#'   and are not eligible for selection.
#' @param duplicates logical, whether the returned dataframe includes duplicates
#' (if yes, then if a pixel is 'purer' for multiple variables then
#' it will show up multiple times)
#'
#' @return A list with:
#' \describe{
#'   \item{data}{Selected rows, #'   columns `selected_var`, `selected_cover_scaled`, and scaled cover columns.}
#'   \item{thresholds}{Quantile thresholds by region and variable.}
#' }
#' @export
select_purer_by_region <- function(dat,
                                   cover_vars,
                                   region_col = "region",
                                   q = 0.9,
                                   min_raw_cover = 0.05,
                                   id_col = NULL,
                                   keep_cols = NULL,
                                   min_sum = 0,
                                   duplicates = FALSE
                                    ) {
  stopifnot(is.data.frame(dat),
            length(cover_vars) >= 1,
            all(cover_vars %in% names(dat)),
            region_col %in% names(dat),
            is.numeric(q), length(q) == 1, q > 0, q < 1,
            is.numeric(min_raw_cover), 
            length(min_raw_cover) == 1, 
            min_raw_cover >= 0, 
            min_raw_cover <= 1,
            all_of(str_detect(cover_vars, 'Cov')))
  scaled_suffix = "_scaled"
  dat2 <- dplyr::as_tibble(dat)
  
  # ensure we have an ID
  if (is.null(id_col)) {
    id_col <- ".row_id"
    dat2 <- dplyr::mutate(dat2, .row_id = dplyr::row_number())
  } else {
    stopifnot(id_col %in% names(dat2))
  }
  
  # choose columns to keep
  if (!is.null(keep_cols)) {
    keep_cols <- unique(c(keep_cols, id_col, region_col, cover_vars))
    dat2 <- dplyr::select(dat2, dplyr::all_of(keep_cols))
  }
  
  # raw cover matrix
  cov_mat <- dat2 %>%
    dplyr::select(dplyr::all_of(cover_vars)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.data.frame()
  
  cov_sum <- rowSums(cov_mat, na.rm = TRUE)
  
  # rescale cover_vars within row to sum to 1 (ignoring bare ground)
  scaled <- cov_mat
  eligible_sum <- is.finite(cov_sum) & cov_sum > min_sum
  scaled[eligible_sum, ] <- scaled[eligible_sum, , drop = FALSE] / cov_sum[eligible_sum]
  scaled[!eligible_sum, ] <- NA_real_
  
  scaled_names <- paste0(cover_vars, scaled_suffix)
  names(scaled) <- scaled_names
  
  dat2 <- dplyr::bind_cols(dat2, dplyr::as_tibble(scaled))
  
  # long form includes BOTH scaled and raw for filtering
  long <- dat2 %>% 
    dplyr::select(dplyr::all_of(c(id_col, region_col,
                                  cover_vars, scaled_names))) %>% 
    pivot_longer(
      cols = all_of(c(cover_vars, scaled_names)),
      names_to  = c("pft", ".value"),
      names_pattern = "^(.*)(Cov.*)$",
      values_drop_na = FALSE
    ) 
  # CONTINUE HERE
  long %>% 
    dplyr::filter(!is.na(Cov),
           !is.na(Cov_scaled)) %>% 
    dplyr::rename(selected_var = pft) %>% 
    dplyr::group_by(.data[[region_col]], .data$selected_var) %>% 
    dplyr::summarise(
      threshold = stats::quantile(.data$Cov_scaled, probs = q, 
                                  na.rm = TRUE, type = 7),
      n_eligible = dplyr::n(),
      .groups = "drop"
    )
  
  long <- dplyr::bind_cols(
    dat2 %>% dplyr::select(dplyr::all_of(c(id_col, region_col))),
    dplyr::as_tibble(cov_mat),
    dat2 %>% dplyr::select(dplyr::all_of(scaled_names))
  ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cover_vars),
      names_to = "selected_var",
      values_to = "selected_cover_raw"
    ) %>%
    dplyr::left_join(
      tidyr::pivot_longer(
        dat2 %>% dplyr::select(dplyr::all_of(c(id_col, region_col, scaled_names))),
        cols = dplyr::all_of(scaled_names),
        names_to = "selected_var_scaled",
        values_to = "selected_cover_scaled"
      ) %>%
        dplyr::mutate(selected_var = sub(paste0(scaled_suffix, "$"), "", .data$selected_var_scaled)) %>%
        dplyr::select(dplyr::all_of(c(id_col, region_col, "selected_var", "selected_cover_scaled"))),
      by = c(id_col, region_col, "selected_var")
    )
  
  # thresholds by region x variable using *scaled* cover, among eligible rows
  thresholds <- long %>%
    dplyr::filter(!is.na(.data$selected_cover_scaled)) %>%
    dplyr::filter(!is.na(.data$selected_cover_raw)) %>%
    dplyr::filter(.data$selected_cover_raw >= min_raw_cover) %>%
    dplyr::group_by(.data[[region_col]], .data$selected_var) %>%
    dplyr::summarise(
      threshold = stats::quantile(.data$selected_cover_scaled, probs = q, 
                                  na.rm = TRUE, type = 7),
      n_eligible = dplyr::n(),
      .groups = "drop"
    )
  
  # select top decile within region x variable (allow duplicates across variables)
  selected_keys <- long %>%
    dplyr::inner_join(thresholds, by = c(region_col, "selected_var")) %>%
    dplyr::filter(!is.na(.data$selected_cover_scaled)) %>%
    dplyr::filter(!is.na(.data$selected_cover_raw)) %>%
    dplyr::filter(.data$selected_cover_raw >= min_raw_cover) %>%
    dplyr::filter(.data$selected_cover_scaled >= .data$threshold) %>%
    dplyr::select(dplyr::all_of(c(id_col, "selected_var", "selected_cover_raw",
                                  "selected_cover_scaled")))
  
  selected_data <- selected_keys %>%
    dplyr::left_join(dat2, by = id_col)
  
  list(
    data = selected_data,
    thresholds = thresholds
  )
}


# testing the function:

# making a toy dataset
set.seed(1)
n <- 60
pfts <- c("treeCov", "shrubCov", "forbCov",
          "c3grassCov", "c4grassCov", "annualCov")
# Simulate arbitrary cover values (not summing to 1)
cover_mat <- matrix(runif(n * length(pfts), 0, 0.6),
                    nrow = n, ncol = length(pfts))
colnames(cover_mat) <- pfts
# Bare ground
bare <- runif(n, 0, 0.8)
# Fake regions
region <- sample(1:3, size = n, replace = TRUE)
# Some arbitrary climate predictors
dat_test <- tibble::tibble(
  region = region,
  tmean = rnorm(n, 10, 5),
  precip = rnorm(n, 500, 100),
  VPD = rnorm(n, 1, 0.3),
  sand = runif(n, 0, 100)
) %>%
  dplyr::bind_cols(as.data.frame(cover_mat)) %>%
  dplyr::mutate(bareGroundCov = bare)

res <- select_purer_by_region(
  dat = dat_test,
  cover_vars = pfts,
  region_col = "region",
  q = 0.9,
  min_raw_cover = 0.05
)
