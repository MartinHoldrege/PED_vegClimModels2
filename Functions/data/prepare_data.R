# various functions for preparing data for modeling and other tasks


# prep dataframe ----------------------------------------------------------

# for preparing version 1 of the biomass dataset, this 
# isn't the final version of the dataset to use
# but more for development/proof of concept
# trim_tree_cov--provide a proportional cover of trees below it's considered
# 0 (i.e. none zero tree predictions are leaking into areas
# where there probably aren't trees)
prepare_d01 <- function(data, cover_suffix, pfts, trim_tree_cov = NULL,
                        trim_shrub_cov = NULL) {
  dat2 <- data %>% 
    # for consistent downstream naming
    rename_with(.fn = \(x) str_replace(x, cover_suffix, 'Cov'),
                .cols = ends_with(cover_suffix)) %>% 
    filter(!is.na(biomass_MgPerHect)) %>% 
    rename(totalBio = biomass_MgPerHect)
  

  dat3 <- dat2
  # make cover sum to 1 -----------------------------------------------------
  # fix shouldn't be needed on 'final' data
  cov_cols <- c(paste0(pfts, "Cov"), "bareGroundCov")
  

  
  dat4 <- dat3 %>%
    mutate(.cov_total = rowSums(across(all_of(cov_cols)), na.rm = TRUE)) %>%
    mutate(across(all_of(cov_cols), ~ .x / .cov_total)) %>%
    select(-.cov_total)
  
  # correcting the totals
  if(!'totalTreeCov' %in% cov_cols) {
    dat4$totalTreeCov <- dat4$needleLeavedTreeCov + dat4$broadLeavedTreeCov
  }
  
  if(!'totalHerbaceousCov' %in% cov_cols) {
    dat4$totalHerbaceousCov <- with(dat4, C3GramCov + C4GramCov + ForbCov)
  }
  
  
  # set.seed(12)
  # dat2 <- sample_n(dat2, size = n_sample)
  if(!is.null(trim_tree_cov)) {
    
    tree_cols <- c("totalTreeCov", "needleLeavedTreeCov", "broadLeavedTreeCov")
    stopifnot(trim_tree_cov >=0,
              trim_tree_cov <= 1,
              tree_cols %in% names(dat3))
    
    # make small "noise' values 0
    dat4[dat4[["totalTreeCov"]] < trim_tree_cov, tree_cols] <- 0
  }
  
  if(!is.null(trim_shrub_cov)) {
    # make small "noise' values 0
    dat4$shrubCov[dat4$shrubCov < trim_shrub_cov] <- 0
  }
  
  # normalize potential predictor variables ---------------------------------
  
  # normalize predictors (all potential predictors, not just those in the model)
  # doing this before sim_bio so all data are on standardized scale
  
  vars_to_scale <- names(dat4)[!names(dat4) %in% c(
    # skip non-numeric / identifier columns
    "biomassSource", "biomassType", "year", 
    "burnedMoreThan20YearsAgo", "region",
    # skip cover columns (these are proportions, not climate)
    paste0(pfts, "Cov"), "bareGroundCov",
    "totalTreeCov", "totalHerbaceousCov",
    # skip response
    "totalBio",
    # skip coordinates (keep raw for spatial analyses)
    "x", "y"
  )]
  
  scale_df <- tibble::tibble(
    variable = vars_to_scale,
    mean = sapply(dat4[vars_to_scale], mean, na.rm = TRUE),
    sd = sapply(dat4[vars_to_scale], sd, na.rm = TRUE)
  )
  
  dat4 <- dat4 |>
    mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.x))))
  
  # add artificial regions --------------------------------------------------
  # based on climate clustering
  
  region_vars <- c("tmean_CLIM", 
                   "precip_CLIM", 
                   "x",
                   "y")
  
  reg <- make_region_kmeans(
    dat = dat4,
    vars = region_vars,
    nstart = 1,
    k = 20,
    seed = 42,
    max_iter = 10000
  )
  
  dat4$region <- reg$data
  
  # 'cleaning' 

  stopifnot(mean(dat4$totalBio == 0) < 0.01) # if there are many zeros that's suspect
  
  # making it a very small number so log is possible
  dat4$totalBio[dat4$totalBio == 0] <- min(dat4$totalBio[dat4$totalBio > 0])/2 
  
  list(data = dat4,
       scale = scale_df)
}


#' Convert standardized (z-score) variables back to original scale
#'
#' @param data Data frame containing standardized columns.
#' @param scale_df Data frame with columns `variable`, `mean`, `sd` — the
#'   standardization parameters used to create the z-scores.
#' @param vars Character vector of column names to back-convert. Default:
#'   all variables present in both `data` and `scale_df`.
#'
#' @return The input data frame with specified columns back on their
#'   original scale.
#' @examples
#' df <- data.frame(tmean = c(-1, 0, 1), precip = c(0.5, -0.5, 0))
#' scale_df <- data.frame(
#'   variable = c("tmean", "precip"),
#'   mean = c(10, 500),
#'   sd = c(5, 200)
#' )
#' unstandardize(df, scale_df)
#' @export
unstandardize <- function(data, scale_df, vars = NULL) {
  stopifnot(all(c("variable", "mean", "sd") %in% names(scale_df)))
  
  if (is.null(vars)) {
    vars <- intersect(names(data), scale_df$variable)
  }
  stopifnot(all(vars %in% names(data)),
            all(vars %in% scale_df$variable))
  
  for (v in vars) {
    row <- scale_df[scale_df$variable == v, ]
    data[[v]] <- data[[v]] * row$sd + row$mean
  }
  
  data
}


#' Back-transform z-scored values in long format
#'
#' Reverse z-score standardization for a long-format data frame where
#' variable names are stored in a column rather than as separate columns.
#' Rows with variable names not found in `scale_df` are dropped with a
#' warning.
#'
#' @param data Data frame with a column of variable names and a column of
#'   standardized values.
#' @param scale_df Data frame with columns `variable`, `mean`, `sd`.
#' @param name_col Character; column in `data` containing variable names.
#'   Default `"name"`.
#' @param value_col Character; column in `data` containing z-scored values.
#'   Default `"mean_value"`.
#'
#' @return The input data frame with `value_col` back on the original scale.
#'   Rows with variable names absent from `scale_df` are excluded.
#' @examples
#' df <- data.frame(
#'   name = c("tmean", "tmean", "precip", "precip"),
#'   mean_value = c(-1, 1, -0.5, 0.5)
#' )
#' scale_df <- data.frame(
#'   variable = c("tmean", "precip"),
#'   mean = c(10, 500),
#'   sd = c(5, 200)
#' )
#' unstandardize_long(df, scale_df)
#' @export
unstandardize_long <- function(data, scale_df,
                               name_col = "name",
                               value_col = "mean_value") {
  stopifnot(
    all(c("variable", "mean", "sd") %in% names(scale_df)),
    name_col %in% names(data),
    value_col %in% names(data)
  )
  
  unmatched <- setdiff(unique(data[[name_col]]), scale_df$variable)
  if (length(unmatched) > 0) {
    warning("Dropping rows with variables not in scale_df: ",
            paste(unmatched, collapse = ", "))
    data <- data[!data[[name_col]] %in% unmatched, ]
  }
  
  data <- dplyr::left_join(data, scale_df, by = setNames("variable", name_col))
  
  data[[value_col]] <- data[[value_col]] * data[["sd"]] + data[["mean"]]
  
  data[["mean"]] <- NULL
  data[["sd"]] <- NULL
  
  data
}
