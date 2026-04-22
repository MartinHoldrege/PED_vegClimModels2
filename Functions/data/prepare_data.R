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


# standardization ---------------------------------------------------------


#' Standardize (z-score) variables in a data frame
#'
#' Subtracts the mean and divides by the standard deviation for each
#' specified variable. Optionally computes and returns the scaling
#' parameters.
#'
#' @param data Data frame containing columns to standardize.
#' @param vars Character vector of column names to standardize. Default:
#'   all numeric columns.
#' @param scale_df Optional data frame with columns `variable`, `mean`, `sd`.
#'   If NULL, scaling parameters are computed from `data`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{data}{The input data frame with specified columns standardized.}
#'     \item{scale_df}{Data frame of scaling parameters used.}
#'   }
#' @examples
#' df <- data.frame(tmean = c(5, 10, 15), precip = c(200, 500, 800))
#' result <- standardize(df, vars = c("tmean", "precip"))
#' result$data
#' result$scale_df
#' unstandardize(result$data, result$scale_df)  # roundtrip
#' @export
standardize <- function(data, vars = NULL, scale_df = NULL) {
  if (is.null(vars)) {
    vars <- names(data)[sapply(data, is.numeric)]
  }
  stopifnot(all(vars %in% names(data)))
  
  if (is.null(scale_df)) {
    scale_df <- purrr::map_dfr(vars, function(v) {
      tibble::tibble(
        variable = v,
        mean = mean(data[[v]], na.rm = TRUE),
        sd = sd(data[[v]], na.rm = TRUE)
      )
    })
  }
  
  stopifnot(all(c("variable", "mean", "sd") %in% names(scale_df)),
            all(vars %in% scale_df$variable))
  
  for (v in vars) {
    row <- scale_df[scale_df$variable == v, ]
    data[[v]] <- (data[[v]] - row$mean) / row$sd
  }
  
  list(data = data, scale_df = scale_df)
}

#' Standardize (z-score) layers of a SpatRaster
#'
#' Applies z-score standardization to specified layers of a SpatRaster
#' using pre-computed scaling parameters. Does not load the full raster
#' into memory — uses lazy raster math.
#'
#' @param rast A `SpatRaster`.
#' @param scale_df Data frame with columns `variable`, `mean`, `sd`.
#'   Variable names must match layer names in `rast`.
#' @param vars Character vector of layer names to standardize. Default:
#'   all layers whose names appear in `scale_df`.
#'
#' @return A `SpatRaster` with specified layers standardized.
#' @examples
#' # r <- read_climate_raster(path)
#' # r_std <- standardize_raster(r, scale_df, vars = c("MAT", "MAP"))
#' @export
standardize_raster <- function(rast, scale_df, vars = NULL) {
  stopifnot(inherits(rast, "SpatRaster"),
            all(c("variable", "mean", "sd") %in% names(scale_df)))
  
  if (is.null(vars)) {
    vars <- intersect(names(rast), scale_df$variable)
  }
  
  missing_in_rast <- setdiff(vars, names(rast))
  missing_in_scale <- setdiff(vars, scale_df$variable)
  if (length(missing_in_rast) > 0) {
    stop("Layers not in raster: ", paste(missing_in_rast, collapse = ", "))
  }
  if (length(missing_in_scale) > 0) {
    stop("Variables not in scale_df: ", paste(missing_in_scale, collapse = ", "))
  }
  
  for (v in vars) {
    row <- scale_df[scale_df$variable == v, ]
    rast[[v]] <- (rast[[v]] - row$mean) / row$sd
  }
  
  rast
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

#' Compute or load cached scale_df for a climate raster
#'
#' Computes mean and SD for specified layers using `terra::global`. Caches
#' the result alongside a fingerprint of the source file. On subsequent
#' calls, returns the cached result if the source file hasn't changed.
#'
#' @param rast A SpatRaster (already renamed to short names).
#' @param vars Character vector of layer names to compute stats for.
#' @param source_path Path to the source raster file (used for fingerprint).
#' @param cache_dir Directory for the cache file. Default: `tempdir()`.
#'
#' @return A tibble with columns `variable`, `mean`, `sd`.
#' @export
compute_scale_df <- function(rast, vars, source_path,
                             cache_dir = file.path(paths$large,
                                                   "Data_processed/cache"),
                             force = FALSE) {
  
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # fingerprint: file size + modification time
  file_info <- file.info(source_path)
  fingerprint <- paste(file_info$size, file_info$mtime,
                       paste(sort(vars), collapse = ","))
  
  cache_path <- file.path(cache_dir, "climate_scale_df.rds")
  
  # check cache
  if (!force && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$fingerprint, fingerprint)) {
      cat("  Using cached scale_df\n")
      return(cached$scale_df)
    }
  }
  
  # compute fresh
  cat("  Computing scale_df from raster (this may take a few minutes)...\n")
  r_sub <- rast[[vars]]
  clim_means <- terra::global(r_sub, "mean", na.rm = TRUE)
  clim_sds <- terra::global(r_sub, "sd", na.rm = TRUE)
  
  scale_df <- tibble::tibble(
    variable = vars,
    mean = clim_means$mean,
    sd = clim_sds$sd
  )
  
  # save cache
  saveRDS(list(fingerprint = fingerprint, scale_df = scale_df), cache_path)
  cat("  Cached scale_df to:", cache_path, "\n")
  
  scale_df
}

# renaming variables ------------------------------------------------------

#' Convert between long and short climate variable names
#'
#' @param x Character vector of variable names to convert.
#' @param to_short Logical; if TRUE (default), convert from long (verbose
#'   raster names) to short names. If FALSE, convert from short to long.
#'
#' @return Character vector of converted names, same length as `x`.
#' @examples
#' climate_name_lookup("tmean_meanAnnAvg_CLIM")
#' climate_name_lookup(c("MAT", "MAP"), to_short = FALSE)
#' @export
climate_name_lookup <- function(x, to_short = TRUE) {
  long_to_short <- c(
    "tmean_meanAnnAvg_CLIM"                   = "MAT",
    "prcp_meanAnnTotal_CLIM"                  = "MAP",
    "PrecipTempCorr_meanAnnAvg_CLIM"          = "PrecipTempCorr",
    "T_warmestMonth_meanAnnAvg_CLIM"          = "T_warmestMonth",
    "T_coldestMonth_meanAnnAvg_CLIM"          = "T_coldestMonth",
    "precip_wettestMonth_meanAnnAvg_CLIM"     = "P_wettestMonth",
    "precip_driestMonth_meanAnnAvg_CLIM"      = "P_driestMonth",
    "precip_Seasonality_meanAnnAvg_CLIM"      = "P_seasonality",
    "aboveFreezing_month_meanAnnAvg_CLIM"     = "frost_free_months",
    "isothermality_meanAnnAvg_CLIM"           = "isothermality",
    "annWaterDeficit_meanAnnAvg_CLIM"         = "WD_mean",
    "annWetDegDays_meanAnnAvg_CLIM"           = "WDD_mean",
    "annVPD_mean_meanAnnAvg_CLIM"             = "VPD_mean",
    "annVPD_max_meanAnnAvg_CLIM"              = "VPD_max",
    "annVPD_min_meanAnnAvg_CLIM"              = "VPD_min",
    "annVPD_max_95percentile_CLIM"            = "VPD_max_p95",
    "annWaterDeficit_95percentile_CLIM"       = "WD_p95",
    "annWetDegDays_5percentile_CLIM"          = "WDD_p05",
    "durationFrostFreeDays_5percentile_CLIM"  = "frost_free_days_p05",
    "durationFrostFreeDays_meanAnnAvg_CLIM"   = "frost_free_days",
    "tmin_meanAnnAvg_CLIM"                    = "T_min",
    "tmax_meanAnnAvg_CLIM"                    = "T_max"
  )
  
  if (to_short) {
    lookup <- long_to_short
  } else {
    lookup <- setNames(names(long_to_short), long_to_short)
  }
  
  missing <- setdiff(x, names(lookup))
  if (length(missing) > 0) {
    stop("Names not found in climate lookup: ",
         paste(missing, collapse = ", "))
  }
  
  unname(lookup[x])
}



# misc --------------------------------------------------------------------

# useful for model's that don't allow zero values
replace_zero <- function(x) {
  stopifnot(all(x >=0))
  
  small <- min(x[x>0])
  
  if(any(x == 0) & small > 1 ) {
    warning('min non zero value > 1, consider other approach for replace 0s')
  }
  x[x == 0] <- small/2
  x
}
