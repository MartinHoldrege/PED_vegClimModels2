# functions for creating quantile ('holdrege') type plots
# (as in Holdrege et al. 2024, Fire Ecology)

#' filter rows by climate variables
#' 
#' @description filters the dataframe by percentiles of the climate variable
#' columns. So the output includes rows corresponding the bottom 2 deciles and 
#' the top two deciles of each climate variable. Note this function has
#' been updated to filter by an arbitrary number of columns (not just climate)
#'
#' @param df dataframe that needs to have MAP, MAT, and prcpPropSum column
#' @param add_mid also add a seperate category of the center 2 deciles of 
#' each climate variable
#' @param filter_vars variables (usually climate variables), to split the others
#' into response and predicted
#'
#' @return dataframe with same columns as df but also filter_var,
#' percentile_category, which give the names of the climate variable filtered 
#' by and the percentile cut-off used that the given row fits in

filter_by_climate <- function(df, add_mid = FALSE,
                              filter_vars) {
  
  # percentile cuttoffs to use, keeping values below low, and above high
  low <- 0.2
  high <- 0.8
  
  
  names(filter_vars) <- filter_vars
  stopifnot(filter_vars %in% names(df))
  
  # fitting empirical cdf's so that percentiles of the climate variables
  # can be computed
  ecdf_list <- map(df[filter_vars], ecdf)
  
  # dataframe, where each column provides the percentiles of a climate variable
  # percentiles correspond to rows in df
  percentiles <- map_dfc(filter_vars, function(var) {
    ecdf_list[[var]](df[[var]])
  })
  
  # only keep rows falling in the low or high category, for each climate var
  df_filtered <- map_dfr(filter_vars, function(var) {
    df_low <- df[percentiles[[var]] < low, ]
    df_low$percentile_category <- paste0("<", low*100, "th")
    df_high <- df[percentiles[[var]] > high, ] 
    df_high$percentile_category <- paste0(">", high*100, "th")
    
    out <- bind_rows(df_low, df_high)
    
    out$percentile_category <- as.factor(out$percentile_category)
    
    # adding seperate category for the middle of the data
    if(add_mid){
      df_mid <- df[percentiles[[var]] < .6 & percentiles[[var]] > .4, ]
      df_mid$percentile_category <- "40th-60th"
      
      # create correct factor order
      levels <- levels(out$percentile_category)
      levels <- c(levels[1], "40th-60th", levels[2])
      # convert to character for binding
      out$percentile_category <- as.character(out$percentile_category)
      out <- bind_rows(out, df_mid)
      out$percentile_category <- factor(out$percentile_category,
                                        levels = levels)
      
    }
    out$filter_var <- var
    
    out
  })
  df_filtered$filter_var <- factor(df_filtered$filter_var, levels = filter_vars)
  df_filtered
}


#' Make long format dataframe with predictor variable becoming a column
#'
#' @param df dataframe (could be output from filter_by_climate)
#' @param response_vars names of response variables
#' @param pred_vars names of predictor variables
#' @param filter_var logical--whether this dataframe also includes 
#' filter_var and percentile_category columns that should be kept
#'
#' @return longform dataframe
predvars2long <- function(df, response_vars, 
                          pred_vars,
                          filter_var = FALSE
) {
  
  select_cols <- c(pred_vars, response_vars)
  
  if(filter_var) {
    select_cols <- c(select_cols, c("filter_var", "percentile_category"))
  }
  

  out <- df[, select_cols] %>% 
    # unname used here b/ https://github.com/tidyverse/tidyr/issues/1481
    pivot_longer(cols = all_of(unname(pred_vars)))
  
  out
}


#' convert longform dataframe to means of predictor quantiles
#'
#' @param df with columns of name (name of predictor variable), 'value'
#' (value of predictor variable), and 1 or more response variable columns.
#' The output of predvars2long() creates a correctly formatted df to use here
#' @param response_vars character vector, names of the response variables
#' @param filter_var logical--whether this dataframe also includes 
#' filter_var and percentile_category columns that should be kept
#' @param return_means logical. if false return the dataframe before
#' means have been calculated for each quantile
#' 
#' @return For each predictor variable calculate the mean of each decile
#' and the corresponding mean (of those same rows) of the response variable
longdf2deciles <- function(df, response_vars, filter_var = FALSE,
                           return_means = TRUE,
                           cut_points = seq(0, 1, 0.01),
                           summary_fun = c('mean', 'median')) {
  
  stopifnot(c("name", "value", response_vars) %in% names(df))
  
  summary_fun <- match.arg(summary_fun) 
  
  f <- switch(summary_fun,
              mean = mean,
              median = median)
  
  group_vars <- 'name'
  if(filter_var) {
    group_vars <- c(group_vars, c("filter_var", "percentile_category"))
  } 
  
  if(!filter_var & 'filter_var' %in% names(df)) {
    warning('dataframe includes a column named filter_var, the
            filter_var argument should probably be set to TRUE')
  }
  
  out0 <- df %>% 
    # the named vector in select was selecting by the names
    # not the vector values!
    select(all_of(group_vars), value, unname(response_vars)) %>% 
    group_by(across(all_of(group_vars))) %>% 
    nest() %>% 
    # empirical cdf
    # calculate the percentile of each data point based on the ecdf
    mutate(percentile = map(data, function(df) ecdf(df$value)(df$value))) %>% 
    # as_tibble() %>% 
    unnest(cols = c("data", "percentile")) %>% 
    group_by(across(all_of(group_vars))) %>% 
    mutate(decile = cut(percentile, cut_points,
                        labels = 1:(length(cut_points) - 1))) %>% 
    # calculate mean of response variables for each decile of each predictor
    # variable
    group_by(across(all_of(c(group_vars, 'decile'))), arrange = FALSE) 
    # lazy_dt() # this speeds up the code ~3x, add back in for speed but
  # requires dtplyr be loaded
  
  
  if(!return_means) {
    return(as_tibble(out0))
  }
  
    out <- out0 %>% 
      summarize(across(unname(response_vars), f),
                mean_value = f(value), # mean of predictor for that decile
                .groups = 'drop')

  
  as_tibble(out)
}

#' wide format data frame, to longformat dataframe summarized to quantiles
#'
#' @description wrapper around predvars2long and longdf2deciles
#'
#' @param df 
#' @param response_vars character vector of response vars
#' @param pre_vars character vector of response vars
#' @param filter_var whether to create a filter var (i.e. filter by the
#' climate variables, so group data into high and low percentiles
#' of each climate variable)
#' @param filter_vars vector of names of columns to use
#' as filter variables (only relevant if filter_var = TRUE)
#' @param add_mid whether to also filter by the middle percentiles of the 
#' climate vars
#' @param cut_points how to group/cut the percentiles before avging. By default
#' compute the average for each percentile
predvars2deciles <- function(df, response_vars, pred_vars,
                             filter_var = FALSE,
                             filter_vars = NULL,
                             add_mid = FALSE,
                             cut_points = seq(0, 1, 0.01),
                             summary_fun = 'mean') {
  
  if (filter_var && is.null(filter_vars)) {
    stop("filter_vars must be specified when filter_var = TRUE.")
  }
  
  stopifnot(
    is.logical(filter_var)
  )
  
  if (filter_var) {
    stopifnot(is.character(filter_vars))
    df <- filter_by_climate(df, add_mid = add_mid,
                            filter_vars = filter_vars)
  }
  
  # longformat df
  long_df <- predvars2long(df, response_vars = response_vars, 
                           pred_vars = pred_vars,
                           filter_var = filter_var)
  # mean of deciles
  out <- longdf2deciles(long_df, response_vars = response_vars,
                        filter_var = filter_var,
                        cut_points = cut_points,
                        summary_fun = summary_fun)
  out
}

#' calculate rmse for decile plot
#'
#' @param df dataframe with name column (containing names of predictor variables)
#' and columns corresponding to yvar and yvar_pred, where yvar is a fire
#' response variable
#' @param yvar string (e.g. mtbs_prop)
#'
#' @return dataframe, giving root mean squared error of quantile level
#' values
rmse4dotplot <- function(df, yvar) {
  df_list <- split(df, df$name) 
  observed <- yvar
  predicted <- paste0(yvar, "_pred")
  rmse_vector <- map_dbl(df_list, function(df){
    squared_error = (df[[observed]] - df[[predicted]])^2
    rmse <- sqrt(mean(squared_error))
  })
  rmse_vector
  
  # convert to dataframe for use in ggplot
  out <- tibble(name = names(rmse_vector),
                rmse = rmse_vector)
  out$rmse <- formatC(out$rmse, digits = 2, format = 'e')
  out
}


#' create dotplot of data summarized to deciles, faceted by predictor variable
#'
#' @param yvar name of the yvar to be plotted (string) (must be present in df) 
#' @param df dataframe longform with 'mean_value' column
#' (i.e. mean value of the predictor variable for the given decile) and 'name' 
#' column which gives the name of the predictor variable
#' @param method string, to be pasted into subtitle (method used to convert
#' polygons to rasters)
#' @param ylab string--y axis label
#' @param add_predicted logical, whether to also add model predicted data
#' to the plot. this requires the dataframe to 
decile_dotplot <- function(yvar, df, ylab = 'response',
                           xlab = "mean of quantile of predictor variable",
                           add_predicted = FALSE, title = NULL,
                           size = 0.75,
                           add_rmse = TRUE,
                           subtitle = NULL) {
  
  if('filter_var' %in% names(df)) {
    stop('filter_var column present, you should used decile_dotplot_filtered()')
  }
  
  caption <- "Each panel shows a different predictor variable"

  g <- ggplot(df, aes(x = .data[['mean_value']], y = .data[[yvar]])) +
    geom_point(aes(color = "Observed", shape = "Observed"),
               size = size) +
    facet_wrap(~name, scales = 'free_x') +
    labs(x = xlab,
         y = ylab,
         subtitle = subtitle,
         title = title) +
    theme(legend.position = 'top',
          legend.title = element_blank()) 
  g
  
  col_values <- c("Observed" = "black")
  shape_values <- c("Observed" = 19)
  
  # whether to also add dots for predicted probability
  # this requires predicted value columns to have the same name as yvar
  # but followed by _pred
  if(add_predicted) {
    yvar_pred <- paste0(yvar, "_pred")
    
    g <- g +
      geom_point(aes(y = .data[[yvar_pred]], color = 'Predicted',
                     shape = 'Predicted'), alpha = 0.5, size = size) 
    col_values <- c("Observed" = "black", "Predicted" = "blue")
    shape_values <- c(shape_values, "Predicted" = 17)
  }
  
  if (add_predicted & add_rmse) {
    rmse_df <- rmse4dotplot(df = df, yvar = yvar)
    
    if(is.factor(df$name)){
      rmse_df$name <- factor(rmse_df$name, levels = levels(df$name))
    }
    caption = paste0(caption, 
                     "\nRMSE of quantile averages shown in each panel")
    g <- g +
      geom_text(data = rmse_df, 
                aes(x = -Inf, y = Inf, label = rmse, hjust = -0.05,
                    vjust = 1.5), size = 3)
    g
  }
  
  out <- g +
    scale_color_manual(name = 'legend', values = col_values) +
    scale_shape_manual(name = 'legend', values = shape_values) +
    labs(caption = caption)
  out
}


#' Compute per-PFT quantile summaries across climate variables
#'
#' For each PFT, bins data by quantiles of each climate variable and computes
#' mean predicted biomass within each bin. Uses `predvars2deciles()` internally.
#'
#' @param fit Fitted cwexp model object (or any object accepted by
#'   `cwexp_mu_by_group`).
#' @param data Data frame with predictor and cover columns.
#' @param formula Model formula.
#' @param cover_cols Character vector of cover column names.
#' @param pred_vars Character vector of climate/predictor variables to bin by.
#' @param weighted Logical. If TRUE, returns cover-weighted per-PFT biomass.
#'   If FALSE, returns potential biomass (cover = 100%).
#' @param cut_points Numeric vector of quantile break points.
#'
#' @return A tibble in long format with columns: `PFT`, `name` (climate var),
#'   `mean_value` (mean or median climate within bin), `decile` (bin id), and `biomass`
#'   (mean predicted biomass within bin).
#' @export
pft_quantile_summary <- function(fit, data, formula, cover_cols, pred_vars,
                                 weighted = TRUE,
                                 cut_points = seq(0, 1, 0.01),
                                 summary_fun = 'median') {
  
  prep <- cwexp_prepare(data = data, formula = formula, cover_cols = cover_cols)
  
  mu_grp <- cwexp_mu_by_group(
    alpha = fit$par$alpha,
    B = fit$par$B,
    X = prep$X,
    C = prep$C,
    weighted = weighted
  )
  
  pft_labels <- stringr::str_remove(cover_cols, "Cov$")
  
  # add per-PFT predictions as columns
  colnames(mu_grp) <- stringr::str_remove(colnames(mu_grp), "Cov$")
  df <- cbind(data[, pred_vars, drop = FALSE], mu_grp)
  
  # predvars2deciles treats the PFT columns as response variables
  quant <- predvars2deciles(
    df = df,
    response_vars = pft_labels,
    pred_vars = pred_vars,
    cut_points = cut_points,
    summary_fun = summary_fun
  )
  
  # pivot PFTs to long format for facet_grid
  quant |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(pft_labels),
      names_to = "PFT",
      values_to = "biomass"
    ) |>
    dplyr::mutate(PFT = factor(PFT, levels = pft_labels))
}


#' Plot per-PFT biomass across climate quantiles
#'
#' Creates a PFT × climate grid plot from the output of
#' `pft_quantile_summary()`. Rows are PFTs, columns are climate variables.
#'
#' @param quant_df Output of `pft_quantile_summary()`.
#' @param title Character; plot title.
#' @param ylab Character; y-axis label.
#' @param point_color Character; color for the dots.
#' @param point_size Numeric; size of the dots.
#' @param scales Character; passed to `facet_grid()`. Default `"free"`.
#'
#' @return A ggplot object.
#' @export
plot_pft_quantiles <- function(quant_df,
                               title = NULL,
                               ylab = "Predicted biomass",
                               point_color = "blue",
                               point_size = 0.5,
                               scales = "free") {
  
  ggplot2::ggplot(quant_df, ggplot2::aes(x = mean_value, y = biomass)) +
    ggplot2::geom_point(size = point_size, color = point_color) +
    ggplot2::facet_grid(PFT ~ name, scales = scales, switch = "x") +
    ggplot2::labs(
      x = NULL,
      y = ylab,
      title = title
    ) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.placement = "outside",
      strip.background.x = ggplot2::element_blank()
    )
}
