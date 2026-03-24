#' Plot CV scores across the lambda path
#'
#' Creates a faceted ggplot showing one or more CV metrics vs lambda,
#' with error bars, selected lambda, and (for loss metrics) the 1-SE
#' threshold.
#'
#' @param score_summary Data frame from `summarize_scores()`, with columns
#'   for `lambda`, metrics, and `{metric}_se`.
#' @param selected_lambda Numeric; the lambda chosen by the selection rule.
#' @param metrics Named character vector where names are column names in
#'   `score_summary` and values are display labels. Default covers
#'   mae_log, rmse_log, mae, and cor.
#' @param higher_is_better Character vector of metric names where larger
#'   values are better (e.g., correlation). The 1-SE threshold line is
#'   only drawn for metrics not in this set.
#' @param x_scale Character; `"raw"` (default), `"log"`, or `"rank"`.
#'   `"rank"` spaces lambdas evenly and labels the axis with actual lambda
#'   values, which helps when small lambdas are clustered together.
#'
#' @return A ggplot object.
#' @examples
#' # using output from run_inner_cv
#' # plot_cv_scores(cv$score_summary, cv$selected$lambda)
#' # evenly spaced:
#' # plot_cv_scores(cv$score_summary, cv$selected$lambda, x_scale = "rank")
#' @export
plot_cv_scores <- function(score_summary,
                           selected_lambda,
                           metrics = c(
                             mae_log  = "MAE (log scale)",
                             rmse_log = "RMSE (log scale)",
                             mae      = "MAE",
                             rmse = 'RMSE',
                             cor      = "Correlation"
                           ),
                           higher_is_better = "cor",
                           x_scale = c("raw", "log", "rank")) {
  
  x_scale <- match.arg(x_scale)
  metric_names <- names(metrics)
  stopifnot(all(metric_names %in% names(score_summary)))
  score_summary$lambda <- signif(score_summary$lambda, 5)
  # --- x-axis mapping ---
  lambdas_sorted <- sort(unique(score_summary$lambda))
  lambda_to_rank <- setNames(seq_along(lambdas_sorted), as.character(lambdas_sorted))
  selected_x <- switch(
    x_scale,
    raw  = selected_lambda,
    log  = selected_lambda,
    rank = unname(lambda_to_rank[as.character(selected_lambda)])
  )
  
  # build long data for faceting
  plot_df <- purrr::imap_dfr(metrics, function(label, m) {
    se_col <- paste0(m, "_se")
    dplyr::tibble(
      lambda = score_summary$lambda,
      x_val  = switch(
        x_scale,
        raw  = score_summary$lambda,
        log  = score_summary$lambda,
        rank = unname(lambda_to_rank[as.character(score_summary$lambda)])
      ),
      value  = score_summary[[m]],
      se     = if (se_col %in% names(score_summary)) score_summary[[se_col]] else 0,
      metric = label
    )
  })
  
  # preserve facet order
  plot_df$metric <- factor(plot_df$metric, levels = unname(metrics))
  
  # 1-SE threshold lines (for loss metrics only)
  threshold_df <- purrr::imap_dfr(metrics, function(label, m) {
    if (m %in% higher_is_better) return(NULL)
    se_col <- paste0(m, "_se")
    y <- score_summary[[m]]
    se <- if (se_col %in% names(score_summary)) score_summary[[se_col]] else 0
    best_idx <- which.min(y)
    dplyr::tibble(
      metric = label,
      threshold = y[best_idx] + se[best_idx]
    )
  })
  if (nrow(threshold_df) > 0) {
    threshold_df$metric <- factor(threshold_df$metric, levels = levels(plot_df$metric))
  }
  
  g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x_val, y = value)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = value - se, ymax = value + se),
      width = 0
    ) +
    ggplot2::geom_vline(xintercept = selected_x,
                        linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(y = NULL)
  
  # x-axis formatting
  if (x_scale == "log") {
    g <- g +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "Lambda (log scale)")
  } else if (x_scale == "rank") {
    # label with actual lambda values, rounded for readability
    label_df <- dplyr::tibble(
      rank = unname(lambda_to_rank),
      label = signif(lambdas_sorted, 2)
    )
    g <- g +
      ggplot2::scale_x_continuous(
        breaks = label_df$rank,
        labels = label_df$label
      ) +
      ggplot2::labs(x = "Lambda (evenly spaced)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = 7))
  } else {
    g <- g + ggplot2::labs(x = "Lambda")
  }
  
  if (nrow(threshold_df) > 0) {
    g <- g +
      ggplot2::geom_hline(
        data = threshold_df,
        ggplot2::aes(yintercept = threshold),
        linetype = "dotted", color = "blue"
      )
  }
  
  g
}
