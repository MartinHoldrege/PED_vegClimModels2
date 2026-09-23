# Helpers for the binomial (classification) cover models.
#
# Each @examples block assigns every argument, so you can run it and then step
# through the function body line by line.


#' Probability cutoff that turns predictions into classes
#'
#' Wraps `PresenceAbsence::optimal.thresholds()`. With method "PredPrev=Obs"
#' the cutoff is the one where the predicted fraction in the class equals the
#' observed fraction in the data used here. Other methods from that function
#' can be named instead.
#'
#' @param obs Observed class, 0 or 1.
#' @param pred Predicted probability, same length as `obs`.
#' @param method Method name, as in the `Method` column that
#'   `optimal.thresholds()` returns.
#' @return Numeric cutoff.
#' @examples
#' set.seed(1)
#' obs <- rbinom(500, size = 1, prob = 0.3)
#' pred <- plogis(rnorm(500, mean = -1 + obs))
#' method <- "PredPrev=Obs"
#' choose_threshold(obs = obs, pred = pred, method = method)
choose_threshold <- function(obs, pred, method) {
  stopifnot(length(obs) == length(pred),
            all(obs %in% c(0, 1)),
            all(pred >= 0 & pred <= 1))
  
  dat <- data.frame(ID = seq_along(obs), obs = obs, pred = pred)
  
  thresholds <- PresenceAbsence::optimal.thresholds(
    DATA = dat,
    threshold = 200,       # number of candidate cutoffs tested
    obs.prev = mean(obs),
    opt.methods = method
  )
  
  stopifnot(method %in% thresholds$Method)
  thresholds[thresholds$Method == method, 2]
}
