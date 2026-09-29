# test_mars_glmnet.R
#
# Overnight test: MARS proposes hinge terms, the lasso (glmnet) selects and
# shrinks them. Main effects only (no squares or logs), hinge x hinge
# interactions allowed (degree = 2).
#
# Cross-validation, one round, same rigor as m01 (not nested):
#   for each environmental fold: MARS on the other folds -> candidate hinge
#   columns -> glmnet on those columns over a fixed lambda grid -> predict the
#   held-out fold at every lambda.
#   lambda: lowest mean out-of-fold deviance over all folds (as m01).
#   final model: MARS on all data, glmnet at that lambda.
# Compared with m01 and m04 (random forest), out of fold, nationally and in
# Nevada and the Sandhills.
#
# Run from the repo root (m01 and m04 must be fit). Runs as a job: prints
# progress with times, and saves
#   test_mars_glmnet_<vc>_progress.rds  after each fold (in case of a crash)
#   test_mars_glmnet_<vc>.rds           results (see "save" at the end)
#   test_mars_glmnet_<vc>.txt           the printed summary


# setup -------------------------------------------------------------------

source("Functions/init.R")
source_functions()
source("Functions/models/cover_specs.R")
library(earth)
library(glmnet)

spec <- cover_specs$classification$forest$m01  # rows, response, folds
test_run <- FALSE  # TRUE: 20,000 random rows, for a quick check

vars <- c("MAT", "MAP", "PrecipTempCorr", "isothermality", "WD_p95",
          "clay_surface", "awc")

# MARS forward pass: nk caps the number of candidate terms; thresh stops it
# early once a new term improves the fit by less than this (earth's default
# 0.001 would stop well short of nk). The lasso does the selecting, so be
# generous with candidates
nk <- 81
thresh <- 1e-5

# one lambda grid for every fold, so errors can be averaged across folds
lambda_grid <- 10^seq(-1, -5, length.out = 50)

out_stem <- paste0("test_mars_glmnet_", opt$vc)
log_time <- \(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# data, as in 02_fit_classification.R -------------------------------------

dat <- read_cover_training(opt$vc)
source_vars <- unique(c(vars, spec$cv$cluster_vars))
dat <- dat[complete.cases(dat[c(spec$response, source_vars)]), ]
if (test_run) dat <- slice_sample(dat, n = 20000)

obs <- as.integer(dat[[spec$response]] > spec$cover_threshold)
x <- as.data.frame(dat[vars])  # original units, so knots are readable

foldid <- make_env_clusters(dat,
                            vars = spec$cv$cluster_vars,
                            iter.max = 500,
                            k = spec$cv$k_clusters,
                            seed = 1,
                            group = dat$cell)$env_cluster
folds <- sort(unique(foldid))
log_time("data:", nrow(dat), "rows,", length(folds), "folds")

# cross-validation --------------------------------------------------------

# out-of-fold predicted probability, one column per lambda
oof <- matrix(NA_real_, nrow(dat), length(lambda_grid))
fold_info <- NULL

for (f in folds) {
  t0 <- Sys.time()
  held_out <- foldid == f
  
  # candidate hinge columns from the training folds only (pmethod = "none":
  # keep every forward-pass term; the lasso does the selecting)
  mars_f <- earth(x[!held_out, ], obs[!held_out], degree = 2, nk = nk,
                  thresh = thresh, pmethod = "none")
  bx_train <- model.matrix(mars_f)[, -1, drop = FALSE]  # drop the intercept
  bx_held <- model.matrix(mars_f, x = x[held_out, ])[, -1, drop = FALSE]
  
  fit_f <- glmnet(bx_train, obs[!held_out], family = "binomial",
                  lambda = lambda_grid)
  oof[held_out, ] <- predict(fit_f, bx_held, type = "response")
  
  fold_info <- bind_rows(fold_info, tibble(
    fold = f,
    n_held_out = sum(held_out),
    n_candidates = ncol(bx_train),
    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  saveRDS(list(oof = oof, fold_info = fold_info),
          paste0(out_stem, "_progress.rds"))
  log_time("fold", f, "done:", ncol(bx_train), "candidates")
}

# lambda: lowest mean out-of-fold deviance --------------------------------

deviance <- function(obs, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -2 * mean(obs * log(p) + (1 - obs) * log(1 - p))
}
auc <- function(obs, p) {
  PresenceAbsence::auc(data.frame(id = seq_along(obs), obs = obs, pred = p),
                       st.dev = FALSE)
}

cv_table <- tibble(lambda = lambda_grid,
                   deviance_oof = map_dbl(seq_along(lambda_grid),
                                          \(j) deviance(obs, oof[, j])))
i_min <- which.min(cv_table$deviance_oof)
lambda_min <- lambda_grid[i_min]
if (i_min %in% c(1, length(lambda_grid))) {
  warning("chosen lambda is at the end of lambda_grid; widen the grid")
}
log_time("chosen lambda:", signif(lambda_min, 3))

# final model -------------------------------------------------------------

mars_all <- earth(x, obs, degree = 2, nk = nk, thresh = thresh,
                  pmethod = "none")
bx_all <- model.matrix(mars_all)[, -1, drop = FALSE]
fit_all <- glmnet(bx_all, obs, family = "binomial", lambda = lambda_grid)

b <- coef(fit_all, s = lambda_min)
terms <- tibble(term = rownames(b), coef = as.vector(b)) |>
  filter(coef != 0)

pred_in <- as.vector(predict(fit_all, bx_all, s = lambda_min,
                             type = "response"))
threshold <- choose_threshold(obs = obs, pred = pred_in,
                              method = spec$threshold_method)
log_time("final model:", nrow(terms) - 1, "terms kept of",
         ncol(bx_all), "candidates")

# compare with m01 and m04, nationally and by region ----------------------

fit_m01 <- read_cover_model("classification", "forest", opt$vc, "m01")
fit_m04 <- read_cover_model("classification", "forest", opt$vc, "m04")

comp <- dat |>
  select(cell, year, x, y, clay_surface) |>
  mutate(observed = obs, p_hybrid = oof[, i_min]) |>
  inner_join(select(fit_m01$predictions, cell, year, p_m01 = pred_oof),
             by = c("cell", "year")) |>
  inner_join(select(fit_m04$predictions, cell, year, p_m04 = pred_oof),
             by = c("cell", "year"))
stopifnot(nrow(comp) == nrow(dat))  # same rows, same order as pred_in

# regions: Nevada, and the Nebraska Sandhills approximated as Nebraska pixels
# with very sandy soils (clay < 5%)
states <- get_or_cache_states(crs = terra::crs(read_mask()))
state <- sf::st_as_sf(comp, coords = c("x", "y"), crs = sf::st_crs(states)) |>
  sf::st_join(states["NAME"]) |>
  pull(NAME)
comp <- comp |>
  mutate(region = case_when(
    state == "Nevada" ~ "Nevada",
    state == "Nebraska" & clay_surface < 5 ~ "Sandhills (NE, clay < 5%)",
    .default = "rest"))

scores <- function(d) {
  tibble(model = c("m01", "m04 (forest)", "MARS + lasso"),
         accuracy = c(mean((d$p_m01 > fit_m01$threshold) == d$observed),
                      mean((d$p_m04 > fit_m04$threshold) == d$observed),
                      mean((d$p_hybrid > threshold) == d$observed)),
         auc = c(auc(d$observed, d$p_m01), auc(d$observed, d$p_m04),
                 auc(d$observed, d$p_hybrid)))
}
by_region <- map(split(comp, comp$region),
                 \(d) mutate(scores(d), region = d$region[1], n = nrow(d))) |>
  list_rbind()
results <- bind_rows(mutate(scores(comp), region = "all", n = nrow(comp)),
                     by_region) |>
  select(region, n, model, accuracy, auc)

# summary -----------------------------------------------------------------

# terms are earth's names: h(x - k) means max(0, x - k), and a * b is a
# product (interaction). Coefficients are on the logit scale, for the
# columns as they are (glmnet reports them unstandardized)
summary_text <- c(
  paste("MARS + lasso test, data", opt$vc, "-", format(Sys.time())),
  paste("rows:", nrow(dat), "| nk:", nk, "| thresh:", thresh,
        "| lambda:", signif(lambda_min, 3),
        "| terms kept:", nrow(terms) - 1, "of", ncol(bx_all)),
  "", "Per fold:", capture.output(print(fold_info, n = Inf)),
  "", "Out-of-fold comparison:", capture.output(print(results, n = Inf)),
  "", "Final model, logit(p) = sum of coef * term:",
  capture.output(print(terms, n = Inf))
)
writeLines(summary_text, paste0(out_stem, ".txt"))
cat(summary_text, sep = "\n")

# save --------------------------------------------------------------------

# small enough to keep: no model objects, which store copies of the data
saveRDS(list(
  results = results,          # accuracy and AUC by model and region
  terms = terms,              # final model: term (earth notation), coef
  cv_table = cv_table,        # out-of-fold deviance per lambda
  lambda_min = lambda_min,
  threshold = threshold,      # cutoff on predicted probability
  fold_info = fold_info,      # per fold: candidates, minutes
  predictions = select(comp, cell, year, region, observed, p_m01, p_m04,
                       p_hybrid_oof = p_hybrid) |>
    mutate(p_hybrid_in = pred_in)
), paste0(out_stem, ".rds"))
log_time("saved", paste0(out_stem, ".rds"), "and .txt")