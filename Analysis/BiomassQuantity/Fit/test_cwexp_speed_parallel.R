# Purpose:
# Benchmark run_inner_cv: sequential vs parallel execution.
# Helps assess whether parallelization is worthwhile at different data sizes.
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

suppressWarnings(source('Functions/init.R'))
source_functions()
library(future)
library(future.apply)

# output ------------------------------------------------------------------

out_path <- "Data_processed/BiomassQuantityData/speed_parallel_results.txt"
sink(out_path, split = TRUE)

cat("run_inner_cv parallel vs sequential benchmark\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Platform:", sessionInfo()$platform, "\n")
cat("R version:", R.version.string, "\n")
cat("Physical cores:", parallel::detectCores(logical = FALSE), "\n")
cat("Logical cores:", parallel::detectCores(), "\n")
cat(paste0(rep("-", 60), collapse = ""), "\n\n")

# setup -------------------------------------------------------------------

dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)
dll_path <- "src/cwexp_lognormal_en_tmb.cpp"

n_folds <- 3
n_lambda <- 5
en_alpha <- 0.5
formula <- totalBio ~ tmean + ppt + vpd + sand

sizes <- c(1e3, 5e3, 1e4)

set.seed(1)

results <- vector("list", length(sizes))

workers <- min(n_folds, parallel::detectCores(logical = FALSE)) 

for (s in seq_along(sizes)) {
  n <- sizes[s]
  cat("=== n =", format(n, big.mark = ","), "===\n")

  dummy <- cwexp_make_dummy_data(n = n)
  dat <- dummy$data
  cover_cols <- dummy$spec$cover_cols

  clust <- make_env_clusters(
    data = dat,
    vars = c("tmean", "ppt", "vpd", "sand"),
    k = 15,
    seed = 1
  )
  dat$env_cluster <- clust$env_cluster

  folds <- make_cluster_folds(
    env_cluster = dat$env_cluster,
    n_folds = n_folds,
    seed = 1
  )

  shared_args <- list(
    data = dat,
    folds = folds,
    en_alpha = en_alpha,
    lambda_max_fun = cwexp_lambda_max_l1_tmb,
    lambda_max_args = list(
      formula = formula,
      cover_cols = cover_cols,
      dll_en = dll_en
    ),
    lambda_path_args = list(
      n_lambda = n_lambda,
      include_zero = TRUE
    ),
    fit_path_fun = cwexp_fit_lambda_path_tmb,
    fit_path_args = list(
      formula = formula,
      cover_cols = cover_cols,
      dll = dll_en,
      include_report = FALSE
    )
  )

  # sequential
  t_seq <- system.time({
    res_seq <- do.call(run_inner_cv, c(shared_args, list(parallel = FALSE)))
  })

  # parallel
  plan(multisession, workers = workers)
  t_par <- system.time({
    res_par <- do.call(run_inner_cv, c(shared_args, list(
      parallel = TRUE,
      dll_path = dll_path
    )))
  })
  plan(sequential)

  # check results match
  lambdas_match <- identical(res_seq$selected$lambda, res_par$selected$lambda)

  results[[s]] <- tibble::tibble(
    n = n,
    seconds_seq = t_seq["elapsed"],
    seconds_par = t_par["elapsed"],
    speedup = t_seq["elapsed"] / t_par["elapsed"],
    lambdas_match = lambdas_match
  )

  cat("  Sequential:", round(t_seq["elapsed"], 2), "s\n")
  cat("  Parallel:  ", round(t_par["elapsed"], 2), "s\n")
  cat("  Speedup:   ", round(t_seq["elapsed"] / t_par["elapsed"], 2), "x\n")
  cat("  Lambdas match:", lambdas_match, "\n\n")
}

# summarize ---------------------------------------------------------------

bench <- dplyr::bind_rows(results)
print(bench)

sink()
cat("Results written to:", out_path, "\n")
