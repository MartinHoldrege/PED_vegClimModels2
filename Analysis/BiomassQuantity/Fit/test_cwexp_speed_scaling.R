# Purpose:
# Benchmark cwexp_fit_tmb fitting time as a function of dataset size.
# Helps estimate feasibility of fitting to large (10^5–10^6) datasets.
#

# dependencies ------------------------------------------------------------

source("Functions/data/simulate_data.R")
source("Functions/models/cwexp_helpers.R")
source("Functions/models/cwexp_lognormal_tmb.R")

# output ------------------------------------------------------------------

out_path <- "Data_processed/BiomassQuantity/Fit/speed_scaling_results.txt"
sink(out_path, split = TRUE)  # split = TRUE also prints to console

cat("cwexp_fit_tmb speed scaling benchmark\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Platform:", sessionInfo()$platform, "\n")
cat("R version:", R.version.string, "\n")
cat(paste0(rep("-", 60), collapse = ""), "\n\n")

# setup -------------------------------------------------------------------

dll_en <- cwexp_tmb_compile("src/cwexp_lognormal_en_tmb.cpp", quiet = TRUE)

sizes <- c(1e3, 1e4, 1e5, 5e5)
formula <- totalBio ~ tmean + ppt + vpd + sand

set.seed(1)

results <- vector("list", length(sizes))

for (i in seq_along(sizes)) {
  n <- sizes[i]
  cat("--- n =", format(n, big.mark = ","), "---\n")
  
  dummy <- cwexp_make_dummy_data(n = n)
  dat <- dummy$data
  cover_cols <- dummy$spec$cover_cols
  
  # time a single elastic-net fit (lambda = 0, so unpenalized but using EN DLL)
  t_unpen <- system.time({
    fit0 <- cwexp_fit_tmb(
      data = dat,
      formula = formula,
      cover_cols = cover_cols,
      dll = dll_en,
      penalty = "elastic_net",
      en_alpha = 0.5,
      lambda = 0,
      include_report = FALSE
    )
  })
  
  # time a penalized fit (moderate lambda)
  t_pen <- system.time({
    fit_pen <- cwexp_fit_tmb(
      data = dat,
      formula = formula,
      cover_cols = cover_cols,
      dll = dll_en,
      penalty = "elastic_net",
      en_alpha = 0.5,
      lambda = 0.01,
      include_report = FALSE
    )
  })
  
  results[[i]] <- tibble::tibble(
    n = n,
    seconds_unpen = t_unpen["elapsed"],
    seconds_pen = t_pen["elapsed"],
    convergence_unpen = fit0$tmb$opt$convergence,
    convergence_pen = fit_pen$tmb$opt$convergence,
    iterations_unpen = fit0$tmb$opt$iterations,
    iterations_pen = fit_pen$tmb$opt$iterations
  )
  
  cat("  unpenalized:", round(t_unpen["elapsed"], 2), "s",
      "(", fit0$tmb$opt$iterations, "iter )\n")
  cat("  penalized:  ", round(t_pen["elapsed"], 2), "s",
      "(", fit_pen$tmb$opt$iterations, "iter )\n\n")
}

# summarize ---------------------------------------------------------------

bench <- dplyr::bind_rows(results)
print(bench)

cat("\nScaling relative to smallest dataset:\n")
bench$relative_unpen <- bench$seconds_unpen / bench$seconds_unpen[1]
bench$relative_pen <- bench$seconds_pen / bench$seconds_pen[1]
print(bench[, c("n", "seconds_unpen", "relative_unpen",
                "seconds_pen", "relative_pen")])

sink()
cat("Results written to:", out_path, "\n")