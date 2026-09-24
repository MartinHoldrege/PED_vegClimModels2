# main_cover.R
#
# Runs the cover modeling pipeline for selected models: fit, predict onto the
# CONUS rasters (current and future climate), and render the evaluation
# report. Data prep (Cover/DataPrep) is run separately.
#
# Each model is cover_specs[[cover_type]][[cover_model]][[vmc]]; the options
# are passed to the scripts on the command line (see params.R).
#
# September, 2026

library(dplyr)
library(purrr)

# what to run -------------------------------------------------------------

run_fit      <- TRUE
run_predict  <- TRUE
run_evaluate <- TRUE

vc <- "c01"  # cover data version

# model versions to run, by family and model
cover_runs <- list(
  classification = list(forest = "m04"),
  cover          = list(),
  proportion     = list()
)

# scripts ------------------------------------------------------------------

# fitting script per family (only classification exists so far)
fit_scripts <- c(
  classification = "Analysis/Cover/Fit/02_fit_classification.R"
)

predict_script <- "Analysis/Cover/Fit/03_predict_raster.R"

# evaluation report per family
eval_rmds <- c(
  classification = "Analysis/Cover/Evaluate/01_classification_diagnostics.Rmd"
)

report_dirs <- c(classification = "Reports/Cover/classification",
                 cover          = "Reports/Cover/cover",
                 proportion     = "Reports/Cover/proportion")

# one row per model to run --------------------------------------------------

runs <- imap(cover_runs, \(models, cover_type) {
  imap(models, \(vmcs, cover_model) {
    tibble(cover_type = cover_type, cover_model = cover_model, vmc = vmcs)
  }) |>
    bind_rows()
}) |>
  bind_rows()

stopifnot(all(runs$cover_type %in% names(fit_scripts)),
          all(runs$cover_type %in% names(eval_rmds)))
print(runs)

# run ------------------------------------------------------------------------

pwalk(runs, \(cover_type, cover_model, vmc) {
  
  cmdargs <- c(paste0("--vc=", vc),
               paste0("--cover_type=", cover_type),
               paste0("--cover_model=", cover_model),
               paste0("--vmc=", vmc))
  
  if (run_fit) {
    callr::rscript(fit_scripts[[cover_type]], cmdargs = cmdargs)
  }
  
  if (run_predict) {
    callr::rscript(predict_script, cmdargs = cmdargs)
  }
  
  if (run_evaluate) {
    out_dir <- report_dirs[[cover_type]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    rmarkdown::render(
      eval_rmds[[cover_type]],
      params = list(cover_type = cover_type, cover_model = cover_model,
                    vc = vc, vmc = vmc),
      output_file = paste0(cover_type, "_", cover_model, "_", vc, "-", vmc,
                           ".html"),
      output_dir = out_dir,
      knit_root_dir = getwd(),  # the Rmd sources Functions/init.R
      envir = new.env()
    )
  }
})
