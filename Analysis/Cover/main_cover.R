# main_cover.R
#
# Runs the cover modeling pipeline for selected models: fit, predict onto the
# CONUS rasters (current and future climate), and render the evaluation
# report. Data prep (Cover/DataPrep) is run separately.
#
# Each model is cover_specs[[cover_type]][[cover_model]][[vmc]]; the options
# are passed to the scripts on the command line (see params.R).
#
# Combined predictions (cover_type 'combined', e.g. combined tree cover) are
# not fit or predicted: they combine the prediction rasters of their component
# models (04_combine_tree.R), so those must already exist. Combined runs are
# done after all other runs.
#
# September, 2026

library(dplyr)
library(purrr)

# what to run -------------------------------------------------------------

run_fit      <- TRUE
run_predict  <- TRUE
run_evaluate <- TRUE

vc <- "c02"  # cover data version (not used by combined runs)

# model versions to run, by family and model
cover_runs <- list(
  classification = list(
    forest    = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    zero_tree = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1')
    ),
  cover          = list(
    tree_forest    = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    tree_nonforest = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1')
  ),
  proportion     = list(),
  # versions of cover_specs$combined
  combined       = list(
    tree = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1')
  )
)

# scripts ------------------------------------------------------------------

# fitting script per family
fit_scripts <- c(
  classification = "Analysis/Cover/Fit/02_fit_classification.R",
  cover          = "Analysis/Cover/Fit/02_fit_cover.R"
)

predict_script <- "Analysis/Cover/Fit/03_predict_raster.R"

# combined predictions (only tree so far)
combine_script <- "Analysis/Cover/Fit/04_combine_tree.R"

# evaluation report per family 
eval_rmds <- c(
  classification = "Analysis/Cover/Evaluate/01_classification_diagnostics.Rmd",
  cover          = "Analysis/Cover/Evaluate/01_cover_diagnostics.Rmd",
  combined       = "Analysis/Cover/Evaluate/01_tree_combined_diagnostics.Rmd"
)

report_dirs <- c(classification = "Reports/Cover/classification",
                 cover          = "Reports/Cover/cover",
                 proportion     = "Reports/Cover/proportion",
                 combined       = "Reports/Cover/combined")

# one row per model to run --------------------------------------------------

runs <- imap(cover_runs, \(models, cover_type) {
  imap(models, \(vmcs, cover_model) {
    tibble(cover_type = cover_type, cover_model = cover_model, vmc = vmcs)
  }) |>
    bind_rows()
}) |>
  bind_rows() |>
  # combined runs use the outputs of the others, so they go last
  arrange(cover_type == "combined")

stopifnot(all(runs$cover_type %in% c(names(fit_scripts), "combined")))
print(runs)


# testing (other scripts to run) ------------------------------------------

# callr::rscript("Analysis/Cover/Fit/test_mars.R")

# run ------------------------------------------------------------------------

pwalk(runs, \(cover_type, cover_model, vmc) {

  is_combined <- cover_type == "combined"
  
  cmdargs <- c(paste0("--vc=", vc),
               paste0("--cover_type=", cover_type),
               paste0("--cover_model=", cover_model),
               paste0("--vmc=", vmc))
  
  if (run_fit && !is_combined) {
    callr::rscript(fit_scripts[[cover_type]], cmdargs = cmdargs)
  }
  
  # for combined runs, "predict" is combining the component predictions
  if (run_predict) {
    script <- if (is_combined) combine_script else predict_script
    callr::rscript(script, cmdargs = cmdargs)
  }
  
  if (run_evaluate && cover_type %in% names(eval_rmds)) {
    out_dir <- report_dirs[[cover_type]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # combined reports have no data version (components can differ)
    if (is_combined) {
      params <- list(vmc = vmc)
      output_file <- paste0(cover_type, "_", cover_model, "_", vmc, ".html")
    } else {
      params <- list(cover_type = cover_type, cover_model = cover_model,
                     vc = vc, vmc = vmc)
      output_file <- paste0(cover_type, "_", cover_model, "_", vc, "-", vmc,
                            ".html")
    }
    
    rmarkdown::render(
      eval_rmds[[cover_type]],
      params = params,
      output_file = output_file,
      output_dir = out_dir,
      knit_root_dir = getwd(),  # the Rmd sources Functions/init.R
      envir = new.env()
    )
  }
})
