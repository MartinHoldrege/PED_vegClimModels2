# main_cover.R
#
# Runs the cover modeling pipeline for selected models: fit, predict onto the
# CONUS rasters (current and future climate), and render the evaluation
# report. Data prep (Cover/DataPrep) is run separately.
#
# Each model is cover_specs[[cover_type]][[cover_model]][[vmc]]; the options
# are passed to the scripts on the command line (see params.R).
#
# Combined predictions (cover_type 'combined': combined tree cover, herbaceous
# cover by functional type) are not fit or predicted: they combine the
# prediction rasters of their component models (04_combine_tree.R,
# 04_combine_herb.R), so those must already exist. Combined runs are done
# after all other runs.
#
# September, 2026

library(dplyr)
library(purrr)

# what to run -------------------------------------------------------------
# Wait before starting (another job is still running)
wait_hours <- 0


run_fit      <- FALSE
run_predict  <- TRUE
run_evaluate <- TRUE

# cover data version(s): vc_default, unless the model is listed in vc_model
# (several versions = each vmc is run on each; not
# used by combined runs)
vc_default <- "c02"
vc_model <- list(
  bare_ground = c("c01", "c02"),
  shrub = c("c01", 'c03'),
  herbaceous = "c01", # "c03",
  needle = "c01", c3 = "c01", c4 = "c01", forb = "c01")

fit_scripts <- c(
  classification = "Analysis/Cover/Fit/02_fit_classification.R",
  cover          = "Analysis/Cover/Fit/02_fit_cover.R",
  proportion     = "Analysis/Cover/Fit/02_fit_cover.R"
)

# model versions to run, by family and model
cover_runs <- list(
  classification = list(
    # forest    = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    # zero_tree = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1')
    ),
  cover          = list(
    # tree_forest    = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    # tree_nonforest = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    # shrub = c('m01.8')
    # herbaceous = c('m01.5', 'm04.1'),
    # m05.5: f5 plus anomalies
    herbaceous = 'm05.5'
    # bare_ground = c('m01.5', 'm04.1')
  ),
  proportion = list(
    c3 = 'm05.5',
    c4 = 'm05.5',
    forb = 'm05.5'
    # needle = c('m01.5', 'm04.1')
    ),
  # versions of cover_specs$combined
  combined       = list(
    # tree = c('m01.4', 'm01.5', 'm01.6', 'm01.7', 'm04.1'),
    herbaceous = 'm05.5',
    all = 'm05.5'
  )
)

# scripts ------------------------------------------------------------------

# fitting script per family
fit_scripts <- c(
  classification = "Analysis/Cover/Fit/02_fit_classification.R",
  cover          = "Analysis/Cover/Fit/02_fit_cover.R",
  proportion     = "Analysis/Cover/Fit/02_fit_cover.R"
)

predict_script <- "Analysis/Cover/Fit/03_predict_raster.R"

# combining script per combined model (cover_specs$combined; tree_all isn't
# written to disk, so has no script)
combine_scripts <- c(
  tree       = "Analysis/Cover/Fit/04_combine_tree.R",
  herbaceous = "Analysis/Cover/Fit/04_combine_herb.R"
)

# evaluation report per family 
eval_rmds <- c(
  classification = "Analysis/Cover/Evaluate/01_classification_diagnostics.Rmd",
  cover          = "Analysis/Cover/Evaluate/01_cover_diagnostics.Rmd",
  proportion     = "Analysis/Cover/Evaluate/01_cover_diagnostics.Rmd"
)

# evaluation report per combined model (none for herbaceous yet)
combined_rmds <- c(
  tree = "Analysis/Cover/Evaluate/01_tree_combined_diagnostics.Rmd",
  all  = "Analysis/Cover/Evaluate/01_combined_all_diagnostics.Rmd"
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
  mutate(vc = map(cover_model, \(m) vc_model[[m]] %||% vc_default)) |>
  tidyr::unnest(vc) |>
  # combined runs use the outputs of the others, so they go last
  arrange(cover_type == "combined")

stopifnot(all(runs$cover_type %in% c(names(fit_scripts), "combined")),
          with(runs, all(cover_model[cover_type == "combined"] %in%
                           c(names(combine_scripts), names(combined_rmds)))))
print(runs)


# delay -------------------------------------------------------------------


message("Sleeping ", wait_hours, " h; started ", format(Sys.time()))
Sys.sleep(wait_hours * 3600)
message("Resuming at ", format(Sys.time()))
# misc other scripts to run (other scripts to run) ------------------------------------------

if(FALSE) {
# temporary: build the c03 training data
  dataprep <- file.path("Analysis/Cover/DataPrep",
                        c(# "06_cover_add-rap.R",
                          # "07_summarise_daymet_climate_points.R",
                          # "08_combine_cover_and_covariates.R",
                          "09_compute_scale_params.R"))
  stopifnot(file.exists(dataprep))
  for (s in dataprep) {
    message('running: ', s)
    callr::rscript(s, cmdargs = "--vc=c03")
  }
}
# run ------------------------------------------------------------------------

pwalk(runs, \(cover_type, cover_model, vmc, vc) {

  is_combined <- cover_type == "combined"
  
  cmdargs <- c(paste0("--vc=", vc),
               paste0("--cover_type=", cover_type),
               paste0("--cover_model=", cover_model),
               paste0("--vmc=", vmc))
  
  if (run_fit && !is_combined) {
    callr::rscript(fit_scripts[[cover_type]], cmdargs = cmdargs)
  }
  
  # for combined runs, "predict" is combining the component predictions; NA
  # if there is nothing to write (e.g. combined all, computed in its report)
  script <- if (is_combined) combine_scripts[cover_model] else predict_script
  if (run_predict && !is.na(script)) {
    callr::rscript(unname(script), cmdargs = cmdargs)
  }
  
  # NA if there is no report (e.g. combined herbaceous)
  rmd <- if (is_combined) combined_rmds[cover_model] else eval_rmds[cover_type]
  
  if (run_evaluate && !is.na(rmd)) {
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
      unname(rmd),
      params = params,
      output_file = output_file,
      output_dir = out_dir,
      knit_root_dir = getwd(),  # the Rmd sources Functions/init.R
      envir = new.env()
    )
  }
})
