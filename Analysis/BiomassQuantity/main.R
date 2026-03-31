# For running other scripts and rendering .rmd docs

source('Functions/constants.R') # for pfts
source('Functions/models/model_specs.R')
source('Functions/general.R')
# params ------------------------------------------------------------------
run_sim <- TRUE # simulate data
run_fit_model <- TRUE
run_model_diagnostics <- TRUE
run_model_diagnostics_sim <- TRUE


# [vs]: data version (s for simulated, d- for real data), 
#       
# p for purer version (sampling purer pixels)
#     p01 more appropriate for the simulated data, (p02 for real data)
# m for model version
#     m01 & m03 matches the formula for simulated data

suffix <- 'd04-p02-m07' # 's07-p01-m03' # 's07-p01-m08'# 'd03-p02-m06' # 's06-p01-m03'

# for rmd's
output_dir <- "Reports/BiomassQuantity/"
knit_root_dir <- getwd()

# functions ---------------------------------------------------------------


render_model_diagnostics <- function(prms) {
  rmarkdown::render(
    "Analysis/BiomassQuantity/Evaluate/01_model_diagnostics.Rmd",
    knit_root_dir = knit_root_dir,
    params = prms,
    output_dir = file.path(output_dir, 'Evaluate'),
    output_file = paste0("01_model_diagnostics_", make_suffix(prms), ".html")
  )
}

# setup -----------------------------------------------------------------

opts_l <- create_opts_l(suffix)
cmdargs <- create_cmdargs(opts_l)


# passed to rmd
prms_model_diagnostics <- list(
  model_version = opts_l$vm,
  data_version = if(opts_l$use_simulated) opts_l$vs else opts_l$vd,
  purer_version = opts_l$vp)

# run scripts -------------------------------------------------------------

# *data -------------------------------------------------------------------

if(run_sim & isTRUE(opts_l$use_simulated)) {

  callr::rscript("Analysis/BiomassQuantity/DataPrep/06_simulate_biomass.R", 
                 cmdargs = cmdargs)
}

# callr::rscript("Analysis/test.R",
#                cmdargs = cmdargs)

# * Fit -------------------------------------------------------------------

if(run_fit_model) {
  callr::rscript("Analysis/BiomassQuantity/Fit/01_fit_model.R", 
                 cmdargs = cmdargs)
}

# Evaluate ----------------------------------------------------------------

if(run_model_diagnostics) {
  #dir.create(file.path(output_dir, 'Evaluate'), recursive = TRUE)
  render_model_diagnostics(prms = prms_model_diagnostics)
}

if(run_model_diagnostics_sim & isTRUE(opts_l$use_simulated)) {

  rmarkdown::render(
    "Analysis/BiomassQuantity/Evaluate/01_evaluate_model_sim-data.Rmd",
    knit_root_dir = knit_root_dir,
    params = prms_model_diagnostics,
    output_dir = file.path(output_dir, 'Evaluate'),
    output_file = paste0("01_evaluate_model_sim-data_", 
                         make_suffix(prms_model_diagnostics), 
                         ".html")
  )
}
