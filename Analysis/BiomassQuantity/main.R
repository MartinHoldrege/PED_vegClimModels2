# For running other scripts and rendering .rmd docs


# params ------------------------------------------------------------------
run_sim <- FALSE # simulate data
run_fit_model <- TRUE
run_model_diagnostics <- TRUE
run_model_diagnostics_sim <- FALSE

opts_l <- list(
  vm = 'm04', # version model, m01 & m03 matches the formula for simulated data
  vs = 's04', # version simulated data
  vd = 'd03', # version data # version d02 has tree cover trimmed at 1%
  # d03 trimmed at 10%
  vp = 'p02', # p01 more appropriate for the simulated data, 
  # p02 better for real (due to size)
  use_simulated = FALSE # use simulated data
)

print(opts_l)

# for rmd's
output_dir <- "Reports/BiomassQuantity/"
knit_root_dir <- getwd()

# functions ---------------------------------------------------------------

make_suffix <- function(params) {
  paste(params$data_version, params$purer_version, params$model_version,
        sep = '-')
}

create_cmdargs <- function(x) {
  x <- lapply(x, as.character)
  list(
    paste0('--vs=', x$vs),
    paste0('--vd=', x$vd),
    paste0('--vp=', x$vp),
    paste0('--vm=', x$vm),
    paste0('--use_simulated=', x$use_simulated)
  )
}

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
