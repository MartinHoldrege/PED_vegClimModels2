# For running other scripts and rendering .rmd docs

source('Functions/constants.R') # for pfts
source('Functions/models/model_specs.R')
source('Functions/general.R')

# params ------------------------------------------------------------------

run_fit_model <- FALSE
run_predict_rasters <- FALSE
run_model_diagnostics <- FALSE
run_spatial_diagnostics <- FALSE
run_model_comparison <- FALSE
run_explore_dat_samp <- TRUE


# with PFTs to fit model to?
fit_woody = TRUE
fit_herb = FALSE

# simulation specific
run_sim <- FALSE # simulate data
run_model_diagnostics_sim <- FALSE

# [vs]: data version (s for simulated, d- for real data), 
#       
# p for purer version (sampling purer pixels)
#     p01 more appropriate for the simulated data, (p02 for real data)
# m for model version
#     m01 & m03 matches the formula for simulated data

# repo updated to only work w/ separate herb/biomass models
# meaning: vd >= d05, vp >= p04, and vm >= 0
suffixes <- c('d06-p04.2-m13', 'd05-p04.2-m13') #,'d05-p04-m12', 'd05-p04.2-m12', 'd05-p04.2-m13')

# pairs of suffixes to compare (model1 vs model2)
comparison_pairs <- list(
  # c("d05-p04-m13", "d05-p04.2-m13"),
  c("d05-p04.2-m13", "d06-p04.2-m13")
)

# for exploring data sampling
suffixes_data <- c('d05-p04.2', 'd06-p04.2') # for plots looking input data 

# for rmd's
output_dir <- "Reports/BiomassQuantity/"
knit_root_dir <- getwd()

# functions ---------------------------------------------------------------


render_model_diagnostics <- function(prms, 
                                     path = "Analysis/BiomassQuantity/Evaluate/01_model_diagnostics_hw.Rmd") {
  rmarkdown::render(
    path, 
    knit_root_dir = knit_root_dir,
    params = prms,
    output_dir = file.path(output_dir, 'Evaluate'),
    output_file = paste0("01_model_diagnostics_", make_suffix(prms), ".html")
  )
}

for (suffix in suffixes) {
  # setup -----------------------------------------------------------------
  
  opts_l <- create_opts_l(suffix)
  vm <- opts_l$vm
  
  # validate that vm and vp keys exist in their respective spec lists
  if(!opts_l$vm %in% names(model_specs)) {
    paste0("vm '", opts_l$vm, "' not found in model_specs")
    
  }
  # e.g. if fitting woody then the woody sublist needs to exist
  stopifnot(!fit_woody | !is.null(model_specs[[vm]]$woody),
            !fit_herb | !is.null(model_specs[[vm]]$herb))
  
  if(!opts_l$vp %in% names(purer_specs)) {
    paste0("vp '", opts_l$vp, "' not found in purer_specs") 
  }
  stopifnot(!fit_woody | !is.null(purer_specs[[opts_l$vp]]$woody),
            !fit_herb | !is.null(purer_specs[[opts_l$vp]]$herb))


  print(opts_l)
  hw_type <- vm >= 'm09' # model types that separately fit models to herb and woody
  
  if(hw_type) {
    model_types0 <- names(model_specs[[opts_l$vm]])
    
    model_types <- vector('character')
    if(fit_herb & ('herb' %in% model_types0))  {
      model_types <- c(model_types, 'herb')
      fit_herb2 <- fit_herb
    } else {
      # some aren't run for herbs b/ they're not specified in model_specs
      fit_herb2 <- FALSE
    }
    if(fit_woody& ('woody' %in% model_types0)) {
      model_types <- c(model_types, 'woody')
      fit_woody2 <- fit_woody
    } else {
      fit_woody2 <- FALSE
    }
    stopifnot(length(model_types) >= 1)
    
    # currently setup for herb & woody  models
  } else {
    model_types <- NULL
  }
  
  cmdargs <- create_cmdargs(opts_l, fit_herb = fit_herb2, 
                            fit_woody = fit_woody2)

  # passed to rmd
  prms_model_diagnostics <- list(
    model_version = opts_l$vm,
    data_version = if(opts_l$use_simulated) opts_l$vs else opts_l$vd,
    purer_version = opts_l$vp
    )
  
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
    p <- if(hw_type) {
      "Analysis/BiomassQuantity/Fit/01_fit_model_hw.R"
    } else {
      "Analysis/BiomassQuantity/Fit/01_fit_model.R"
    }
    callr::rscript(p, cmdargs = cmdargs)
  }
  

# predict -----------------------------------------------------------------

if(run_predict_rasters & hw_type) {
  lapply(model_types, \(model_type) {
      callr::rscript('Analysis/BiomassQuantity/Fit/02_predict_rasters.R',
                     cmdargs = create_cmdargs(opts_l, model_type = model_type))
  })

}  
  
  # Evaluate ----------------------------------------------------------------
  
if(run_model_diagnostics) {
  #dir.create(file.path(output_dir, 'Evaluate'), recursive = TRUE)
  if(hw_type) {
    lapply(model_types, function(model_type) {
      prms_model_diagnostics$model_type <- model_type
      render_model_diagnostics(prms = prms_model_diagnostics)
    })
  } else {
    render_model_diagnostics(
       path = "Analysis/BiomassQuantity/Evaluate/01_model_diagnostics.Rmd",
       prms = prms_model_diagnostics)
  }
}
  
if(run_spatial_diagnostics & hw_type) {
    lapply(model_types, function(model_type) {
      prms <- prms_model_diagnostics
      prms$model_type <- model_type
      rmarkdown::render(
        "Analysis/BiomassQuantity/Evaluate/01_spatial_diagnostics.Rmd", 
        knit_root_dir = knit_root_dir,
        params = prms,
        output_dir = file.path(output_dir, 'Evaluate'),
        output_file = paste0("01_spatial_diagnostics_", make_suffix(prms), ".html")
      )
  })
}
  
  
if(run_model_diagnostics_sim & isTRUE(opts_l$use_simulated) & !hw_type) {

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

}

# data prep exploration -------------------------------------------------------

if(run_explore_dat_samp) {

  for(suffix_data in suffixes_data) {
    
    for(model_type in c('woody', 'herb')[c(fit_woody, fit_herb)]) {
        
  
      # adding a dummy model version to get function to work
      opts_i <- create_opts_l(paste0(suffix_data, '-m00')) 
      
      prms_i <- list(
        data_version = if(opts_i$use_simulated) opts_i$vs else opts_i$vd,
        purer_version = opts_i$vp,
        model_type = model_type)
      
      
      rmarkdown::render(
        "Analysis/BiomassQuantity/DataPrep/08_explore_data_sampling.Rmd",
        knit_root_dir = knit_root_dir,
        params = prms_i,
        output_dir = file.path(output_dir, 'DataPrep'),
        output_file = paste0("08_explore_data_sampling_", model_type, '_', 
                             suffix_data, 
                             ".html")
      )
    }
  }
}
# model comparison --------------------------------------------------------
  
if (run_model_comparison) {
  for (pair in comparison_pairs) {
    opts1 <- create_opts_l(pair[1])
    opts2 <- create_opts_l(pair[2])
    
    # compare within each model type
    model_types_1 <- names(model_specs[[opts1$vm]])
    model_types_2 <- names(model_specs[[opts2$vm]])
    shared_types <- intersect(model_types_1, model_types_2)
    
    for (mt in shared_types) {
      if ((mt == "herb" & !fit_herb) | (mt == "woody" & !fit_woody)) next
    }
    model_types_2 <- names(model_specs[[opts2$vm]])
    shared_types <- intersect(model_types_1, model_types_2)
    
    for (mt in shared_types) {
      if ((mt == "herb" & !fit_herb) | (mt == "woody" & !fit_woody)) next
      
      prms <- list(
        model_type = mt,
        vd1 = opts1$vd, vp1 = opts1$vp, vm1 = opts1$vm,
        vd2 = opts2$vd, vp2 = opts2$vp, vm2 = opts2$vm,
        n_sample_raster = 50000
      )
      
      out_file <- paste0("01_model_comparison_", mt, "_",
                         pair[1], "_vs_", pair[2], ".html")
      
      rmarkdown::render(
        "Analysis/BiomassQuantity/Evaluate/01_model_comparison.Rmd",
        knit_root_dir = knit_root_dir,
        params = prms,
        output_dir = file.path(output_dir, "Evaluate", 'comparison'),
        output_file = out_file
      )
    }
  }
}
