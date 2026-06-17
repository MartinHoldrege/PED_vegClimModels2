# For running other scripts and rendering .rmd docs

source('Functions/constants.R') # for pfts
source('Functions/models/model_specs.R')
source('Functions/general.R')

# params ------------------------------------------------------------------

run_fit_model <- FALSE
run_predict_rasters <- FALSE
run_model_diagnostics <- FALSE
run_spatial_diagnostics <- TRUE
run_model_comparison <-  TRUE
run_explore_dat_samp <- FALSE

# if TRUE, only predict with the cover source used for fitting;
# if FALSE (default), RAP-trained models also predict with model cover
predict_native_only <- FALSE


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
suffixes <- c('d08-p04.2-m13', 'd06-p04.2-m13') #, , 'd05-p04.2-m13','d05-p04-m12', 'd05-p04.2-m12', 'd05-p04.2-m13')

# comparison tuples: each element is a vector of equal length.
# Each index defines one comparison (suffix1 vs suffix2, cs1 vs cs2).
# cs means 'cover source' 
comparison_pairs <- list(
  suffix1 = c("d06-p04.2-m13"),
  suffix2 = c("d08-p04.2-m13"),
  cs1     = c("model"       ),
  cs2     = c("model"      )
)

# for exploring data sampling
suffixes_data <- c(#'d05-p04.2', 'd06-p04.2', 
  'd07-p04.2') # for plots looking input data 

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

native_cover <- data_specs[[opts_l$vd]]$cover_source

if (predict_native_only || native_cover == "model") {
  predict_cover_sources <- native_cover
} else {
  # RAP-trained models also predict with model cover
  predict_cover_sources <- c("rap", "model")
}
  
if(run_predict_rasters & hw_type) {

  
  for (model_type in model_types) {
    for (cs in predict_cover_sources) {
      callr::rscript('Analysis/BiomassQuantity/Fit/02_predict_rasters.R',
                     cmdargs = create_cmdargs(opts_l, model_type = model_type,
                                              cover_source = cs))
    }
  }
}  
  
  # Evaluate ----------------------------------------------------------------
  
if(run_model_diagnostics) {
  #dir.create(file.path(output_dir, 'Evaluate'), recursive = TRUE)
  if(hw_type) {
    for (model_type in model_types) {
      for (cs in predict_cover_sources) {
        prms <- prms_model_diagnostics
        prms$model_type <- model_type
        prms$cover_source <- cs
        
        out_suffix <- make_suffix(prms)
        if (cs != 'rap') out_suffix <- paste0(out_suffix, "_cov-", cs)
        
        rmarkdown::render(
          "Analysis/BiomassQuantity/Evaluate/01_model_diagnostics_hw.Rmd",
          knit_root_dir = knit_root_dir,
          params = prms,
          output_dir = file.path(output_dir, 'Evaluate'),
          output_file = paste0("01_model_diagnostics_", out_suffix, ".html")
        )
      }
    }
  } else {
    render_model_diagnostics(
       path = "Analysis/BiomassQuantity/Evaluate/01_model_diagnostics.Rmd",
       prms = prms_model_diagnostics)
  }
}
  
if(run_spatial_diagnostics & hw_type) {
  for (model_type in model_types) {
    for (cs in predict_cover_sources) {
      prms <- prms_model_diagnostics
      prms$model_type <- model_type
      prms$cover_source <- cs
      
      out_suffix <- make_suffix(prms)
      if (cs != "rap") out_suffix <- paste0(out_suffix, "_cov-", cs)
      
      rmarkdown::render(
        "Analysis/BiomassQuantity/Evaluate/01_spatial_diagnostics.Rmd", 
        knit_root_dir = knit_root_dir,
        params = prms,
        output_dir = file.path(output_dir, 'Evaluate'),
        output_file = paste0("01_spatial_diagnostics_", out_suffix, ".html")
      )
    }
  }
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
  n_comparisons <- length(comparison_pairs$suffix1)
  
  for (i in seq_len(n_comparisons)) {
    opts1 <- create_opts_l(comparison_pairs$suffix1[i])
    opts2 <- create_opts_l(comparison_pairs$suffix2[i])
    cs1_i <- comparison_pairs$cs1[i]
    cs2_i <- comparison_pairs$cs2[i]
    
    # compare within each model type
    model_types_1 <- names(model_specs[[opts1$vm]])
    model_types_2 <- names(model_specs[[opts2$vm]])
    shared_types <- intersect(model_types_1, model_types_2)
    
    for (mt in shared_types) {
      if ((mt == "herb" & !fit_herb) | (mt == "woody" & !fit_woody)) next
      
      prms <- list(
        model_type = mt,
        vd1 = opts1$vd, vp1 = opts1$vp, vm1 = opts1$vm, cs1 = cs1_i,
        vd2 = opts2$vd, vp2 = opts2$vp, vm2 = opts2$vm, cs2 = cs2_i,
        n_sample_raster = 50000
      )
      
      # build output filename reflecting cover sources when they differ
      tag1 <- comparison_pairs$suffix1[i]
      tag2 <- comparison_pairs$suffix2[i]
      if (cs1_i != "rap") tag1 <- paste0(tag1, "_cov-", cs1_i)
      if (cs2_i != "rap") tag2 <- paste0(tag2, "_cov-", cs2_i)
      
      out_file <- paste0("01_model_comparison_", mt, "_",
                         tag1, "_vs_", tag2, ".html")
      
      rmarkdown::render(
        "Analysis/BiomassQuantity/Evaluate/01_model_comparison.Rmd",
        knit_root_dir = knit_root_dir,
        params = prms,
        output_dir = file.path(output_dir, "Evaluate", "comparison"),
        output_file = out_file
      )
    }
  }
}
