# misc functions for reading/writing data

#' read the clean data file for start of modeling pipeline
#' (generally output of DataPrep/07_make_analysis_ready.R)
#'
#' @param opt list of options (init.R needs to have been run)
#' @param root where data is, by default needs init.R to be run
#'
#' @returns
#' dataframe
read_analysis_ready <- function(opt, root = paths$large) {
  stopifnot(is.list(opt),
            c('use_simulated', 'vs', 'vd') %in% names(opt),
            dir.exists(root)
            )
  
  
  v <- if(opt$use_simulated) opt$vs else opt$vd
  
  p <- file.path(paths$large, 
                'Data_processed/BiomassQuantityData/analysis_ready', 
                paste0('biomass_', v, '.rds'))
  readRDS(p)
}
