# misc functions for reading/writing data

#' read the clean data file for start of modeling pipeline
#' (generally output of DataPrep/07_make_analysis_ready.R)
#'
#' @param opt list of options (init.R needs to have been run)
#' @param vd vrsion of data, starts with s (simulated) or d (not simulated)
#' @param root where data is, by default needs init.R to be run
#'
#' @returns
#' dataframe
read_analysis_ready <- function(opt = NULL, vd = NULL, root = paths$large) {
  stopifnot(is.list(opt) | is.character(vd),
            is.null(opt) | c('use_simulated', 'vs', 'vd') %in% names(opt),
            dir.exists(root),
            !(is.null(opt) & is.null(vd))
            )
  
  if(is.null(opt)) {
    use_simulated <- str_detect(vd, '^s')
    v <- vd
  } else {
    use_simulated <- opt$use_simulated
    v <- if(opt$use_simulated) opt$vs else opt$vd
  }
  
  if(use_simulated) {
    
    p <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated',
                    paste0('simBiomass_', v, '.rds'))
    data <- readRDS(p)$data
  } else if(v == 'd01') {
    data <- read_prepare_d01(root = root)$data
  } else if(v == 'd02') {
    data <- read_prepare_d01(root = root, trim_tree_cov = 0.01)$data
  } else if(v == 'd03') {
    data <- read_prepare_d01(root = root, trim_tree_cov = 0.1)$data
  } else{
    p <- file.path(paths$large, 
                   'Data_processed/BiomassQuantityData/analysis_ready', 
                   paste0('biomass_', v, '.rds'))
    data <- readRDS(p)
  }
  
  data
}

read_prepare_d01 <- function(root = paths$large,
                             trim_tree_cov = NULL) {
  p <- file.path(root, 'Data_processed/BiomassQuantityData/', 
                 'GEDIbiomass_modeledCover_clim_soils.csv')
  
  dat1 <- read_csv(p)
  prepare_d01(dat1, cover_suffix = "Cover_rel",pfts = const$pfts,
              trim_tree_cov = trim_tree_cov)
} 
