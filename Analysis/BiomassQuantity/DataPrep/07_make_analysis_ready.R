# Purpose: 
# Create a single dataframe that is used in downstream 
# model fitting. Basically the goal here is to create
# the same formatted dataframe no matter the 
# source (simulated data, actual data, biomass source etc.).
# Develop as needed to keep inputs from different sources as consistent as
# possible

# Author: Martin Holdrege

# Started February 27, 2026

# dependencies ------------------------------------------------------------

source('Functions/init.R')

# params ------------------------------------------------------------------

# read in data ------------------------------------------------------------

if(opt$use_simulated) {
  p <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated', 
                  paste0('simBiomass_', opt$vs, '.csv'))
  dat1 <- read_csv(p)
  v <- opt$vs
} else {
  stop('code not updated')
}


# data prep ---------------------------------------------------------------
# develop code here as needed
dat3 <- dat1
# write file --------------------------------------------------------------

p_out <- file.path(paths$large, 'Data_processed/BiomassQuantityData/analysis_ready', 
                   paste0('biomass_', v, '.rds'))

saveRDS(dat3, p_out)
