# Purpose:
# Initial fitting model to biomass data
# this may expand, new scripts created etc.
# for starters this is for testing on simulated data


# dependencies ------------------------------------------------------------

source('Functions/init.R')
source('Functions/read_write.R')
source('Functions/models/model_specs.R')

# params ------------------------------------------------------------------


# when matures migrate to defining in 'model_specs.R'
spec <- list(
  
)

# read in data ------------------------------------------------------------

# reading output from DataProp/07_make_analysis_ready.R
dat1 <- read_analysis_ready(opt)


# find purer pixels -------------------------------------------------------
# model will be fit on subset of pixels that are purer (in terms
# of cover functional groups)

select_purer <- function() {
  
}

