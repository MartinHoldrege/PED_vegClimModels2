# Purpose: 
# simulate biomass of individual functional types
# to create a fake dataset to test models against
# The predicted cover, predictor variables are 
# real. we're simulating total biomass and
# we have 'true' (simulated) functional group biomass that we can test against
# (normally this isn't 'seen')

# Author Martin Holdrege
# Script Started January 26, 2026

# dependencies ------------------------------------------------------------

source('Functions/init.R')
source_functions()

# params ------------------------------------------------------------------

pfts <- const$pfts
# version of simulation (this script can be expanded 
# to be able to create multiple versions)
vs <- opt$vs
cover_suffix <- 'Cover_rel' # end of col names that denote cover
n_sample <- 1e5

# noise parameters
sigma_pft <- 0.2    # per-PFT linear predictor noise, multiplicative
sigma_obs <- 0.4     # log-scale observation noise on totalBio, noise is proportional
sigma_region <- 0.4  # region-level correlated bias on linear predictor, multiplicative

# s01 was part of earlier code and can't be recreated without
# going back to an older commit
trim_tree_cov <- NULL
if(vs == 's02') {
  intercepts <- c(20, 148, 148, 12, 7, 4) # s02 intercepts
} else if (vs == 's03') {
  intercepts <- c(2, 148, 148, 1.2, 0.7, 0.4) # s03 intercepts
} else if (vs == 's04') {
  intercepts <- c(2, 148, 148, 1.2, 0.7, 0.4) 
  trim_tree_cov <- 0.01 # small tree covers become zero
}

# read in data ------------------------------------------------------------

# file output by
# 'Analysis/BiomassQuantity/DataPrep/DataToShareWithMartin.R'
p <- file.path(paths$large, 'Data_processed/BiomassQuantityData/', 
               'GEDIbiomass_modeledCover_clim_soils.csv')

dat1 <- read_csv(p)

# prep dataframe ----------------------------------------------------------

dat4_l <- prepare_d01(data = dat1, cover_suffix = cover_suffix, pfts = pfts,
                      trim_tree_cov = trim_tree_cov)

dat4 <- dat4_l$data
nrow(dat4)
set.seed(12)
dat4 <- sample_n(dat4, size = n_sample)

# define coefficients -----------------------------------------------------

pred_vars1 <- c("tmean_CLIM", 
                "precip_CLIM", 
                "PrecipTempCorr_CLIM", 
                "sand")

# interactions, touples define the interactions
inter <- list(
  # c('tmean_CLIM', 'precip_CLIM'),
  NULL
)
inter <- purrr::keep(inter, \(x) !is.null(x))
pred_vars2 <- c(as.list(pred_vars1), inter)
n <- length(pfts)


sd <- intercepts*0.2 # making climate effects proportional to intercept

rcoef <- function() {
  rnorm(n = n, mean = 0,  sd = sd)
}

set.seed(123)
coefs1 <- map(pred_vars2, function(var) {
  coef <- rcoef()
  names(coef) <- pfts
  list(var = var,
       coef = coef)
})

names(intercepts) <- const$pfts

# simulated data -----------------------------------------------------------
# right now biomass of each group perfectly sums to the total model biomass
# in real data that need not 

sim <- sim_bio(data = dat4,
               coefs = coefs1,
               intercepts = intercepts,
               pred_vars = pred_vars1,
               inter = inter,
               sigma_pft = sigma_pft,
               sigma_obs = sigma_obs,
               sigma_region = sigma_region,
               region = dat4$region,
               # already normalized
               normalize = FALSE) 

with(sim$data, cor(totalMu, totalBio))

sim$scale <- dat4_l$scale

# write files -------------------------------------------------------------


# save full sim object (includes true parameters for validation)
p3 <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated',
                paste0('simBiomass_', vs, '.rds'))

saveRDS(sim, p3)
