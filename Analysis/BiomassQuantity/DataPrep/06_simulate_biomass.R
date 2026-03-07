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
source('Functions/data/simulate_data.R')
source('Functions/grouping/environmental_clustering.R')

# params ------------------------------------------------------------------

pfts <- const$pfts
# version of simulation (this script can be expanded 
# to be able to create multiple versions)
vs <- opt$vs
cover_suffix <- 'Cover_rel' # end of col names that denote cover
n_sample <- 1e5
# read in data ------------------------------------------------------------

# file output by
# 'Analysis/BiomassQuantity/DataPrep/DataToShareWithMartin.R'
p <- file.path(paths$large, 'Data_processed/BiomassQuantityData/', 
               'GEDIbiomass_modeledCover_clim_soils.csv')

dat1 <- read_csv(p)

# prep dataframe ----------------------------------------------------------

dat2 <- dat1 %>% 
  # for consistent downstream naming
  rename_with(.fn = \(x) str_replace(x, cover_suffix, 'Cov'),
              .cols = ends_with(cover_suffix)) %>% 
  filter(!is.na(biomass_MgPerHect)) %>% 
  rename(totalBio = biomass_MgPerHect)

set.seed(12)
dat2 <- sample_n(dat2, size = n_sample)

dat3 <- dat2

# make cover sum to 1 -----------------------------------------------------

cov_cols <- c(paste0(pfts, "Cov"), "bareGroundCov")

dat4 <- dat3 %>%
  mutate(.cov_total = rowSums(across(all_of(cov_cols)), na.rm = TRUE)) %>%
  mutate(across(all_of(cov_cols), ~ .x / .cov_total)) %>%
  select(-.cov_total)


# add artificial regions --------------------------------------------------
# based on climate clustering

region_vars <- c("tmean_CLIM", 
                "precip_CLIM", 
                "PrecipTempCorr_CLIM",
                "x",
                "y")

reg <- make_region_kmeans(
  dat = dat4,
  vars = region_vars,
  nstart = 1,
  k = 10,
  seed = 42
)

dat4$region <- reg$data

# define coefficients -----------------------------------------------------

pred_vars1 <- c("tmean_CLIM", 
                "precip_CLIM", 
                "PrecipTempCorr_CLIM", 
                "VPD_mean", 
                "sand")

# interactions, touples define the interactions
inter <- list(
  c('tmean_CLIM', 'precip_CLIM')
)

pred_vars2 <- c(as.list(pred_vars1), inter)
n <- length(pfts)

rcoef <- function() {
  rnorm(n = n, mean = 0,  sd = 0.4)
}

set.seed(123)
coefs1 <- map(pred_vars2, function(var) {
  coef <- rcoef()
  names(coef) <- pfts
  list(var = var,
       coef = coef)
})

intercepts <- c(3, 5, 5, 2.5, 2, 1.5)
names(intercepts) <- const$pfts

# simulated data -----------------------------------------------------------
# right now biomass of each group perfectly sums to the total model biomass
# in real data that need not 

sim <- sim_bio(data = dat4,
               coefs = coefs1,
               intercepts = intercepts,
               pred_vars = pred_vars1,
               inter = inter,
               sigma = 0.1,
               normalize = TRUE) 

# data = dat4
# coefs = coefs1
# intercepts = intercepts
# pred_vars = pred_vars1
# inter = inter
# sigma = 0.1
dat5 <- bind_cols(sim, select(dat4, -totalBio))

# write files -------------------------------------------------------------

p2 <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated', 
                paste0('simBiomass_', vs, '.csv'))

write_csv(dat5, p2)
