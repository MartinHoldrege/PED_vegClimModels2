# Purpose:
# Create exploratory figures of the simulated biomass data
# to understand what the simulations actually did


# dependencies ------------------------------------------------------------

source('Functions/init.R')
source('Functions/figFuns.R')
library(patchwork)
# params ------------------------------------------------------------------

vs <- opt$vs # version of simulated data
pfts <- const$pfts

# read in data ------------------------------------------------------------

p1 <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated', 
                paste0('simBiomass_', vs, '.csv'))

sim1 <- read_csv(p1)

# prep data ---------------------------------------------------------------

set.seed(123)
sim2 <- sample_n(sim1, size = 1e4)

cols <-  paste0(rep(pfts, 2), rep(c('Cov', 'Bio'), each = length(pfts)))
stopifnot(cols %in% names(sim2))

sim_long1 <- sim2 %>%
  pivot_longer(
    cols = all_of(cols),
    names_to = c("pft", ".value"),
    names_pattern = "(.*)(Cov|Bio)"
  )  %>%
  rename(
    cover   = Cov,
    biomass = Bio
  )

# create plots ------------------------------------------------------------

g1 <- ggplot(sim_long1, aes(cover, biomass)) +
   geom_point(alpha = 0.2) +
   geom_smooth(se = FALSE) +
   facet_wrap(~pft) +
  labs(y = 'simulated biomass')
 
g2 <- ggplot(sim_long1, aes(cover, biomass/biomass_MgPerHect)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  facet_wrap(~pft) +
  labs(y = 'Proportion of total biomass')

g <- g1 + g2

ggsave(
  filename = paste0('Figures/BiomassQuantity/simulated/',
                    'cov_vs_simBio_', vs, '.png'),
  plot = g,
  width = 14,
  height = 6
)
