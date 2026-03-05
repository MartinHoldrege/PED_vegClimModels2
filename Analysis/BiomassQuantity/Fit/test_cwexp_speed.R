# checking speed for fitting cwexp model
library(TMB)
source('Functions/init.R')
source('Functions/models/cwexp_helpers.R')
source('Functions/models/cwexp_lognormal.R')
source('Functions/models/cwexp_lognormal_tmb.R')
source('Functions/data/simulate_data.R')

# create fake data
dummy <- cwexp_make_dummy_data(n = 1e5)
dat <- dummy$data
cover_cols <- dummy$spec$cover_cols


# 1) compile + load
 dll <- cwexp_tmb_compile("src/cwexp_lognormal_tmb.cpp")
#dyn.load(TMB::dynlib('src/cwexp_lognormal_tmb'))

 # 2) fit

formula <- totalBio ~ tmean + ppt + vpd + sand

# compare speeds
results <- microbenchmark::microbenchmark(
  cwexp_fit_tmb(
    data = dat,
    formula = formula ,
    cover_cols = cover_cols,
    dll = 'cwexp_lognormal_tmb',
    control = list(iter.max = 200, eval.max = 200)
  ),
  # fitting same model but not using tmb
  cwexp_fit_lognormal(
    data = dat,
    formula = formula,
    cover_cols = cover_cols,
    control = list(maxit = 200)
  ),
  times = 1
)
results

results <- microbenchmark::microbenchmark(
  cwexp_fit_tmb(
    data = dat,
    formula = formula ,
    cover_cols = cover_cols,
    dll = 'cwexp_lognormal_tmb',
    control = list(iter.max = 200, eval.max = 200)
  ),
  # fitting same model but not using tmb
  cwexp_fit_tmb(
    data = dat,
    formula = formula ,
    cover_cols = cover_cols,
    dll = 'cwexp_lognormal_tmb',
    include_gradient = FALSE,
    control = list(iter.max = 200, eval.max = 200)
  ),
  times = 1
)
results
