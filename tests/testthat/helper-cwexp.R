# testhat automatically sources all 'helper-....R' files first
# to run all tests
# the following in the console: testthat::test_dir("tests/testthat")
library(testthat)
library(dplyr)

# Project root
root <- here::here() # needed b/ testthat often treats wd as tests/testthat

# Source project code
source(file.path(root, "Functions", "data", "simulate_data.R"))
source(file.path(root, "Functions", "models", "cwexp_helpers.R"))
source(file.path(root, "Functions", "models", "cwexp_lognormal_tmb.R"))


set.seed(1)
cwexp_dummy <- cwexp_make_dummy_data(n = 1000)

cwexp_dll <- cwexp_tmb_compile(file.path(root, "src/cwexp_lognormal_tmb.cpp"), 
                               quiet = TRUE)
cwexp_dll_en <- cwexp_tmb_compile(file.path(root, "src/cwexp_lognormal_en_tmb.cpp"), 
                                  quiet = TRUE)

cwexp_test_formula <- totalBio ~ tmean + ppt + vpd + sand
