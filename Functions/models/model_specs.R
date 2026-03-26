# purpose: create specifications for various model types

# build list as new model version are developed

dll_path1 <- "src/cwexp_lognormal_en_tmb.cpp"

# version that used log1p() in likelihood
dll_path2 <- "src/cwexp_lognormal_en_tmb2.cpp"

model_specs <- list(
  'm01' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM",
                  "sand"),
    inter = NULL,
    dll_path = dll_path1
  ),
    'm02' = list(
    pred_vars = c("tmean_CLIM",
                   "precip_CLIM",
                   "PrecipTempCorr_CLIM"),
    inter = NULL,
    dll_path = dll_path1
  ),
  'm03' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM",
                  "sand"),
    inter = NULL,
    dll_path = dll_path2
  ),
  'm04' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM"),
    inter = NULL,
    dll_path = dll_path2
  )
)

# for sampling purer pixels
# this is version_purer or 'vp'
purer_specs = list(
  'p01' = list(
    q = 0.9,
    min_raw_cover = 0.05
  ),
  'p02' = list(
    q = 0.98,
    min_raw_cover = 0.05
  )
)
