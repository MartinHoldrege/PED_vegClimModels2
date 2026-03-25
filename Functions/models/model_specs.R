# purpose: create specifications for various model types

# build list as new model version are developed
model_specs <- list(
  'm01' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM",
                  "sand"),
    
    inter = NULL
  ),
    'm02' = list(
    pred_vars = c("tmean_CLIM",
                   "precip_CLIM",
                   "PrecipTempCorr_CLIM"),
    
    inter = NULL
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
