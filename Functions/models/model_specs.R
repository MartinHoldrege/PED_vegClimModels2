# purpose: create specifications for various model types

# build list as new model version are developed

dll_path1 <- "src/cwexp_lognormal_en_tmb.cpp"

# version that used log1p() in likelihood
dll_path2 <- "src/cwexp_lognormal_en_tmb2.cpp"

pfts1 <- c("shrub", "needleLeavedTree", "broadLeavedTree", "C3Gram", "C4Gram", 
           "Forb")

pfts2 <- c("shrub", 'totalTree', 'totalHerbaceous')

model_specs <- list(
  'm01' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM",
                  "sand"),
    inter = NULL,
    dll_path = dll_path1,
    pfts = pfts1
  ),
    'm02' = list(
    pred_vars = c("tmean_CLIM",
                   "precip_CLIM",
                   "PrecipTempCorr_CLIM"),
    inter = NULL,
    dll_path = dll_path1,
    pfts = pfts1
  ),
  'm03' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM",
                  "sand"),
    inter = NULL,
    dll_path = dll_path2,
    pfts = pfts1
  ),
  'm04' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM"),
    inter = NULL,
    dll_path = dll_path2,
    pfts = pfts1
  ),
  'm05' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM"),
    inter = NULL,
    dll_path = dll_path2,
    pfts = pfts2
  ),
  'm06' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM"),
    inter = NULL,
    # testing version that limits alpha to >=0
    dll_path = "src/cwexp_lognormal_en_tmb3.cpp",
    pfts = pfts2
  ),
  'm07' = list(
    pred_vars = c("tmean_CLIM",
                  "precip_CLIM",
                  "PrecipTempCorr_CLIM"),
    inter = NULL,
    # testing version that limits alpha to >=0
    dll_path = dll_path2,
    pfts = pfts2,
    fix_alpha_pfts = c('totalHerbaceous'),
    fix_alpha_filter = list(
      # columns to require near-zero for the subset
      exclude_cols = c('totalTreeCov', 'shrubCov'),
      max_cover = c(0.1, 0.1)  # threshold for "absent"
    )),
    'm08' = list(
      pred_vars = c("tmean_CLIM",
                    "precip_CLIM",
                    "PrecipTempCorr_CLIM",
                    "sand"),
      inter = NULL,
      # testing version that limits alpha to >=0
      dll_path = dll_path2,
      pfts = pfts1,
      fix_alpha_pfts = c("C3Gram", "C4Gram", "Forb"),
      fix_alpha_filter = list(
        # columns to require near-zero for the subset
        exclude_cols = c('totalTreeCov', 'shrubCov'),
        # usually fit to data where cover is pre-trimmed
        max_cover = c(0.01, 0.01)  # threshold for "absent"
    ))
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
  ),
  'p03' = list(
    q = 0,
    min_raw_cover = 0.05
  )
)
