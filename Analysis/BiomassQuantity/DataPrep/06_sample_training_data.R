## Purpose:
# Create training datasets for the herbaceous and woody biomass models
# by sampling CONUS rasters. Produces two data frames:
#   1. Herbaceous: RAP herb biomass ~ climate + herb cover
#   2. Woody: GEDI biomass ~ climate + shrub/tree cover
#
# Author: Martin Holdrege

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()

# params ------------------------------------------------------------------
vd <- "d07" # uses RAP cover, for code development/testing

run_herb <- TRUE
if(vd == 'd05') {
  n_sample <- 500000
  k_regions <- 20
} else if (vd == 'd06') {
  n_sample <- NULL # grab entire dataset
  k_regions <- 'L3' # use level 3 ecoregions as regions
} else if (vd == 'd07') {
  run_herb <- FALSE # this dataset is only different for woody
  n_sample <- NULL
  k_regions <- 'L2'
} 

seed <- 42
lcmap_threshold <- 0.9  # fraction of pixel that must be "keepable" land cover
fire_threshold <- 0.9 # unburned fraction, for masking woody

# output version tag


# load rasters ------------------------------------------------------------

rasters <- load_conus_rasters(
  pred_vars = NULL, # by default calculate mean/sd for all climate vars
  lcmap_threshold = lcmap_threshold,
  fire_threshold = fire_threshold,
  force_scale = FALSE,
  epa_lev = if(is.character(k_regions)) k_regions else NULL
)

r_climate_subset <- rasters$climate
r_cover <- rasters$cover
r_herb_bio <- rasters$herb_bio
r_gedi <- rasters$gedi
mask_lcmap <- rasters$mask_lcmap
mask_fire <- rasters$mask_fire
scale_df <- rasters$scale_df

pred_vars <- scale_df$variable
# === HERBACEOUS training data =============================================

# stack: response + cover + climate
r_herb_stack <- c(
  r_herb_bio,
  r_cover[["totalHerbaceousCov"]],
  r_climate_subset
)

if(k_regions == 'L3') {
  r_herb_stack <- c(r_herb_stack, rasters$region)
}

# mask to valid LCMAP pixels
r_herb_stack <- terra::mask(r_herb_stack, mask_lcmap, maskvalues = c(NA, 0))

# sample
set.seed(seed)

if(is.numeric(n_sample)) {
  
  df_herb <- terra::spatSample(
    r_herb_stack,
    size = n_sample,
    method = "random",
    na.rm = TRUE,
    xy = TRUE
  ) 
} else if(is.null(n_sample)) {
  df_herb <- terra::as.data.frame(r_herb_stack, xy = TRUE) |> 
    drop_na()
}

# rename response for model compatibility
df_herb <- dplyr::rename(df_herb, totalBio = totalHerbaceousBio) |> 
  mutate(totalBio = replace_zero(totalBio))

# standardize climate predictors using CONUS-wide scale_df
herb_std <- standardize(df_herb, vars = pred_vars, scale_df = scale_df)
df_herb <- herb_std$data

cat("  Final herbaceous dataset:", nrow(df_herb), "rows,",
    ncol(df_herb), "cols\n")

# === WOODY training data ==================================================

# stack: response + covers + climate
r_woody_stack <- c(
  r_gedi,
  r_cover[["totalTreeCov"]],
  r_cover[["totalShrubCov"]],
  r_climate_subset
)

if(is.character(k_regions)) {
  r_woody_stack <- c(r_woody_stack, rasters$region)
}
# mask to valid LCMAP pixels
r_woody_stack <- terra::mask(r_woody_stack, mask_lcmap, maskvalues = c(NA, 0))
r_woody_stack <- terra::mask(r_woody_stack, mask_fire, maskvalues = c(NA, 0))
# sample

set.seed(seed + 1)  # different seed from herbaceous


if(is.numeric(n_sample)) {
  
  df_woody <- terra::spatSample(
    r_woody_stack,
    size = n_sample,
    method = "random",
    na.rm = TRUE,
    xy = TRUE
  )
} else if(is.null(n_sample)) {
  df_woody <- terra::as.data.frame(r_woody_stack, xy = TRUE) |> 
    drop_na()
}
# rename response
df_woody <- dplyr::rename(df_woody, totalBio = totalWoodyBio) |> 
  mutate(totalBio = replace_zero(totalBio))

# standardize climate predictors using same CONUS-wide scale_df
woody_std <- standardize(df_woody, vars = pred_vars, scale_df = scale_df)
df_woody <- woody_std$data


# add regions -------------------------------------------------------------
# artificial regions for sampling purer cells



if(is.numeric(k_regions)) {
  region_vars <- c('x', 'y', 'MAT', 'MAP')
  
  df_woody$region <- suppressWarnings(make_region_kmeans(
    dat = df_woody,
    vars = region_vars,
    nstart = 5,
    k = k_regions,
    max_iter = 10000
  ))$data
  
  if(run_herb) {
    df_herb$region <- suppressWarnings(make_region_kmeans(
      dat = df_herb,
      vars = region_vars,
      nstart = 5,
      k = k_regions,
      max_iter = 10000
    ))$data
  }

}


# save outputs -------------------------------------------------------------
out_dir <- file.path(paths$large,
                     "Data_processed/BiomassQuantityData/analysis_ready")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if(run_herb) {

# herbaceous
herb_out <- list(
  data = df_herb,
  scale = scale_df,
  pred_vars = pred_vars,
  cover_cols = "totalHerbaceousCov",
  response = "totalBio",
  n_sample = n_sample,
  seed = seed,
  lcmap_threshold = lcmap_threshold,
  type = "herbaceous"
)

p_herb_out <- file.path(out_dir, paste0("biomass_herb_sample_", vd, ".rds"))
saveRDS(herb_out, p_herb_out)
cat("  Saved herbaceous:", p_herb_out, "\n")
}


# woody
woody_out <- list(
  data = df_woody,
  scale = scale_df,
  pred_vars = pred_vars,
  cover_cols = c("totalTreeCov", "totalShrubCov"),
  response = "totalBio",
  n_sample = n_sample,
  seed = seed + 1,
  lcmap_threshold = lcmap_threshold,
  fire_threshold = fire_threshold,
  type = "woody"
)

p_woody_out <- file.path(out_dir, paste0("biomass_woody_sample_", vd, ".rds"))
saveRDS(woody_out, p_woody_out)
cat("  Saved woody:", p_woody_out, "\n")

cat("\nDone.\n")

