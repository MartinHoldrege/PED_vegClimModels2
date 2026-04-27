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

n_sample <- 500000
seed <- 42
lcmap_threshold <- 0.9  # fraction of pixel that must be "keepable" land cover
fire_threshold <- 0.9 # unburned fraction, for masking woody

# output version tag
vd <- "d05" # uses RAP cover, for code development/testing

# load rasters ------------------------------------------------------------

rasters <- load_conus_rasters(
  pred_vars = NULL, # by default calculate mean/sd for all climate vars
  lcmap_threshold = lcmap_threshold,
  fire_threshold = fire_threshold,
  force_scale = FALSE
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

# mask to valid LCMAP pixels
r_herb_stack <- terra::mask(r_herb_stack, mask_lcmap, maskvalues = 0)

# sample
set.seed(seed)
df_herb <- terra::spatSample(
  r_herb_stack,
  size = n_sample,
  method = "random",
  na.rm = TRUE,
  xy = TRUE
)

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

# mask to valid LCMAP pixels
r_woody_stack <- terra::mask(r_woody_stack, mask_lcmap, maskvalues = 0)
r_woody_stack <- terra::mask(r_woody_stack, mask_fire, maskvalues = 0)
# sample

set.seed(seed + 1)  # different seed from herbaceous
df_woody <- terra::spatSample(
  r_woody_stack,
  size = n_sample,
  method = "random",
  na.rm = TRUE,
  xy = TRUE
)

# rename response
df_woody <- dplyr::rename(df_woody, totalBio = totalWoodyBio) |> 
  mutate(totalBio = replace_zero(totalBio))

# standardize climate predictors using same CONUS-wide scale_df
woody_std <- standardize(df_woody, vars = pred_vars, scale_df = scale_df)
df_woody <- woody_std$data


# add regions -------------------------------------------------------------
# artificial regions for sampling purer cells

region_vars <- c('x', 'y', 'MAT', 'MAP')

df_woody$region <- suppressWarnings(make_region_kmeans(
  dat = df_woody,
  vars = region_vars,
  nstart = 5,
  k = 20,
  max_iter = 10000
))$data

df_herb$region <- suppressWarnings(make_region_kmeans(
  dat = df_herb,
  vars = region_vars,
  nstart = 5,
  k = 20,
  max_iter = 10000
))$data

# save outputs -------------------------------------------------------------

out_dir <- file.path(paths$large,
                     "Data_processed/BiomassQuantityData/analysis_ready")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
