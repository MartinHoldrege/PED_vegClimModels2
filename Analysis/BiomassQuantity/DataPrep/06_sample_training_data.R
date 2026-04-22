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
# which climate variables to include (short names)
pred_vars_herb <- c("MAT", "MAP", "PrecipTempCorr", "WD_mean", "VPD_mean",
                    "P_seasonality", "T_coldestMonth", "T_warmestMonth", 
                    "frost_free_days", "P_wettestMonth")

pred_vars_woody <- pred_vars_herb

# output version tag
vd <- "d05" # uses RAP cover, for code development/testing

# file paths --------------------------------------------------------------

years <- '2000-2023'

p_climate <- file.path(
  paths$large,
  "Data_processed/BiomassQuantityData",
  "DayMetData_allCONUS_2023ClimateValues_raster.tif"
)

p_lcmap <- file.path(
  paths$large,
  "Data_processed/masks",
  "LCMAP_fracKeep_1000m.tif"
)

p_fire <- file.path(
  paths$large,
  "Data_processed/masks",
  paste0("MTBS_fracUnburned_", years, "_1000m.tif")
)

p_cover <- file.path(
  paths$large,
  "Data_processed/CoverData/rap",
  paste0("RAP_v3_cover_", years, "_1000m.tif")
)

p_herb_bio <- file.path(
  paths$large,
  "Data_processed/BiomassQuantityData/rap",
  paste0("RAP_v3_herbaceousAGB_mask-Lcmap", lcmap_threshold*100, "_", years, "_1000m.tif")
)

p_gedi <- file.path(
  paths$large,
  "Data_processed/BiomassQuantityData",
  "GEDI_biomassRaster_onDayMetGrid.tif"
)

# load rasters (lazy — no data in memory) ---------------------------------

r_climate <- read_climate_raster(p_climate)

r_lcmap <- terra::rast(p_lcmap)
r_fire <- terra::rast(p_fire)

r_cover <- terra::rast(p_cover)

r_herb_bio <- terra::rast(p_herb_bio)
names(r_herb_bio) <- "totalHerbaceousBio"

r_gedi <- terra::rast(p_gedi)
names(r_gedi) <- "totalWoodyBio"


# align rasters -----------------------------------------------------------

aligned <- align_raster_extents(
  rast_list = list(
    climate = r_climate,
    lcmap = r_lcmap,
    fire = r_fire,
    cover = r_cover,
    herb_bio = r_herb_bio,
    gedi = r_gedi
  )
)

r_climate  <- aligned$climate
r_lcmap    <- aligned$lcmap
r_fire     <- aligned$fire
r_cover    <- aligned$cover
r_herb_bio <- aligned$herb_bio
r_gedi     <- aligned$gedi

# create masks ------------------------------------------------------------

# binary LCMAP mask (1 where fraction keepable >= threshold)
lcmap_lyr <- paste0('fracKeep_gte', lcmap_threshold*100) 
if(lcmap_lyr %in% names(r_lcmap)) {
  mask_lcmap <- r_lcmap[[lcmap_lyr]]
} else {
  message('lcmap threshold layer missing,... computing')
  mask_lcmap <- r_lcmap[["fracKeep"]] >= lcmap_threshold
}

fire_lyr <- paste0('fracUnburned_gte', fire_threshold*100) 
if(fire_lyr %in% names(r_fire)) {
  mask_fire <- r_fire[[fire_lyr]]
} else {
  message('fire threshold layer missing,... computing')
  mask_fire <- r_fire[["fracUnburned"]] >= fire_threshold
}

# select climate layers and compute scale_df -------------------------------

all_pred_vars <- unique(c(pred_vars_herb, pred_vars_woody))
r_climate_subset <- r_climate[[all_pred_vars]]

# compute standardization parameters from the full CONUS climate raster

scale_df <- compute_scale_df(rast = r_climate_subset, 
                 vars = all_pred_vars, 
                 source_path = p_climate,
                 force = FALSE)


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
herb_std <- standardize(df_herb, vars = pred_vars_herb, scale_df = scale_df)
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
woody_std <- standardize(df_woody, vars = pred_vars_woody, scale_df = scale_df)
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
  pred_vars = pred_vars_herb,
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
  pred_vars = pred_vars_woody,
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
