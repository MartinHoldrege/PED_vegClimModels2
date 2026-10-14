# misc functions for reading/writing data

#' read the clean data file for start of modeling pipeline
#' (generally output of DataPrep/07_make_analysis_ready.R)
#'
#' @param opt list of options (init.R needs to have been run)
#' @param vd vrsion of data, starts with s (simulated) or d (not simulated)
#' @param root where data is, by default needs init.R to be run
#'
#' @returns
#' dataframe
read_analysis_ready <- function(opt = NULL, vd = NULL, root = paths$large,
                                only_dataframe = TRUE) {
  stopifnot(is.list(opt) | is.character(vd),
            is.null(opt) | c('use_simulated', 'vs', 'vd') %in% names(opt),
            dir.exists(root),
            !(is.null(opt) & is.null(vd))
            )
  
  if(is.null(opt)) {
    use_simulated <- str_detect(vd, '^s')
    v <- vd
  } else {
    use_simulated <- opt$use_simulated
    v <- if(opt$use_simulated) opt$vs else opt$vd
  }
  
  if(use_simulated) {
    
    p <- file.path(paths$large, 'Data_processed/BiomassQuantityData/simulated',
                    paste0('simBiomass_', v, '.rds'))
    obj <- readRDS(p) 
  } else if(v == 'd01') {
    obj <- read_prepare_d01(root = root) 
  } else if(v == 'd02') {
    obj <- read_prepare_d01(root = root, trim_tree_cov = 0.01) 
  } else if(v == 'd03') {
    obj <- read_prepare_d01(root = root, trim_tree_cov = 0.1) 
  } else if(v == 'd04') {
    obj <- read_prepare_d01(root = root, trim_tree_cov = 0.1,
                             trim_shrub_cov = 0.1) 
  } else{
    stop('function not updated to read in provided data version')
  }
  
  if(only_dataframe) {
    return(obj$data)
  }
  obj
}

read_analysis_ready_hw <- function(vd, model_type, only_dataframe = FALSE) {
  
  p_data <- file.path(
    paths$large,
    "Data_processed/BiomassQuantityData/analysis_ready",
    paste0("biomass_", model_type, "_sample_", vd, ".rds")
  )
  stopifnot(file.exists(p_data))
  data_obj <- readRDS(p_data)
  
  if(only_dataframe) {
    return(data_obj$data)
  }
  data_obj
}

read_prepare_d01 <- function(root = paths$large,
                             trim_tree_cov = NULL,
                             trim_shrub_cov = NULL,
                             pfts = const$pfts1) {
  p <- file.path(root, 'Data_processed/BiomassQuantityData/', 
                 'GEDIbiomass_modeledCover_clim_soils.csv')
  
  dat1 <- read_csv(p)
  prepare_d01(dat1, cover_suffix = "Cover_rel",pfts = pfts,
              trim_tree_cov = trim_tree_cov,
              trim_shrub_cov = trim_shrub_cov)
} 

#' Read a climate raster (current or projected) and rename layers
#'
#' @param scenario One of "current" (Daymet), "BNU-ESM", or "IPSL-CM5A-MR".
#' @param path Optional explicit path; overrides `scenario` for the climate file.
#' @param path_soil Path to soil AWC raster, or NULL to skip.
#' @return SpatRaster with short-named climate layers (and soil if provided).
read_climate_raster <- function(
    scenario = c("current", "BNU-ESM", "IPSL-CM5A-MR"),
    path = NULL,
    path_soil = file.path(paths$large, 
              "./Data_processed/soils/", 
              "soil_covariates_solus100_1000m.tif")
    ) {
  scenario <- match.arg(scenario)
  
  if (is.null(path)) {
    path <- if (scenario == "current") {
      file.path(paths$large, "Data_processed/WallToWallClimateData", 
                "DaymetClimateData_1991-2020_CLIM.tif")
    } else {
      file.path(paths$large, "Data_processed/WallToWallClimateData",
                paste0("ForecastedClimateData_", scenario, "_rcp85_CLIM.tif"))
    }
  }
  
  r <- terra::rast(path)
  names(r) <- climate_name_lookup(names(r))
  
  if(!is.null(path_soil)) {
    r_soil <- terra::rast(path_soil)
    names(r_soil) <- stringr::str_replace(names(r_soil), '_cm', '')
    names(r_soil)[names(r_soil) == "AWHC"] <- climate_name_lookup("AWHC")
    r <- c(r, r_soil)
  }
  r
}

#' Load and align CONUS rasters for biomass modeling
#'
#' Loads climate, cover, biomass response, and mask rasters, aligns
#' extents and CRS, and applies masks. All operations are lazy (no data
#' loaded into memory). The scale_df computation is cached via
#' `compute_scale_df()`.
#'
#' @param pred_vars Character vector of climate variable short names.
#' @param lcmap_threshold Numeric; LCMAP fraction keep threshold.
#' @param fire_threshold Numeric; MTBS fraction unburned threshold.
#' @param root Root path for large files.
#' @param years Character; year range for time-varying products.
#' @param force_scale Logical; if TRUE, recompute scale_df even if cached.
#'
#' @return A named list with elements: `climate` (subset, unstandardized),
#'   `cover`, `herb_bio`, `gedi`, `mask_lcmap`, `mask_fire`, `scale_df`,
#'   `paths`, `params`.
#' @export
load_conus_rasters <- function(pred_vars = NULL,
                               lcmap_threshold = 0.9,
                               fire_threshold = 0.9,
                               root = paths$large,
                               force_scale = FALSE,
                               epa_lev = 'L3',
                               cover_source = 'rap'
                               ) {

  years = "2000-2023"
  # file paths
  p_climate <- file.path(root, "Data_processed/BiomassQuantityData",
                         "DayMetData_allCONUS_2023ClimateValues_raster.tif")
  p_lcmap <- file.path(root, "Data_processed/masks",
                       "LCMAP_fracKeep_1000m.tif")
  p_fire <- file.path(root, "Data_processed/masks",
                      paste0("MTBS_fracUnburned_", years, "_1000m.tif"))
  
  p_herb_bio <- file.path(root, "Data_processed/BiomassQuantityData/rap",
                          paste0("RAP_v3_herbaceousAGB_mask-Lcmap",
                                 lcmap_threshold * 100, "_", years, "_1000m.tif"))
  p_gedi <- file.path(root, "Data_processed/BiomassQuantityData",
                      "GEDI_biomassRaster_onDayMetGrid.tif")
  p_fracNotForest <- file.path(root, "Data_processed/CoverData/rap",
    paste0("RAP_v3_fracNotForest_lt3_mask-lcmap50_2019-2023_1000m.tif"))
  
  
  all_paths <- c(p_climate, p_lcmap, p_fire, 
                 p_herb_bio, p_gedi, p_fracNotForest)
  
  missing <- all_paths[!file.exists(all_paths)]
  if (length(missing) > 0) {
    stop("Missing raster files:\n  ", paste(missing, collapse = "\n  "))
  }
  
  cat("Loading CONUS rasters...\n")
  
  r_climate <- read_climate_raster(p_climate)
  r_lcmap <- terra::rast(p_lcmap)
  r_fire <- terra::rast(p_fire)
  
  r_region <- load_ecoregion_raster(epa_lev = epa_lev)
  
  r_herb_bio <- terra::rast(p_herb_bio)
  names(r_herb_bio) <- "totalHerbaceousBio"
  
  r_gedi <- terra::rast(p_gedi)
  names(r_gedi) <- "totalWoodyBio"
  
  r_cover <- load_cover(cover_source = cover_source, root = root)
  
  r_zero_tree = terra::rast(p_fracNotForest) > 0.5 # places we're considering 0 trees
  # align
  aligned <- align_raster_extents(rast_list = list(
    climate = r_climate, lcmap = r_lcmap, fire = r_fire,
    cover = r_cover, herb_bio = r_herb_bio, gedi = r_gedi,
    zero_tree = r_zero_tree, region = r_region)
  )
  
  # masks
  lcmap_lyr <- paste0("fracKeep_gte", lcmap_threshold * 100)
  mask_lcmap <- if (lcmap_lyr %in% names(aligned$lcmap)) {
    aligned$lcmap[[lcmap_lyr]]
  } else {
    message("LCMAP threshold layer missing, computing...")
    aligned$lcmap[["fracKeep"]] >= lcmap_threshold
  }
  
  fire_lyr <- paste0("fracUnburned_gte", fire_threshold * 100)
  mask_fire <- if (fire_lyr %in% names(aligned$fire)) {
    aligned$fire[[fire_lyr]]
  } else {
    message("Fire threshold layer missing, computing...")
    aligned$fire[["fracUnburned"]] >= fire_threshold
  }
  
  if(is.null(pred_vars)) {
    pred_vars <- names(aligned$climate)
  }
  # climate subset + cached scale_df
  r_climate_subset <- aligned$climate[[pred_vars]]
  
  scale_df <- compute_scale_df(
    rast = r_climate_subset,
    vars = pred_vars,
    source_path = p_climate,
    force = force_scale
  )
  
  cat("  Done loading CONUS rasters.\n")
  
  list(
    climate = r_climate_subset,
    cover = aligned$cover,
    herb_bio = aligned$herb_bio,
    gedi = aligned$gedi,
    zero_tree = aligned$zero_tree,
    region = aligned$region,
    mask_lcmap = mask_lcmap,
    mask_fire = mask_fire,
    scale_df = scale_df,
    paths = list(climate = p_climate, lcmap = p_lcmap, fire = p_fire,
                 herb_bio = p_herb_bio, gedi = p_gedi),
    params = list(lcmap_threshold = lcmap_threshold,
                  fire_threshold = fire_threshold,
                  pred_vars = pred_vars,
                  years = years,
                  cover_source = cover_source)
  )
}

load_fit <- function(suffix) {
  p <- file.path(paths$large, "Data_processed/BiomassQuantityData/Fit",
                 paste0("fitted_model_", suffix, ".rds"))
  stopifnot(file.exists(p))
  readRDS(p)
}

load_pred_raster <- function(suffix, cover_source = "rap") {
  # created in Fit/02_predict_rasters.R
  if (cover_source != "rap") {
    suffix <- paste0(suffix, "_cov-", cover_source)
  }
  p <- file.path(paths$large, "Data_processed/BiomassQuantityData/Predictions",
                 paste0("predicted_biomass_", suffix, ".tif"))
  stopifnot(file.exists(p))
  terra::rast(p)
}

#' Compare predicted biomass to BIGMAP AGB (treed areas)
#'
#' Loads the FIA BIGMAP AGB raster (aggregated to 1 km), masks both predicted
#' and BIGMAP rasters to treed areas (plus LCMAP and fire masks), samples
#' pixels, and computes summary metrics.
#'
#' @param r_pred Single-layer SpatRaster of predicted biomass (e.g., the
#'   tree cover-weighted layer).
#' @param rasters List returned by `load_conus_rasters()`, providing
#'   `zero_tree`, `mask_lcmap`, and `mask_fire`.
#' @param n_sample Integer; number of pixels to sample for scatter/metrics.
#' @param seed Integer; RNG seed for reproducible sampling.
#' @param root Character; root path for large files.
#'
#' @return A list with:
#'   - `metrics`: tibble from `cv_metrics_df()`
#'   - `sample`: data frame with columns `predicted` and `bigmap`
#'   - `r_pred_masked`: masked prediction raster
#'   - `r_bigmap_masked`: masked BIGMAP raster
compare_to_bigmap <- function(r_pred, rasters,
                              n_sample = 50000, seed = 8174,
                              root = paths$large) {

  p_bigmap <- file.path(root, "Data_processed", "BiomassQuantityData",
                        "BIGMAP_AGB-Mgha_2018_TOTAL_1000m.tif")
  stopifnot(file.exists(p_bigmap))
  r_bigmap <- terra::rast(p_bigmap)

  # align all rasters to r_pred extent (extents may differ across runs)
  r_bigmap <- terra::crop(r_bigmap, r_pred)
  mask_treed <- terra::crop(rasters$zero_tree, r_pred) == 0
  mask_lcmap <- terra::crop(rasters$mask_lcmap, r_pred)
  mask_fire <- terra::crop(rasters$mask_fire, r_pred)

  mask_all <- function(r) {
    r |>
      terra::mask(mask_treed, maskvalues = 0) |>
      terra::mask(mask_lcmap, maskvalues = 0) |>
      terra::mask(mask_fire, maskvalues = 0)
  }

  r_pred_masked <- mask_all(r_pred)
  r_bigmap_masked <- mask_all(r_bigmap)

  r_stack <- c(r_pred_masked, r_bigmap_masked)
  names(r_stack) <- c("predicted", "bigmap")

  set.seed(seed)
  df <- terra::spatSample(r_stack, size = n_sample,
                          method = "random", na.rm = TRUE, xy = FALSE)

  metrics <- cv_metrics_df(df, observed = "bigmap", predicted = "predicted")

  list(
    metrics = metrics,
    sample = df,
    r_pred_masked = r_pred_masked,
    r_bigmap_masked = r_bigmap_masked
  )
}

#' Load modelled (GLM) cover raster
#'
#' Loads a CONUS-wide raster of modelled vegetation cover by PFT.
#' Returns a SpatRaster with band names matching the RAP cover convention
#' (totalTreeCov, totalShrubCov, totalHerbaceousCov).
#'
#' @param cover_source Character; identifier for the cover model version
#'   (e.g., "glm_v01"). Used to construct the file path.
#' @param cover_cols Character vector of required band names. The output
#'   raster will have exactly these bands.
#' @param root Character; root path for large files.
#'
#' @return A SpatRaster with named layers matching `cover_cols`.
#' @export
load_modeled_cover <- function(root = paths$large) {
  
  p <- file.path(root, "Data_processed/CoverData/finalAbsoluteCoverData.tif")
  
  stopifnot(file.exists(p))
  
  r <- terra::rast(p)
  
  lookup <- c(
    "absTotalTree_CONUS" = "totalTreeCov",
    "absShrub_CONUS" = "totalShrubCov",
    "absTotalHerb_CONUS" = "totalHerbaceousCov"
  )
  
  stopifnot(names(lookup) %in% names(r))
  r <- r[[names(lookup)]]
  names(r) <- lookup[names(r)]
  r
}

load_cover <- function(cover_source = c('rap', 'model'),
                       root = paths$large) {
  cover_source = match.arg(cover_source)
  
  if(cover_source == 'model') {
    return(load_modeled_cover(root = root))
  } else {
    years = "2000-2023"
    p_cover_herb <- file.path(root, "Data_processed/CoverData/rap",
                              paste0("RAP_v3_cover_", years, "_1000m.tif"))
    # using woody cover from the gedi dataset period
    p_cover_woody <- file.path(root, "Data_processed/CoverData/rap",
                               paste0("RAP_v3_cover_", '2019-2023', "_1000m.tif"))
    
    stopifnot(file.exists(p_cover_herb, p_cover_woody))
    r_cover_herb <- terra::rast(p_cover_herb)$totalHerbaceousCov
    r_cover_woody <- terra::rast(p_cover_woody)[[c('totalTreeCov', 'totalShrubCov')]]
    r_cover <- c(r_cover_herb, r_cover_woody)/100 # RAP cover is %
  }
  r_cover
}

# raster showing areas we'd like to classify as having zero trees
# zeroTree is the binary layer (1 if no trees, 0 otherwise); fracZeroTree the
# fraction (0-1) of the cell's natural-land 30 m pixels with < 3% tree cover
read_zero_tree_raster <- function(layer = c('zeroTree', 'fracZeroTree'), 
                                  root = paths$large) {
  
  # file output by "Analysis/Cover/DataPrep/03_rap_frac-zero-tree.js"
  path <- file.path(paths$large, 'Data_processed/CoverData/rap',
                    'RAP_v3_fracZeroTree_lt3_2021-2021_1000m.tif')
  r <- terra::rast(path)
  layer <- match.arg(layer)
  lyr <- stringr::str_subset(names(r), layer)
  stopifnot(length(lyr) == 1)
  r <- r[[lyr]]
  names(r) <- layer
  r
}


#' Load EPA L3 ecoregion raster on the daymet grid
#'
#' @param root Root path for large files.
#' @return A single-layer SpatRaster with integer ecoregion IDs.
#' @export
load_ecoregion_raster <- function(
    epa_lev = NULL, 
    root = paths$large) {
  # file created in "Cover/DataPrep/01_rasterize_ecoregions.R"
  if (is.null(epa_lev)) epa_lev <- "L3"
  epa_lev <- match.arg(epa_lev, choices = c("L3", "L2"))
  p <- file.path(root, "Data_processed/regions",
                 paste0("EPA_", epa_lev, "_ecoregion_daymet_1000m.tif"))
  stopifnot(file.exists(p))
  r <- terra::rast(p)
  names(r) <- 'region'
  r
}

#' Load the EPA ecoregion ID lookup for the rasterized ecoregions
#'
#' `eco_id` is an arbitrary integer assigned when the shapefile was
#' rasterized, so it is only meaningful alongside this table. Always join it
#' rather than storing bare IDs downstream.
#'
#' @param epa_lev Ecoregion level, "L3" (default) or "L2".
#' @param root Root path for large files.
#' @return Tibble with eco_id, eco_code, eco_name.
#' @export
load_ecoregion_lookup <- function(epa_lev = NULL, root = paths$large) {
  # file created in "Cover/DataPrep/01_rasterize_ecoregions.R"
  if (is.null(epa_lev)) epa_lev <- "L3"
  epa_lev <- match.arg(epa_lev, choices = c("L3", "L2"))
  
  p <- file.path(root, "Data_processed/regions",
                 paste0("EPA_", epa_lev, "_ecoregion_lookup.csv"))
  stopifnot(file.exists(p))
  
  out <- readr::read_csv(p, show_col_types = FALSE)
  
  code_col <- if (epa_lev == "L3") "US_L3CODE" else "NA_L2CODE"
  name_col <- if (epa_lev == "L3") "US_L3NAME" else "NA_L2NAME"
  stopifnot(all(c("eco_id", code_col, name_col) %in% names(out)))
  
  out |>
    dplyr::select(eco_id,
                  eco_code = dplyr::all_of(code_col),
                  eco_name = dplyr::all_of(name_col)) |>
    dplyr::arrange(eco_id)
}

#' Load sagebrush biome mask raster
#'
#' Returns a SpatRaster with 1 inside the sagebrush biome, NA outside,
#' on the daymet 1km grid. Caches the reprojected result; recomputes
#' if the source file changes. Input file is from the holdrege et al 2024
#' data release (associated w/ fire ecology pub)
#'
#' @param path Path to the source raster.
#' @param cache_dir Directory for cached output.
#' @param force Logical; if TRUE, recompute even if cache exists.
#' @return A single-layer SpatRaster on the daymet grid.
#' @export
load_sagebrush_mask <- function(
    path = "../cheatgrass_fire/data_processed/data_publication/cell_number_ids.tif",
    cache_dir = file.path(paths$large, "Data_processed/cache"),
    force = FALSE) {
  
  if (!file.exists(path)) {
    stop("Source file not found: ", path,
         "\n  See load_sagebrush_mask() for details.")
  }
  
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  file_info <- file.info(path)
  fingerprint <- paste(file_info$size, file_info$mtime)
  
  cache_tif <- file.path(cache_dir, "sagebrush_mask_daymet.tif")
  cache_meta <- file.path(cache_dir, "sagebrush_mask_daymet.rds")
  
  if (!force && file.exists(cache_tif) && file.exists(cache_meta)) {
    cached <- readRDS(cache_meta)
    if (identical(cached$fingerprint, fingerprint)) {
      return(terra::rast(cache_tif))
    }
  }
  
  r1 <- terra::rast(path)
  # daymet_grid defined in mapping.R
  r1 <- terra::project(r1, daymet_grid, method = "near")
  r1 <- terra::trim(r1)
  r_mask <- terra::ifel(r1 > 0, 1, NA)
  names(r_mask) <- "sagebrush_mask"
  
  terra::writeRaster(r_mask, cache_tif, overwrite = TRUE, datatype = "INT1U")
  saveRDS(list(fingerprint = fingerprint), cache_meta)
  
  r_mask
}

read_mask <- function() {
  # file created in 
  # Analysis/Cover/00_create_snap_raster.R
  # this is the reference grid/mask to be used downstream
  # cell values are integers (or NA)
  terra::rast(file.path(paths$large, "Data_processed", 'masks', 
              "daymet_conus_snap_1000m.tif"))
}

#' Available-area mask: LCMAP and 20-yr unburned in a given year
#'
#' 1 where the cell passes the LCMAP mask (>= 90% keepable) and was >= 90%
#' unburned in the 20 years up to `year` (layer year_<year> of the per-year
#' fire mask); NA elsewhere. For maps and CONUS samples. RAP sampling (06)
#' uses the 2011-2023 mean fire mask instead.
#'
#' @param year Year of the fire mask layer.
#' @param root Root path for large files.
#' @return Single-layer SpatRaster on the snap grid.
read_available_mask <- function(year = 2023, root = paths$large) {
  mask_dir <- file.path(root, "Data_processed/masks")
  snap <- read_mask()
  lcmap <- terra::rast(file.path(mask_dir, "LCMAP_fracKeep_gte90_1000m.tif")) |>
    align_raster(snap)
  fire <- terra::rast(file.path(mask_dir,
                                "MTBS_fracUnburned_gte90_20yr_2000-2024_1000m.tif"))
  lyr <- paste0("year_", year)
  stopifnot(lyr %in% names(fire))
  fire <- align_raster(fire[[lyr]], snap)
  
  out <- terra::ifel(!is.na(snap) & lcmap == 1 & fire == 1, 1, NA)
  names(out) <- paste0("available_", year)
  out
}

#' Read the cover training data (cover, climate and soils by pixel-year)
#'
#' Reads the table written by Cover/DataPrep/08_combine_cover_and_covariates.R.
#' The _3yr columns are dropped. Climate columns are renamed to the short names
#' in `climate_name_lookup()`: the _CLIM suffix is dropped and _3yrAnom is kept
#' (e.g. tmean_meanAnnAvg_CLIM -> MAT, tmean_meanAnnAvg_3yrAnom ->
#' MAT_3yrAnom). AWHC becomes awc.
#'
#' @param vc Cover data version (e.g. "c01"), usually `opt$vc`.
#' @param group NULL for all rows, or one cover group, to keep only the rows
#'   where that group's response is observed (not NA): tree, shrub,
#'   herbaceous, bare_ground (cover, %); needle, broad (share of tree); forb,
#'   c3_grass, c4_grass (share of herbaceous).
#' @param normalize Logical; z-score the climate and soil columns with the
#'   CONUS-wide means and sds, and divide the anomalies by their sd across
#'   training pixel-years (see `read_scale_params()`). For versions without
#'   anomalies (c02), only the climate and soils parameters are used, and
#'   anomaly columns are an error.
#' @param as_sf Logical; return an sf data frame of cell-centre points in the
#'   CRS of `read_mask()`. The x and y columns are kept.
#' @param root Root path for large files.
#' @return Tibble, or sf data frame if `as_sf = TRUE`.
read_cover_training <- function(vc, group = NULL, normalize = FALSE,
                                as_sf = FALSE, root = paths$large) {
  group_cols <- c(tree = "cov_tree", shrub = "cov_shrub",
                  herbaceous = "cov_herbaceous",
                  bare_ground = "cov_bare_ground",
                  needle = "frac_needle", broad = "frac_broad",
                  forb = "frac_forb", c3_grass = "frac_c3",
                  c4_grass = "frac_c4")
  
  p <- file.path(root, "Data_processed/cover_combined",
                 paste0("cover_clim_soils_", vc, ".csv"))
  stopifnot(file.exists(p))
  if (!is.null(group)) group <- match.arg(group, names(group_cols))
  
  out <- readr::read_csv(p, show_col_types = FALSE) |>
    dplyr::select(-dplyr::ends_with("_3yr"))
  
  # short climate names, keeping the _3yrAnom suffix;
  # climate_name_lookup() errors on any name it doesn't know
  old <- stringr::str_subset(names(out), "_(CLIM|3yrAnom)$") |>
    setdiff("n_years_CLIM")
  long   <- stringr::str_replace(old, "_3yrAnom$", "_CLIM")
  suffix <- stringr::str_extract(old, "_3yrAnom$") |> tidyr::replace_na("")
  new    <- paste0(climate_name_lookup(long), suffix)
  names(out)[match(old, names(out))] <- new
  names(out)[names(out) == "AWHC"] <- climate_name_lookup("AWHC")
  
  if (normalize) {
    # versions without anomalies (c02: 2021 RAP-only sample, 05_rap2021_*.R)
    # have no anomaly scale params; use the shared climate and soils ones
    if (vc %in% c("c02")) {
      anom_cols <- stringr::str_subset(names(out), "_3yrAnom$")
      if (length(anom_cols) > 0) {
        stop(vc, " is not expected to have anomaly columns, but has: ",
             paste(anom_cols, collapse = ", "))
      }
      sp <- read_scale_params(root = root)
    } else {
      sp <- read_scale_params(vc, root = root)
    }
    soil_vars <- c("soilDepth", "clay_surface", "clay", "sand", "coarse",
                   "carbon", "awc")
    vars <- c(new, soil_vars)
    # every predictor must have parameters, so none is left unscaled
    stopifnot(all(vars %in% names(out)), all(vars %in% sp$variable))
    out <- standardize(out, vars = vars, scale_df = sp)$data
  }
  
  if (!is.null(group)) {
    out <- dplyr::filter(out, !is.na(.data[[group_cols[[group]]]]))
  }
  
  if (as_sf) {
    out <- sf::st_as_sf(out, coords = c("x", "y"),
                        crs = terra::crs(read_mask()), remove = FALSE)
  }
  out
}

#' Read the predictor scaling parameters
#'
#' Written by Cover/DataPrep/09_compute_scale_params.R. Climate normals and
#' soils: mean and sd across CONUS cells (shared by all models). Anomalies
#' (if `vc` is given): mean 0 and the sd across training pixel-years, so
#' scaling divides by the sd and keeps 0 as the normal.
#'
#' @param vc Cover data version for the anomaly parameters, or NULL for the
#'   climate and soils parameters only.
#' @param root Root path for large files.
#' @return Tibble with columns `variable`, `mean`, `sd` only (the `scale_df`
#'   format of `standardize()`, `standardize_raster()`, `unstandardize()` and
#'   `unstandardize_long()`). The CSVs also hold `n`, `source` and, for
#'   anomalies, `mean_obs`.
read_scale_params <- function(vc = NULL, root = paths$large) {
  dir <- file.path(root, "Data_processed/scale_params")
  p <- file.path(dir, "scale_params_climate_soils.csv")
  if (!is.null(vc)) {
    p <- c(p, file.path(dir, paste0("scale_params_anomalies_", vc, ".csv")))
  }
  stopifnot(file.exists(p))
  
  out <- purrr::map(p, readr::read_csv, show_col_types = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::select(variable, mean, sd)
  stopifnot(!anyDuplicated(out$variable))
  out
}


#' Read a fitted cover model
#'
#' @inheritParams cover_model_path
#' @return The list saved by the fitting script (e.g. fit, lambda, threshold,
#'   clusters, config, data).
#' @examples
#' cover_type <- "classification"
#' cover_model <- "forest"
#' vc <- "c01"
#' vmc <- "m01"
#' root <- paths$large
#' mod <- read_cover_model(cover_type = cover_type, cover_model = cover_model,
#'                         vc = vc, vmc = vmc, root = root)
read_cover_model <- function(cover_type, cover_model, vc, vmc,
                             root = paths$large) {
  p <- cover_model_path(cover_type, cover_model, vc, vmc, root = root)
  if (!file.exists(p)) stop("fitted model not found: ", p)
  readRDS(p)
}

#' Prediction raster of one component of a combined model
#'
#' @param component One row of a combined spec's `components` tibble
#'   (`cover_specs$combined`): `cover_type`, `cover_model`, `vc`, `vmc`.
#' @param scenario Climate scenario: "current", "BNU-ESM" or "IPSL-CM5A-MR".
#' @return SpatRaster, one layer: % cover, or % of the total for proportion
#'   models (the `cover` layer of the prediction raster).
#' @examples
#' component <- cover_specs$combined$herbaceous$m01.5$components[1, ]
#' scenario <- "current"
#' read_component_pred(component = component, scenario = scenario)
read_component_pred <- function(component, scenario) {
  stopifnot(nrow(component) == 1)
  path <- if (component$cover_type == "combined") {
    combined_path(component$cover_model, component$vmc, "prediction",
                  scenario)
  } else {
    cover_pred_path(component$cover_type, component$cover_model,
                    component$vc, component$vmc, scenario)
  }
  terra::rast(path)[["cover"]]
}

#' Needleleaf and broadleaf tree cover
#'
#' For one version of `cover_specs$combined$tree_all`: reads the combined tree
#' cover (04_combine_tree.R) and the needleleaf share (proportion/needle,
#' 03_predict_raster.R) for one scenario, and splits tree cover between
#' needleleaf (share) and broadleaf (100 - share). Computed on the fly, not
#' written to disk.
#'
#' @param vmc Version in `cover_specs$combined$tree_all`, e.g. "m01.5".
#' @param scenario Climate scenario: "current", "BNU-ESM" or "IPSL-CM5A-MR".
#' @param specs The `cover_specs` list (Functions/models/cover_specs.R).
#' @return SpatRaster with layers `share_needle`, `share_broad` (% of tree
#'   cover) and `cover_needle`, `cover_broad` (% cover).
#' @examples
#' vmc <- "m01.5"
#' scenario <- "current"
#' specs <- cover_specs
#' tree_leaf_cover(vmc = vmc, scenario = scenario, specs = specs)
tree_leaf_cover <- function(vmc, scenario, specs = cover_specs) {
  spec <- specs$combined$tree_all[[vmc]]
  if (is.null(spec)) stop("no spec for cover_specs$combined$tree_all$", vmc)
  comp <- spec$components
  
  tree <- read_component_pred(comp[comp$component == "tree", ], scenario)
  needle <- read_component_pred(comp[comp$component == "needle", ], scenario)
  shares <- c(needle, 100 - needle)
  names(shares) <- c("needle", "broad")
  
  # the two shares sum to 100 already, so share_sum is dropped
  out <- split_cover(shares, tree)
  out[[setdiff(names(out), "share_sum")]]
}

#' All final cover layers for one version of the combined "all" model
#'
#' For one version of `cover_specs$combined$all`: tree cover and its
#' needleleaf part (`tree_all`), shrub cover, and herbaceous cover and its
#' C3, C4 and forb parts (`herbaceous`), for one scenario. Read from the saved
#' predictions, or computed on the fly (needleleaf), so the layers can be used
#' in reports or figures without repeating the calculations.
#'
#' @param vmc Version in `cover_specs$combined$all`, e.g. "m01.5".
#' @param scenario Climate scenario: "current", "BNU-ESM" or "IPSL-CM5A-MR".
#' @param specs The `cover_specs` list (Functions/models/cover_specs.R).
#' @return SpatRaster, all in %: `cover_tree`, `cover_needle`,
#'   `share_needle` (% of tree cover), `cover_shrub`, `cover_herbaceous`,
#'   `share_c3`, `share_c4`, `share_forb` (% of herbaceous cover, summing to
#'   100), `cover_c3`, `cover_c4`, `cover_forb`.
#' @examples
#' vmc <- "m01.5"
#' scenario <- "current"
#' specs <- cover_specs
#' read_combined_all(vmc = vmc, scenario = scenario, specs = specs)
read_combined_all <- function(vmc, scenario, specs = cover_specs) {
  spec <- specs$combined$all[[vmc]]
  if (is.null(spec)) stop("no spec for cover_specs$combined$all$", vmc)
  comp <- spec$components
  tree_vmc <- comp$vmc[comp$component == "tree_all"]
  herb_vmc <- comp$vmc[comp$component == "herbaceous"]
  
  # total tree and herbaceous cover: the components the splits were made from
  tree_comp <- specs$combined$tree_all[[tree_vmc]]$components
  herb_comp <- specs$combined$herbaceous[[herb_vmc]]$components
  tree <- read_component_pred(tree_comp[tree_comp$component == "tree", ],
                              scenario)
  herb <- read_component_pred(herb_comp[herb_comp$component == "herbaceous", ],
                              scenario)
  shrub <- read_component_pred(comp[comp$component == "shrub", ], scenario)
  
  leaf <- tree_leaf_cover(tree_vmc, scenario, specs = specs)
  herb_types <- terra::rast(combined_path("herbaceous", herb_vmc, "prediction",
                                          scenario))
  
  names_leaf <- c("cover_needle", "share_needle", "cover_broad", "share_broad")
  names_herb_type <- c("share_c3", "share_c4", "share_forb",
                       "cover_c3", "cover_c4", "cover_forb")
  out <- c(tree, leaf[[names_leaf]], shrub, herb,
           herb_types[[names_herb_type]])
  names(out) <- c("cover_tree", names_leaf,
                  "cover_shrub", "cover_herbaceous",
                  names_herb_type)
  out
}


# paths -------------------------------------------------------------------



#' Paths for a combined tree cover prediction (04_combine_tree.R)
#'
#' No data version in the name: the component models can come from different
#' versions (see `cover_specs$combined$tree`).
#'
#' @param vmc Combination version, e.g. "m01".
#' @param what "prediction" (raster, one per scenario) or "blend" (csv of the
#'   component models and blending parameters).
#' @param scenario Climate scenario, for `what = "prediction"`.
#' @param root Root path for large files.
#' @return File path (character).
#' @examples
#' vmc <- "m01"
#' what <- "prediction"
#' scenario <- "current"
#' root <- paths$large
#' combined_tree_path(vmc = vmc, what = what, scenario = scenario, root = root)
combined_tree_path <- function(vmc, what = c("prediction", "blend"),
                               scenario = NULL, root = paths$large) {
  what <- match.arg(what)
  dir <- file.path(root, "Data_processed", "CoverData")
  if (what == "blend") {
    return(file.path(dir, "Fit", paste0("combined_tree_", vmc, "_blend.csv")))
  }
  stopifnot(!is.null(scenario))
  file.path(dir, "Predictions",
            paste0("combined_tree_", vmc, "_", scenario, ".tif"))
}

#' Path to a fitted cover model file
#'
#' The file name is defined only here, so the fitting, prediction and
#' diagnostics scripts can't drift apart.
#'
#' @param cover_type Model family: "classification", "cover" or "proportion".
#' @param cover_model Model within the family, e.g. "forest".
#' @param vc Cover data version, e.g. "c01".
#' @param vmc Cover model version, e.g. "m01".
#' @param root Root path for large files.
#' @return File path (character).
#' @examples
#' cover_type <- "classification"
#' cover_model <- "forest"
#' vc <- "c01"
#' vmc <- "m01"
#' root <- paths$large
#' cover_model_path(cover_type = cover_type, cover_model = cover_model,
#'                  vc = vc, vmc = vmc, root = root)
cover_model_path <- function(cover_type, cover_model, vc, vmc,
                             root = paths$large) {
  file.path(root, "Data_processed", "CoverData", "Fit",
            paste0(cover_type, "_", cover_model, "_", vc, "-", vmc, ".rds"))
}

#' Path to a cover model's prediction raster
#'
#' As written by 03_predict_raster.R; the file name is defined only here and
#' in that script.
#'
#' @inheritParams cover_model_path
#' @param scenario Climate scenario: "current", "BNU-ESM" or "IPSL-CM5A-MR".
#' @return File path (character).
#' @examples
#' cover_type <- "classification"
#' cover_model <- "forest"
#' vc <- "c02"
#' vmc <- "m01.0"
#' scenario <- "current"
#' root <- paths$large
#' cover_pred_path(cover_type = cover_type, cover_model = cover_model,
#'                 vc = vc, vmc = vmc, scenario = scenario, root = root)
cover_pred_path <- function(cover_type, cover_model, vc, vmc, scenario,
                            root = paths$large) {
  file.path(root, "Data_processed", "CoverData", "Predictions",
            paste0(cover_type, "_", cover_model, "_", vc, "-", vmc, "_",
                   scenario, ".tif"))
}

# Replaces combined_tree_path() in read_write.R. File names for tree are
# unchanged, so existing callers and outputs still work.

#' Paths for combined predictions (04_combine_*.R)
#'
#' No data version in the name: the component models can come from different
#' versions (see `cover_specs$combined`).
#'
#' @param cover_model Combined model, e.g. "tree" or "herbaceous".
#' @param vmc Combination version, e.g. "m01.5".
#' @param what "prediction" (raster, one per scenario), "blend" (tree: csv of
#'   the component models and blending parameters) or "components" (csv of
#'   the component models).
#' @param scenario Climate scenario, for `what = "prediction"`.
#' @param root Root path for large files.
#' @return File path (character).
#' @examples
#' cover_model <- "herbaceous"
#' vmc <- "m01.5"
#' what <- "prediction"
#' scenario <- "current"
#' root <- paths$large
#' combined_path(cover_model = cover_model, vmc = vmc, what = what,
#'               scenario = scenario, root = root)
combined_path <- function(cover_model, vmc,
                          what = c("prediction", "blend", "components"),
                          scenario = NULL, root = paths$large) {
  what <- match.arg(what)
  dir <- file.path(root, "Data_processed", "CoverData")
  stem <- paste0("combined_", cover_model, "_", vmc)
  if (what != "prediction") {
    return(file.path(dir, "Fit", paste0(stem, "_", what, ".csv")))
  }
  stopifnot(!is.null(scenario))
  file.path(dir, "Predictions", paste0(stem, "_", scenario, ".tif"))
}

# combined tree cover (04_combine_tree.R); see combined_path()
combined_tree_path <- function(vmc, what = c("prediction", "blend"),
                               scenario = NULL, root = paths$large) {
  combined_path("tree", vmc, match.arg(what), scenario, root)
}

# downloading files -----------------------------------------

#' Download a file from Drive if it is newer than the local copy
#'
#' @param drive_row One-row tibble from drive_ls with modifiedTime column.
#' @param local_dir Local directory to save to.
#' @return Invisible path to local file, or NULL if skipped.
download_if_newer <- function(drive_row, local_dir) {
  local_path <- file.path(local_dir, drive_row$name)
  
  if (file.exists(local_path)) {
    local_mtime <- file.mtime(local_path)
    drive_mtime <- as.POSIXct(drive_row$modifiedTime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
    if (!is.na(local_mtime) && local_mtime >= drive_mtime) {
      message("Skipping (local is current): ", drive_row$name)
      return(invisible(NULL))
    }
  }
  
  message("Downloading: ", drive_row$name)
  drive_download(file = drive_row$id, path = local_path, overwrite = TRUE)
  invisible(local_path)
}


#' Path to a fitted cover model file
#'
#' The file name is defined only here, so the fitting, prediction and
#' diagnostics scripts can't drift apart.
#'
#' @param cover_type Model family: "classification", "cover" or "proportion".
#' @param cover_model Model within the family, e.g. "forest".
#' @param vc Cover data version, e.g. "c01".
#' @param vmc Cover model version, e.g. "m01".
#' @param root Root path for large files.
#' @return File path (character).
#' @examples
#' cover_type <- "classification"
#' cover_model <- "forest"
#' vc <- "c01"
#' vmc <- "m01"
#' root <- paths$large
#' cover_model_path(cover_type = cover_type, cover_model = cover_model,
#'                  vc = vc, vmc = vmc, root = root)
cover_model_path <- function(cover_type, cover_model, vc, vmc,
                             root = paths$large) {
  file.path(root, "Data_processed", "CoverData", "Fit",
            paste0(cover_type, "_", cover_model, "_", vc, "-", vmc, ".rds"))
}

#' Read a fitted cover model
#'
#' @inheritParams cover_model_path
#' @return The list saved by the fitting script (e.g. fit, lambda, threshold,
#'   clusters, config, data).
#' @examples
#' cover_type <- "classification"
#' cover_model <- "forest"
#' vc <- "c01"
#' vmc <- "m01"
#' root <- paths$large
#' mod <- read_cover_model(cover_type = cover_type, cover_model = cover_model,
#'                         vc = vc, vmc = vmc, root = root)
read_cover_model <- function(cover_type, cover_model, vc, vmc,
                             root = paths$large) {
  p <- cover_model_path(cover_type, cover_model, vc, vmc, root = root)
  if (!file.exists(p)) stop("fitted model not found: ", p)
  readRDS(p)
}

# observed biomass datasets -----------------------------------------------

# aboveground biomass of live shrubs from Wright et al 2015
read_prepare_wright <- function() {
  
  dir <- file.path(paths$large, "Data_raw/BiomassDataSources",
                   "wright_sagebrush_fuels", "Data")
  p <- file.path(dir, "Pre-burn biomass.csv")
  dat <- read_csv(p)
  
  cover <- read_csv(file.path(dir, "Cover and Height.csv"))
  
  cover2 <- cover |> 
    # Pre-burn vegetation cover: sage, other shrub
    rename(shrubCover = p_sa_os) |> 
    select(siteid, shrubCover)
  
  
  dat2 <- dat |> 
    # Pre-burn biomass of live other shrub
    # plus Pre-burn biomass of live sage
    # units are Mg/ha
    mutate(shrubBiomass = b_los + b_ls) |> 
    select(siteid, shrubBiomass) |> 
    left_join(cover2, by = 'siteid') |>
    mutate(shrubBiomass_potential = shrubBiomass/(shrubCover/100))
  
  site <- read_csv(file.path(dir, "Site variables.csv"))
  
  # Merge with site data
  dat3 <- site |> 
    select(siteID, lat, long) |> 
    # source CSV stores western longitudes as positive values
    mutate(long = -abs(long)) |> 
    left_join(dat2, by = c("siteID" = "siteid")) |> 
    sf::st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE) |> 
    rename(site = siteID) |> 
    select(-long, -lat)
  
  dat3
}


#' Read and prepare Sevilleta SEV-182 shrub biomass data
#'
#' Loads the per-quadrat biomass CSV, filters to shrub species, and joins
#' hardcoded site coordinates extracted from the package EML metadata.
#'
#' @return An sf data frame of shrub biomass observations with site
#'   coordinates (WGS84).
#' @noRd
read_prepare_sevilleta <- function() {
  dir <- file.path(paths$large, "Data_raw/BiomassDataSources",
                   "sevilleta_agb")
  
  dat <- read_csv(file.path(dir, 'sev182_NPP_core_biomass.csv'))
  
  # site coordinates copied from package EML metadata
  # (geographicCoverage blocks in eml_metadata.xml)
  sites <- tibble::tribble(
    ~site,           ~lon,      ~lat,
    "core_black",    -106.736,  34.3331,
    "core_blue",     -106.631,  34.3348,
    "core_creosote", -106.736,  34.3331
  )
  
  shrub_biomass <- dat |>
    filter(FunctionalGroup == 'shrub') |>
    # sum across shrub species within each quadrat-season
    group_by(site, year, season, web, plot, quad) |>
    summarise(biomass_g_m2 = sum(biomass.BM, na.rm = TRUE), 
              cover = sum(cover, na.rm = TRUE),
              .groups = "drop") |>
    # take peak (max across seasons) within each quadrat-year
    group_by(site, year, web, plot, quad) |>
    summarise(peak_g_m2 = max(biomass_g_m2, na.rm = TRUE),
              cover = max(cover, na.rm = TRUE), 
              .groups = "drop") |>
    # average across quadrats within site-year
    group_by(site, year) |>
    summarise(
      # this file indicates quadrats are 1m x 1m:
      # https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sev.129.258322&contentType=application/xml
      # so units are g/m2
      peak_g_m2  = mean(peak_g_m2, na.rm = TRUE),
      peak_Mg_ha = peak_g_m2 * 0.01,
      shrubCover = mean(cover, na.rm = TRUE),
      n_quadrats = n(),
      .groups    = "drop"
    ) |> 
    rename(shrubBiomass = peak_Mg_ha) |> 
    mutate(shrubBiomass_potential = shrubBiomass/(shrubCover/100)) |> 
    select(-peak_g_m2) |> 
    left_join(sites, by = "site") |> 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
    select(-lat, -lon)
  shrub_biomass
}


#' Read and prepare Jornada Basin LTER NPP shrub biomass data
#'
#' Loads the per-quadrat biomass CSV (knb-lter-jrn.210011001), filters to
#' shrub and sub-shrub species (`form %in% c("SHRUB", "S-SHR")`), and
#' aggregates to site-year means following the same logic as
#' `read_prepare_sevilleta()`. Per-site coordinates are not available in
#' the package metadata; a single centroid of the bounding box is used for
#' all sites.
#'
#' @return An sf data frame of shrub biomass observations with centroid
#'   coordinates (WGS84).
read_prepare_jornada <- function() {
  dir <- file.path(paths$large, "Data_raw/BiomassDataSources",
                   "jornada_agb")

  dat <- read_csv(file.path(dir, "quadrat_plant_volume_and_biomass_data.csv"))

  # bounding box from EML metadata (eml_metadata.xml):
  # west = -106.865, east = -106.713, south = 32.488, north = 32.669
  # Per-site coordinates unavailable; use centroid for all sites
  bbox_lon <- mean(c(-106.865, -106.713))
  bbox_lat <- mean(c(32.488, 32.669))

  shrub_biomass <- dat |>
    filter(form %in% c("SHRUB", "S-SHR")) |>
    # sum across shrub species within each quadrat-season
    group_by(site, zone, year, season, quad) |>
    summarise(biomass_g_m2 = sum(biomass, na.rm = TRUE),
              cover = sum(cum_cover, na.rm = TRUE),
              .groups = "drop") |>
    # take peak (max across seasons) within each quadrat-year
    group_by(site, zone, year, quad) |>
    summarise(peak_g_m2 = max(biomass_g_m2, na.rm = TRUE),
              cover = max(cover, na.rm = TRUE),
              .groups = "drop") |>
    # average across quadrats within site-year
    group_by(site, zone, year) |>
    summarise(
      # quadrats are 1m x 1m so biomass is g/m2
      peak_g_m2  = mean(peak_g_m2, na.rm = TRUE),
      peak_Mg_ha = peak_g_m2 * 0.01,
      shrubCover = mean(cover, na.rm = TRUE),
      n_quadrats = n(),
      .groups    = "drop"
    ) |>
    rename(shrubBiomass = peak_Mg_ha) |>
    mutate(shrubBiomass_potential = shrubBiomass / (shrubCover / 100)) |>
    select(-peak_g_m2) |>
    mutate(lon = bbox_lon, lat = bbox_lat) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
    select(-lon, -lat)

  shrub_biomass
}

# read in field data of shrub biomass observations
read_shrub_biomass <- function() {
  dat_l <- list(wright = read_prepare_wright(),
       jornada = read_prepare_jornada(),
       sevilleta = read_prepare_sevilleta())
  
  dat <- map(dat_l, sf::st_transform, crs = crs_daymet) |> 
    bind_rows(.id = 'study')
  dat
}

