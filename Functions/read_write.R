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

# renames layers to shorter names
read_climate_raster <- function(
    path = file.path(paths$large, "Data_processed/BiomassQuantityData", 
                     "DayMetData_allCONUS_2023ClimateValues_raster.tif"),
    
    # created in 01_soils_calculate_awc
    path_soil = file.path(paths$large, 
              "./Data_processed/soils/", 
              "awc_SOLUS100_1000m.tif")
    ) {
  r <- terra::rast(path)
  names(r) <- climate_name_lookup(names(r))
  
  if(!is.null(path_soil)) {
    r_soil <- terra::rast(path_soil)
    names(r_soil) <- stringr::str_replace(names(r_soil), '_cm', '')
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
                               force_scale = FALSE) {

  years = "2000-2023"
  # file paths
  p_climate <- file.path(root, "Data_processed/BiomassQuantityData",
                         "DayMetData_allCONUS_2023ClimateValues_raster.tif")
  p_lcmap <- file.path(root, "Data_processed/masks",
                       "LCMAP_fracKeep_1000m.tif")
  p_fire <- file.path(root, "Data_processed/masks",
                      paste0("MTBS_fracUnburned_", years, "_1000m.tif"))
  p_cover_herb <- file.path(root, "Data_processed/CoverData/rap",
                       paste0("RAP_v3_cover_", years, "_1000m.tif"))
  # using woody cover from the gedi dataset period
  p_cover_woody <- file.path(root, "Data_processed/CoverData/rap",
                            paste0("RAP_v3_cover_", '2019-2023', "_1000m.tif"))
  p_herb_bio <- file.path(root, "Data_processed/BiomassQuantityData/rap",
                          paste0("RAP_v3_herbaceousAGB_mask-Lcmap",
                                 lcmap_threshold * 100, "_", years, "_1000m.tif"))
  p_gedi <- file.path(root, "Data_processed/BiomassQuantityData",
                      "GEDI_biomassRaster_onDayMetGrid.tif")
  p_fracNotForest <- file.path(root, "Data_processed/CoverData/rap",
    paste0("RAP_v3_fracNotForest_lt3_mask-lcmap50_2019-2023_1000m.tif"))
  
  
  all_paths <- c(p_climate, p_lcmap, p_fire, p_cover_herb, p_cover_woody, 
                 p_herb_bio, p_gedi, p_fracNotForest)
  
  missing <- all_paths[!file.exists(all_paths)]
  if (length(missing) > 0) {
    stop("Missing raster files:\n  ", paste(missing, collapse = "\n  "))
  }
  
  cat("Loading CONUS rasters...\n")
  
  r_climate <- read_climate_raster(p_climate)
  r_lcmap <- terra::rast(p_lcmap)
  r_fire <- terra::rast(p_fire)
  r_cover_herb <- terra::rast(p_cover_herb)$totalHerbaceousCov
  r_cover_woody <- terra::rast(p_cover_woody)[[c('totalTreeCov', 'totalShrubCov')]]
  r_cover <- c(r_cover_herb, r_cover_woody)/100 # RAP cover is %
  r_region <- load_ecoregion_raster()
  
  r_herb_bio <- terra::rast(p_herb_bio)
  names(r_herb_bio) <- "totalHerbaceousBio"
  
  r_gedi <- terra::rast(p_gedi)
  names(r_gedi) <- "totalWoodyBio"
  
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
                 cover_herb = p_cover_herb,
                 cover_woody = p_cover_woody, 
                 herb_bio = p_herb_bio, gedi = p_gedi),
    params = list(lcmap_threshold = lcmap_threshold,
                  fire_threshold = fire_threshold,
                  pred_vars = pred_vars,
                  years = years)
  )
}

load_fit <- function(suffix) {
  p <- file.path(paths$large, "Data_processed/BiomassQuantityData/Fit",
                 paste0("fitted_model_", suffix, ".rds"))
  stopifnot(file.exists(p))
  readRDS(p)
}

load_pred_raster <- function(suffix) {
  # created in Fit/02_predict_rasters.R
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

  mask_treed <- rasters$zero_tree == 0

  mask_all <- function(r) {
    r |>
      terra::mask(mask_treed, maskvalues = 0) |>
      terra::mask(rasters$mask_lcmap, maskvalues = 0) |>
      terra::mask(rasters$mask_fire, maskvalues = 0)
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
load_modelled_cover <- function(cover_source,
                                cover_cols,
                                root = paths$large) {
  
  p <- file.path(root,
                 "Data_processed/CoverData/modelled",
                 paste0("modelled_cover_", cover_source, ".tif"))
  
  if (file.exists(p)) {
    r <- terra::rast(p)
    
    # rename bands if needed to match RAP convention
    # (assumes same order if names differ)
    if (!all(cover_cols %in% names(r))) {
      if (terra::nlyr(r) == length(cover_cols)) {
        stop("need to implement renaming code")
        # names(r) <- cover_cols
      } else {
        stop("Modelled cover raster has ", terra::nlyr(r),
             " layers but ", length(cover_cols), " cover_cols specified.")
      }
    } else {
      r <- r[[cover_cols]]
    }
    
    return(r)
    
  } else {
    # PLACEHOLDER: generate fake modelled cover from RAP + noise
    message("Modelled cover file not found: ", p,
            "\n  Using noised RAP cover as placeholder.")
    
    r_rap <- load_conus_rasters(pred_vars = "MAT")$cover[cover_cols]
    
    set.seed(123)
    r_out <- terra::rast(lapply(cover_cols, function(col) {
      r_layer <- r_rap[[col]]
      # add spatially correlated noise (~10% of value)
      noise <- r_layer * 0.1 * (terra::init(r_layer, fun = runif) - 0.5) * 2
      r_noisy <- r_layer + noise
      # clamp to valid range
      r_noisy <- terra::clamp(r_noisy, lower = 0, upper = 100)
      names(r_noisy) <- col
      r_noisy
    }))
    
    return(r_out)
  }
}

#' Load EPA L3 ecoregion raster on the daymet grid
#'
#' @param root Root path for large files.
#' @return A single-layer SpatRaster with integer ecoregion IDs.
#' @export
load_ecoregion_raster <- function(
    epa_lev = c('L3', 'L2'), 
    root = paths$large) {
  # file created in "DataPrep/01_rasterize_ecoregions.R"
  lev <- arg.match(epa_lev)
  p <- file.path(root, "Data_processed/regions",
                 paste0("EPA_", lev, "_ecoregion_daymet_1000m.tif"))
  stopifnot(file.exists(p))
  r <- terra::rast(p)
  names(r) <- 'region'
  r
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
