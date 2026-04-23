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
    path = file.path(paths$large, "Data_processed/BiomassQuantityData/", "
                     DayMetData_allCONUS_2023ClimateValues_raster.tif")
    ) {
  r <- terra::rast(path)
  names(r) <- climate_name_lookup(names(r))
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
                               years = "2000-2023",
                               force_scale = FALSE) {
  
  # file paths
  p_climate <- file.path(root, "Data_processed/BiomassQuantityData",
                         "DayMetData_allCONUS_2023ClimateValues_raster.tif")
  p_lcmap <- file.path(root, "Data_processed/masks",
                       "LCMAP_fracKeep_1000m.tif")
  p_fire <- file.path(root, "Data_processed/masks",
                      paste0("MTBS_fracUnburned_", years, "_1000m.tif"))
  p_cover <- file.path(root, "Data_processed/CoverData/rap",
                       paste0("RAP_v3_cover_", years, "_1000m.tif"))
  p_herb_bio <- file.path(root, "Data_processed/BiomassQuantityData/rap",
                          paste0("RAP_v3_herbaceousAGB_mask-Lcmap",
                                 lcmap_threshold * 100, "_", years, "_1000m.tif"))
  p_gedi <- file.path(root, "Data_processed/BiomassQuantityData",
                      "GEDI_biomassRaster_onDayMetGrid.tif")
  
  all_paths <- c(p_climate, p_lcmap, p_fire, p_cover, p_herb_bio, p_gedi)
  missing <- all_paths[!file.exists(all_paths)]
  if (length(missing) > 0) {
    stop("Missing raster files:\n  ", paste(missing, collapse = "\n  "))
  }
  
  cat("Loading CONUS rasters...\n")
  
  r_climate <- read_climate_raster(p_climate)
  r_lcmap <- terra::rast(p_lcmap)
  r_fire <- terra::rast(p_fire)
  r_cover <- terra::rast(p_cover)/100 # RAP cover is %
  
  r_herb_bio <- terra::rast(p_herb_bio)
  names(r_herb_bio) <- "totalHerbaceousBio"
  
  r_gedi <- terra::rast(p_gedi)
  names(r_gedi) <- "totalWoodyBio"
  
  # align
  aligned <- align_raster_extents(rast_list = list(
    climate = r_climate, lcmap = r_lcmap, fire = r_fire,
    cover = r_cover, herb_bio = r_herb_bio, gedi = r_gedi
  ))
  
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
    mask_lcmap = mask_lcmap,
    mask_fire = mask_fire,
    scale_df = scale_df,
    paths = list(climate = p_climate, lcmap = p_lcmap, fire = p_fire,
                 cover = p_cover, herb_bio = p_herb_bio, gedi = p_gedi),
    params = list(lcmap_threshold = lcmap_threshold,
                  fire_threshold = fire_threshold,
                  pred_vars = pred_vars,
                  years = years)
  )
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

