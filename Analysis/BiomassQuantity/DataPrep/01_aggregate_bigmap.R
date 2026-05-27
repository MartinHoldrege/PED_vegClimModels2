# Description: take the FIA bigmap dataset and average to the daymet grid
# started May 27, 2026
source('Functions/init.R')
source_functions()

r_agb_30m <- terra::rast(file.path(paths$large, "Data_raw", "BiomassDataSources",
                                   "BIGMAP_AGB_2018_SPCD0000_TOTAL",
                                   "Hosted_AGB_0000_2018_TOTAL_11172024101136.tif"))

r_agb_tac <- terra::project(r_agb_30m, daymet_grid, method = "average") # tons per acre

r_agb <- r_agb_tac * 2.2417 # to Mg/ha

filename <- file.path(paths$large, 'Data_processed', 'BiomassQuantityData',
                      'BIGMAP_AGB-Mgha_2018_TOTAL_1000m.tif')

terra::writeRaster(r_agb, filename = filename)