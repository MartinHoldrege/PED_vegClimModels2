/*
Export GeoTIFFs of cover, biomass, and mask assets created in upstream scripts
to Google Drive for downstream use in R modeling pipeline.

Input assets created in:
  - LCMAP fracKeep: 01_create_lcmap_mask.js
  - MTBS fracUnburned: 01_mtbs_everBurned.js, 02_mtbs_fracUnburned.js
  - RAP cover: 02_rap_mean.js
  - RAP herbaceous biomass: 02_rap_mean.js

Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------

var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------

var yearStart = 2000;
var yearEnd = 2023;
var maskCutoffLcmap = 0.9;
var maskCutoffFire = 0.9;
var driveFolder = 'PED_vegClimModels2';

// export toggles
var exportLcmapMask = false;
var exportFireMask = false;
var exportRapCover = false;
var exportRapBiomass = false;
var exportfracNotForest = true

// LCMAP mask (fracKeep + binary) -------------------

  // created in 01_create_lcmap_mask.js
  var fracKeep = ee.Image(fg.pathAsset + 'LCMAP_fracKeep_daymet')
    .rename('fracKeep');
  var fracKeepBinary = fracKeep
    .gte(maskCutoffLcmap)
    .selfMask()
    .rename('fracKeep_gte' + maskCutoffLcmap*100);
    
  var lcmapStack = fracKeep.addBands(fracKeepBinary).toFloat();
if (exportLcmapMask) {
  var lcmapFileName = 'LCMAP_fracKeep' + fg.resLabel;
  Export.image.toDrive({
    image: lcmapStack,
    description: lcmapFileName,
    folder: driveFolder,
    fileNamePrefix: lcmapFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });
}

// Fire mask (fracUnburned + binary) ----------------

  // created in 01_mtbs_everBurned.js, 02_mtbs_fracUnburned.js
  var fracUnburned = ee.Image(fg.pathAsset + 'fire/MTBS_fracUnburned_' +
    yearStart + '-' + yearEnd + fg.resLabel)
    .rename('fracUnburned');
  var fracUnburnedBinary = fracUnburned
    .gte(maskCutoffFire)
    .selfMask()
    .rename('fracUnburned_gte' + maskCutoffFire*100);
  
  var fireStack = fracUnburned.addBands(fracUnburnedBinary).toFloat();

if (exportFireMask) {
  var fireFileName = 'MTBS_fracUnburned_' + yearStart + '-' + yearEnd + fg.resLabel;
  Export.image.toDrive({
    image: fireStack,
    description: fireFileName,
    folder: driveFolder,
    fileNamePrefix: fireFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });
}

// RAP cover ----------------------------------------
if (exportRapCover) {
  // created in 02_rap-cover_mean.js
  var rapCoverFileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd + fg.resLabel;
  var rapCover = ee.Image(fg.pathAsset + 'rap/' + rapCoverFileName);

  Export.image.toDrive({
    image: rapCover,
    description: rapCoverFileName,
    folder: driveFolder,
    fileNamePrefix: rapCoverFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });
}

// RAP herbaceous biomass ---------------------------
if (exportRapBiomass) {
  // created in 02_rap-cover_mean.js
  var rapBioAssetName = 'RAP_v3_herbaceousAGB_' + yearStart + '-' + yearEnd + fg.resLabel;
  var rapBioFileName = 'RAP_v3_herbaceousAGB_mask-Lcmap' + maskCutoffLcmap * 100 + '_' +
      yearStart + '-' + yearEnd + fg.resLabel;
  var rapBio = ee.Image(fg.pathAsset + 'rap/' + rapBioAssetName)
      .updateMask(fracKeepBinary);

  Export.image.toDrive({
    image: rapBio,
    description: rapBioFileName,
    folder: driveFolder,
    fileNamePrefix: rapBioFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });

}

if(exportfracNotForest) {
  var notForestAssetName = 'RAP_v3_fracNotForest_2019-2023' + fg.resLabel;
  var notForestFileName = 'RAP_v3_fracNotForest_mask-lcmap' + maskCutoffLcmap * 100 
    + '-fire' + maskCutoffFire*100 + '_2019-2023' + fg.resLabel;
    
  // created in 03_rap_fracNotForest.js
  var notForest =   ee.Image(fg.pathAsset + 'rap/' + notForestAssetName)
      .updateMask(fracKeepBinary)
      .updateMask(fracUnburnedBinary);
      
  Map.addLayer(notForest, {min: 0, max: 1, palette: 'white,black'}, 'fracNotForest', false);
      
  Export.image.toDrive({
    image: notForest,
    description: notForestFileName,
    folder: driveFolder,
    fileNamePrefix: notForestFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });
}



