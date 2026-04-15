/*
fraction of daymet cell that has burned in given time period
Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------

var yearStart = 2000;
var yearEnd = 2023;

// created in 01_mtbs_everburned.js
var everBurned = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/fire/MTBS_everBurned_30m_' + yearStart + '-' + yearEnd);

//  proces -----
var fracUnburned = everBurned.not()
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 2e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform
  })
  .rename('fracUnburned');

// export -------------------------------------------

var fileName = 'MTBS_fracUnburned_' + yearStart + '-' + yearEnd + fg.resLabel;

Export.image.toAsset({
  image: fracUnburned,
  description: fileName,
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/fire/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});