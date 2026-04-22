/*
Fraction of unburned, natural-land 30m pixels per daymet cell 
that have <1% tree cover.

Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStartRap = 2019;
var yearEndRap = 2023;

// read in data -------------------------------------
// created in 03_rap_notForest.js (masked by fire + LCMAP)
var notForest30 = ee.Image(fg.pathAsset + 'rap/RAP_v3_tree-lt1_masked_' +
  yearStartRap + '-' + yearEndRap + '_30m');

// process ------------------------------------------
var fracNotForest = notForest30
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 2e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform
  })
  .rename('fracNotForest');

// visualize ----------------------------------------
Map.addLayer(fracNotForest, {min: 0, max: 1, palette: ['white', 'black']},
  'frac not forest (unburned, natural land)', false);

// export -------------------------------------------
var fileName = 'RAP_v3_fracNotForest_' +
  yearStartRap + '-' + yearEndRap + fg.resLabel;

Export.image.toAsset({
  image: fracNotForest,
  description: fileName,
  assetId: fg.pathAsset + 'rap/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});