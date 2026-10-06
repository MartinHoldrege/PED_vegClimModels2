/*
Fraction of natural-land 30m pixels per daymet cell 
that have <x% tree cover. Cells are masked with the 1km fire mask
(02_fire_fracUnburned.js) and LCMAP mask, as in 03_rap_sample.js.

Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStartRap = 2021;
var yearEndRap = 2021;

var cutoffs = [3];
var zeroTreeFracCutoff = 0.9; // cuttoff for faction of 1km 

var export2Drive = true;
var export2Asset = false;
// 1km masks
var cellKeep = fg.lcmapMaskBinary().and(fg.fireMaskYear(yearEndRap));

for (var i = 0; i < cutoffs.length; i++) {
  var cutoff = cutoffs[i];
  // read in data -------------------------------------
  // created in 03_rap_no-tree.js (masked by LCMAP)
  var notForest30 = ee.Image(fg.pathAsset + 'rap/RAP_v3_tree-lt' + cutoff + '_masked_' +
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
    .updateMask(cellKeep)
    .rename('fracNotForest');
    
  var zeroTree = fracNotForest
    .gte(zeroTreeFracCutoff)
    .rename('zeroTree_' + zeroTreeFracCutoff*100);
  
  fracNotForest = fracNotForest.addBands(zeroTree);
  
  // visualize ----------------------------------------
  Map.addLayer(fracNotForest, {min: 0, max: 1, palette: ['white', 'black']},
    'frac not forest (unburned, natural land)', false);
  Map.addLayer(zeroTree.selfMask(), {palette: ['black']}, 'zero tree area', false);
  // export -------------------------------------------
  var fileName = 'RAP_v3_fracNotForest_lt' + cutoff + '_' +
    yearStartRap + '-' + yearEndRap + fg.resLabel;

  if(export2Drive) {
    fg.exportDrive(fracNotForest, fileName, 'PED_vegClimModels2');
  }
  
  if(export2Asset) {
    fg.exportAsset(fracNotForest, fileName, fg.pathAsset + 'rap/');
  }
  
}
