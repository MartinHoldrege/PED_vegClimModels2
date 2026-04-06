/*

get rap cover, for the years of the biomass datasets, so can
filter pure herbceous pixels for fitting the intercepts to those data

Author: Martin Holdrege

Started: April 6, 2026

*/

// dependencies -------------------------------------

var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStart = 2019; // GEDI v2.1 data collected from 2019--2023
var yearEnd = 2023;
var fracKeepThreshold = 0.9;

// read in data -------------------------------------

var fileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd  + fg.resLabel;
var cover = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + fileName);

// created in 01_create_lcmap_mask
var fracKeep = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/LCMAP_fracKeep_daymet');

// process ------------------------------------------


// apply LCMAP mask
var mask = fracKeep.gte(fracKeepThreshold);
var coverMasked = cover.updateMask(mask);


// visualize ----------------------------------------
 
// for speed 
Map.addLayer(fracKeep.gte(fracKeepThreshold), {min: 0, max: 1, palette: ['white', 'black']},
  'mask, ' + fracKeepThreshold, false);
Map.addLayer(coverMasked.select('totalTreeCov'), {min: 0, max: 60, palette: ['white', 'darkgreen']},
  'tree cover', false);
Map.addLayer(coverMasked.select('totalShrubCov'), {min: 0, max: 60, palette: ['white', 'brown']},
  'shrub cover', false);
Map.addLayer(coverMasked.select('totalHerbaceousCov'), {min: 0, max: 60, palette: ['white', 'gold']},
  'herbaceous cover', false);

// exploration --------------------------------------

var herbThreshold = 10;
var treeCutoff = 2;
var shrubCutoff = 5;

var herbDominant = coverMasked.select('totalTreeCov').lt(treeCutoff)
  .and(coverMasked.select('totalShrubCov').lt(shrubCutoff))
  .and(coverMasked.select('totalHerbaceousCov').gte(herbThreshold))
  .selfMask();

Map.addLayer(herbDominant, {palette: ['lime']}, 'herb-dominant pixels');


// export -------------------------------------------
var maskName = 'mask-' + fracKeepThreshold;
var fileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd + '_' + maskName + '_daymet';

// Export.image.toDrive({
//   image: coverDaymet,
//   description: fileName,
//   folder: 'PED_vegClimModels2',
//   fileNamePrefix: fileName,
//   crs: fg.crs,
//   crsTransform: fg.crsTransform,
//   region: fracKeep.geometry(),
//   maxPixels: 1e12,
//   fileFormat: 'GeoTIFF'
// });