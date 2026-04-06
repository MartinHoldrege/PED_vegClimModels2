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

var rap = ee.ImageCollection('projects/rangeland-analysis-platform/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));

// created in 01_create_lcmap_mask
var fracKeep = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/LCMAP_fracKeep_daymet');
print(fracKeep.geometry())
// process ------------------------------------------

print(rap.first())
// mean cover across years
var rapMean = rap
  .select(['TRE', 'SHR', 'AFG', 'PFG'])
  .mean()
  .setDefaultProjection(rap.first().projection());

// bands of interest
var tree = rapMean.select('TRE').rename('totalTreeCov');
var shrub = rapMean.select('SHR').rename('totalShrubCov');
var totalHerbaceous = rapMean.select('AFG').add(rapMean.select('PFG')).rename('totalHerbaceousCov');

var coverStack = tree
  .addBands(shrub)
  .addBands(totalHerbaceous);

// apply LCMAP mask
var mask = fracKeep.gte(fracKeepThreshold);
var coverMasked = coverStack.updateMask(mask);

// reproject to daymet grid
var coverDaymet = coverMasked
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 2e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform
  });

// visualize ----------------------------------------
 
// for speed 
Map.addLayer(fracKeep.gte(fracKeepThreshold), {min: 0, max: 1, palette: ['white', 'black']},
  'mask, ' + fracKeepThreshold);
Map.addLayer(coverMasked.select('totalTreeCov'), {min: 0, max: 60, palette: ['white', 'darkgreen']},
  'tree cover');
Map.addLayer(coverMasked.select('totalShrubCov'), {min: 0, max: 60, palette: ['white', 'brown']},
  'shrub cover');
Map.addLayer(coverMasked.select('totalHerbaceousCov'), {min: 0, max: 60, palette: ['white', 'gold']},
  'herbaceous cover');

// exploration --------------------------------------

var herbThreshold = 5;
var woodyCutoff = 5;

var herbDominant = coverDaymet.select('totalTreeCov').lt(woodyCutoff)
  .and(coverDaymet.select('totalShrubCov').lt(woodyCutoff))
  .and(coverDaymet.select('totalHerbaceousCov').gte(herbThreshold))
  .selfMask();

Map.addLayer(herbDominant, {palette: ['lime']}, 'herb-dominant pixels');


// export -------------------------------------------
var scale = Math.abs(fg.crsTransform[0]);
print(scale)
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