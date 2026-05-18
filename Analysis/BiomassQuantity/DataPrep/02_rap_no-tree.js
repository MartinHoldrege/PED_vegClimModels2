/*
Determine which pixels have no or close to no tree cover (based on RAP),
masked to exclude burned (MTBS 2000-2023) and developed/cropland/water
(LCMAP 2021) pixels.

Author: Martin Holdrege
Started: April 6, 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStartRap = 2019;
var yearEndRap = 2023;
var yearStartFire = 2000;
var yearEndFire = 2023;

var cutoffs = [3, 5, 10];
// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStartRap, yearEndRap, 'year'));

// created in 01_mtbs_everBurned.js
var everBurned = ee.Image(fg.pathAsset + 'fire/MTBS_everBurned_30m_' +
  yearStartFire + '-' + yearEndFire)
  .unmask(0);

// process ------------------------------------------
var rapMean = rap.select(['TRE']).mean();

for (var i = 0; i < cutoffs.length; i++) {
    
  var cutoff = cutoffs[i];
  var notForest = rapMean
    .lt(ee.Image(cutoff)); // <1% cover
  
  // mask: unburned AND natural land
  var keepMask = everBurned.not().and(fg.lcmapMask);
  var notForestMasked = notForest
    .updateMask(keepMask)
    .toByte()
    .rename('TRE_lt' + cutoff);
  
  // visualize ----------------------------------------
  Map.addLayer(rapMean, {min: 0, max: 5, palette: 'white,green'}, 'tree cov', false);
  Map.addLayer(notForestMasked, {min: 0, max: 1, palette: 'white,black'}, '<1% trees (masked)', false);
  
  // export -------------------------------------------
  var fileName = 'RAP_v3_tree-lt' + cutoff + '_masked_' +
    yearStartRap + '-' + yearEndRap + '_30m';
  
  var policy = {};
  policy['TRE_lt' + cutoff] = 'mode';
  
  Export.image.toAsset({
    image: notForestMasked,
    description: fileName,
    assetId: fg.pathAsset + 'rap/' + fileName,
    crs: fg.crs,
    scale: 30,
    region: fg.region,
    maxPixels: 1e12,
    pyramidingPolicy: policy
  });
}