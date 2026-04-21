/*
determine which pixels have no or close to no tree cover (based on RAP)

Author: Martin Holdrege
Started: April 6, 2026
*/

// params -------------------------------------------

var exportCov = false;
var exportBio = false; 

var yearStart = 2019; // years corresponding to cover model training data
var yearEnd = 2023;

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');



// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));

// process ------------------------------------------
var rapMean = rap
  .select(['TRE'])
  .mean();


Map.addLayer(rapMean, {min:0, max: 5, palette: 'white,green'}, 'tree cov', false)

Map.addLayer(rapMean.lte(ee.Image(1)).selfMask(), {min:0, max: 1, palette: 'black'}, '<=1% trees', false)
Map.addLayer(rapMean.lt(ee.Image(1)).selfMask(), {min:0, max: 1, palette: 'black'}, '<1% trees', false)
Map.addLayer(rapMean.lte(ee.Image(0.5)).selfMask(), {min:0, max: 1, palette: 'black'}, '1/2% trees', false)
Map.addLayer(rapMean.eq(ee.Image(0)).selfMask(), {min:0, max: 1, palette: 'black'}, '0% trees', false)


var notForest = rapMean
  .lt(ee.Image(1)) // less than 1% tree cover
  .selfMask()
  .toByte()
  .rename('TRE_lt1');

// export -------------------------------------------
var fileName = 'RAP_v3_tree-lt1_' + yearStart + '-' + yearEnd  + '_30m';


Export.image.toAsset({
  image: notForest.toByte(),
  description: fileName,
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + fileName,
  crs: fg.crs,
  scale: 30,
  region: fg.region,
  maxPixels: 1e12,
  pyramidingPolicy: {'TRE_lt1': 'mode'}
});

