/*
Compute fraction of each daymet cell unburned (MTBS 2000-2023).
Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStart = 2000;
var yearEnd = 2023;

// read in data -------------------------------------
var mtbs = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1')
  .filter(ee.Filter.stringContains('system:index', 'CONUS'))
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .map(function(img) {
    // in some cases band is 'Burn_Severity' 
    return img.select([0]).rename('Severity');
  });

// process ------------------------------------------
var bandNames = mtbs.map(function(img) {
  return img.set('bands', img.bandNames());
});

// for each year: 1 = burned (classes 2-5), 0 = unburned (or low) (0, 1, 6)
var burned = mtbs.map(function(img) {
  var severity = img.select('Severity');
  return severity.gte(2).and(severity.lte(5))
    .copyProperties(img, ['system:time_start']);
});

// sum across years, then collapse to binary (ever burned = 1)
var everBurned = burned.sum().gte(1);

// fraction unburned = 1 - everBurned
var unburned = everBurned.not()
  .setDefaultProjection(mtbs.first().select('Severity').projection());

// aggregate to daymet grid
var fracUnburned = unburned
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

// visualize ----------------------------------------
// Map.addLayer(fracUnburned, {min: 0, max: 1, palette: ['red', 'white', 'green']}, 'fraction unburned');
// Map.addLayer(fracUnburned.gte(0.9).selfMask(), {palette: ['blue']}, '>=90% unburned', false);

// export -------------------------------------------
var fileName = 'MTBS_fracUnburned_' + yearStart + '-' + yearEnd + fg.resLabel;

Export.image.toAsset({
  image: fracUnburned,
  description: fileName,
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});