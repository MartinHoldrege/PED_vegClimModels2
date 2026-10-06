/*
Determine which pixels have no or close to no tree cover (based on RAP),
masked to exclude burned (MTBS, 20 years up to and including yearEndRap) and
developed/cropland/water (LCMAP 2021) pixels.

Author: Martin Holdrege
Started: April 6, 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStartRap = 2021;
var yearEndRap = 2021;
var windowLength = 20; // fire window (years), inclusive of yearEndRap

var cutoffs = [3, 5, 10];
// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStartRap, yearEndRap, 'year'));

var mtbs = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1')
  .filter(ee.Filter.stringContains('system:index', 'CONUS'))
  .filter(ee.Filter.calendarRange(yearEndRap - windowLength + 1, yearEndRap, 'year'))
  .map(function(img) {
    // band is named 'Burn_Severity' in some years, so select by position
    var severity = img.select([0]);
    return severity.gte(2).and(severity.lte(5));
  });

// burned in window; unmask so unmapped areas count as unburned
var everBurned = mtbs.max().unmask(0);

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
