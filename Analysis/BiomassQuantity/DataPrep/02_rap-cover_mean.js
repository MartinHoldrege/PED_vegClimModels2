/*
Get RAP v3 cover on the daymet grid, for the years of the biomass datasets.
Unmasked asset for flexible downstream use.

Author: Martin Holdrege
Started: April 6, 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
// var yearStart = 2019; // years corresponding to GEDI dataset
// var yearEnd = 2023;

var yearStart = 2010; // years corresponding to spawn dataset
var yearEnd = 2010;

// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rangeland-analysis-platform/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));

// process ------------------------------------------
var rapMean = rap
  .select(['TRE', 'SHR', 'AFG', 'PFG'])
  .mean()
  .setDefaultProjection(rap.first().projection());

var coverStack = rapMean.select('TRE').rename('totalTreeCov')
  .addBands(rapMean.select('SHR').rename('totalShrubCov'))
  .addBands(rapMean.select('AFG').add(rapMean.select('PFG')).rename('totalHerbaceousCov'));

var coverDaymet = coverStack
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 2e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform
  });

// export -------------------------------------------

var fileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd + fg.resLabel;

Export.image.toAsset({
  image: coverDaymet,
  description: fileName,
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});