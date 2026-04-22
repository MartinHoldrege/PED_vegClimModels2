/*
Get RAP v3 cover on the daymet grid, for the years of the biomass datasets.
Unmasked asset for flexible downstream use.

Author: Martin Holdrege
Started: April 6, 2026
*/

// params -------------------------------------------
// var yearStart = 2019; // years corresponding to GEDI dataset
// var yearEnd = 2023;

var exportCov = false;
var exportBio = true; 

var yearStart = 2000; // years corresponding to cover model training data
var yearEnd = 2023;

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');



// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));

// process ------------------------------------------
var rapMean = rap
  .select(['TRE', 'SHR', 'AFG', 'PFG'])
  .mean();

var coverStack = rapMean.select('TRE').rename('totalTreeCov')
  .addBands(rapMean.select('SHR').rename('totalShrubCov'))
  .addBands(rapMean.select('AFG').add(rapMean.select('PFG')).rename('totalHerbaceousCov'));
  
// herbaceous biomass (Mg/ha)
var herbBiomass = fg.rapHerbBiomass(yearStart, yearEnd);

// export -------------------------------------------
var fileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd  + fg.resLabel;

if(exportCov) {
  Export.image.toAsset({
    image: coverStack,
    description: fileName,
    assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + fileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12
  });
}

if(exportBio) {
  var biomassFileName = 'RAP_v3_herbaceousAGB_' + yearStart + '-' + yearEnd + fg.resLabel;
  Export.image.toAsset({
    image: herbBiomass,
    description: biomassFileName,
    assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + biomassFileName,
    crs: fg.crs,
    crsTransform: fg.crsTransform,
    region: fg.region,
    maxPixels: 1e12
  });
}

