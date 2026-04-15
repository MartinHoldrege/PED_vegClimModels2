/*
determeine whether a grid cell burned in MTBS 2000-2023, at daymet resolution.
this is an intermediate needed for not getting reprojection issues
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

var everBurned = burned.max(); // 1 if burned in any year

var mtbsProj = mtbs.first().projection();

Map.addLayer(everBurned.selfMask(), {min: 0, max: 1, palette: ['white', 'black']}, 'ever burned', false);

Export.image.toAsset({
  image: everBurned.rename('everBurned'),
  description: 'MTBS_everBurned_30m_' + yearStart + '-' + yearEnd,
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/fire/MTBS_everBurned_30m_' + yearStart + '-' + yearEnd,
  crs: mtbsProj.crs().getInfo(),
  scale: 30,
  region: fg.region,
  maxPixels: 1e12
});


