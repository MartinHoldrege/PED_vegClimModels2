/*
For each year 2000-2024, whether a pixel burned in the preceding 20 years
(year Y window = Y-19 through Y, inclusive). MTBS starts in 1984, so the
earliest target years have a shorter window.

Exported on the daymet CRS at 100m. The base level is binary (burned/unburned),
with mean pyramiding so that coarser levels are the burned proportion. This is
what lets the 1km proportion be computed in a later step without a CONUS-scale
reduceResolution from 30m.

Author: Martin Holdrege
Started: August 2026
*/

// params -------------------------------------------
var yearStart = 2000;
var yearEnd = 2024;
var windowLength = 20;    // years, inclusive of the target year
var mtbsYearStart = 1984; // first year of the MTBS record
var scaleOut = 100;       // m

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -------------------------------------
var mtbs = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1')
  .filter(ee.Filter.stringContains('system:index', 'CONUS'))
  .map(function(img) {
    // band is named 'Burn_Severity' in some years, so select by position
    return img.select([0]).rename('Severity');
  });

print('MTBS years available', mtbs.aggregate_array('system:index'));

// process ------------------------------------------
// 1 = burned (severity classes 2-5), 0 = unburned, low, or non-mapping (0, 1, 6)
var burned = mtbs.map(function(img) {
  var severity = img.select('Severity');
  return severity.gte(2).and(severity.lte(5))
    .copyProperties(img, ['system:time_start']);
});

// one band per target year: burned any time in that year's window
var stack = null;
var nYears = [];
var yearsUsed = [];

for (var year = yearStart; year <= yearEnd; year++) {
  // clamp to the start of the MTBS record so the counts below are honest
  var yearFirst = Math.max(year - windowLength + 1, mtbsYearStart);
  var window = burned.filter(ee.Filter.calendarRange(yearFirst, year, 'year'));

  // unmask so unmapped areas count as unburned rather than staying masked
  var band = ee.Image(window.max()).unmask(0).rename('year_' + year);
  stack = (stack === null) ? band : stack.addBands(band);

  nYears.push(year - yearFirst + 1);
  yearsUsed.push(yearFirst + '-' + year);
}

// band metadata. GEE has no per-band properties, so these lists are
// positionally parallel to bandNames(): band i covers yearsUsed[i].
stack = stack.toFloat()
  .set({
    windowLength: windowLength,
    scale_m: scaleOut,
    bandYearRange: yearsUsed,
    bandNYears: nYears,
    severityClasses: '2-5 = burned',
    mtbsYearStart: mtbsYearStart
  });

print('band names', stack.bandNames());
print('years per band', nYears);

// visualize ----------------------------------------
Map.addLayer(stack.select(String(yearEnd)).selfMask(),
  {min: 0, max: 1, palette: ['white', 'black']},
  'burned ' + yearsUsed[yearsUsed.length - 1], false);

// export -------------------------------------------
var fileName = 'MTBS_burned' + windowLength + 'yr_' + scaleOut + 'm_' +
  yearStart + '-' + yearEnd;

Export.image.toAsset({
  image: stack,
  description: fileName,
  assetId: fg.pathAsset + 'fire/' + fileName,
  crs: fg.crs,
  scale: scaleOut,
  region: fg.region,
  maxPixels: 1e12,
  // coarser pyramid levels are then the burned proportion
  pyramidingPolicy: {'.default': 'mean'}
});