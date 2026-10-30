/*
Zero-tree status at the cell-years of the cover training data (c04-c06), each
for its own year. Used as the zero-tree model response and for the
c_needle/c_broad zeroing in 08_combine_cover_and_covariates.R.

Same definition as the 2021 layer (02_rap_no-tree.js and
03_rap_frac-zero-tree.js), applied to each year:
  - 30 m: 1 where RAP tree cover < cutoff, masked to pixels unburned in that
    year's 20-year window and kept by the LCMAP mask; resampled onto a 30 m
    grid in the daymet CRS, as the 02 export is
  - 1 km: mean of that over the snap cell (fracZeroTree), masked by the 1 km
    LCMAP mask and that year's fire mask; zeroTree = 1 where
    fracZeroTree > zeroTreeFracCutoff
Values are sampled at the cell centres, so GEE only computes the sampled
cells (no CONUS-wide 30 m export).

Input: cover_cell_years_c04-c06.csv (06_build_cover_datasets.R), uploaded as
a table asset to <fg.pathAsset>cover/cover_cell_years_c04-c06, with
longitude and latitude as the X and Y columns.

Output (Drive): RAP_v3_zeroTree_lt3_cell-years_c04-c06.csv, downloaded with
04_download_gee_output.R

Author: Martin Holdrege
Started: October 2026
*/

// params -------------------------------------------
var cutoff = 3;               // % tree cover, as in 02_rap_no-tree.js
var zeroTreeFracCutoff = 0.9; // as in 03_rap_frac-zero-tree.js
var windowLength = 20;        // fire window (years), inclusive of the year
var yearStart = 2000;
var yearEnd = 2024;
var vcLabel = 'c04-c06';
var driveFolder = 'PED_vegClimModels2';

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -------------------------------------
var cellYears = ee.FeatureCollection(fg.pathAsset + 'cover/cover_cell_years_' +
  vcLabel);

// checks the upload and that the year filter below works on its year column
print('cell-years in the upload', cellYears.size());
print('cell-years in 2021', cellYears.filter(ee.Filter.eq('year', 2021)).size());

var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');

// 1 = burned (severity classes 2-5), as in 02_rap_no-tree.js
var burned = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1')
  .filter(ee.Filter.stringContains('system:index', 'CONUS'))
  .map(function(img) {
    // band is named 'Burn_Severity' in some years, so select by position
    var severity = img.select([0]);
    return severity.gte(2).and(severity.lte(5))
      .copyProperties(img, ['system:time_start']);
  });

// process ------------------------------------------

// fracZeroTree and zeroTree on the snap grid, for one year
var zeroTreeYear = function(year) {
  var tree = ee.Image(rap.filter(ee.Filter.calendarRange(year, year, 'year'))
    .first()).select('TRE');

  // burned in the window; unmask so unmapped areas count as unburned
  var everBurned = burned
    .filter(ee.Filter.calendarRange(year - windowLength + 1, year, 'year'))
    .max()
    .unmask(0);

  var notTree30 = tree.lt(cutoff)
    .updateMask(everBurned.not().and(fg.lcmapMask))
    .reproject({crs: fg.crs, scale: 30});

  var cellKeep = fg.lcmapMaskBinary().and(fg.fireMaskYear(year));

  var frac = notTree30
    .reduceResolution({reducer: ee.Reducer.mean(), bestEffort: true,
                       maxPixels: 2e3})
    .reproject({crs: fg.crs, crsTransform: fg.crsTransform})
    .updateMask(cellKeep)
    .rename('fracZeroTree');

  return frac.addBands(frac.gt(zeroTreeFracCutoff).rename('zeroTree'));
};

// sample each year's cell-years from that year's layer; masked cells give an
// empty value
var sampled = ee.FeatureCollection([]);
for (var year = yearStart; year <= yearEnd; year++) {
  var pts = cellYears.filter(ee.Filter.eq('year', year));
  var s = zeroTreeYear(year).reduceRegions({
    collection: pts,
    reducer: ee.Reducer.first(),
    crs: fg.crs,
    crsTransform: fg.crsTransform
  });
  sampled = sampled.merge(s);
}

// visualize ----------------------------------------
Map.addLayer(zeroTreeYear(2021).select('zeroTree').selfMask(),
  {palette: ['black']}, 'zero tree 2021 (should match 03 layer)', false);

// export -------------------------------------------
var fileName = 'RAP_v3_zeroTree_lt' + cutoff + '_cell-years_' + vcLabel;

Export.table.toDrive({
  collection: sampled,
  description: fileName,
  folder: driveFolder,
  fileFormat: 'CSV',
  selectors: ['cell', 'year', 'fracZeroTree', 'zeroTree']
});
