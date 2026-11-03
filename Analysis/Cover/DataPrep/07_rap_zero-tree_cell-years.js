/*
Zero-tree status at the cell-years of the cover training data (c04-c06), each
for its own year. Used as the zero-tree model response and for the
c_needle/c_broad zeroing in 08_combine_cover_and_covariates.R.

The 2021 layer's definition (02_rap_no-tree.js and 03_rap_frac-zero-tree.js),
applied to each year, but computed at `scale` m rather than 30 m:
  - `scale` m: 1 where RAP tree cover < cutoff, masked to pixels unburned in
    that year's 20-year window and kept by the LCMAP mask. Above 30 m, GEE
    reads RAP's averaged pyramid levels, so the threshold applies to mean
    tree cover over that larger pixel. The burned layer is the 100 m asset
    from 01_fire_everBurned_windows.js (one band per year)
  - 1 km: mean of that over the snap cell (fracZeroTree), over each cell's
    square; zeroTree = 1 where fracZeroTree > zeroTreeFracCutoff. Unlike 03,
    no 1 km masks: the cell-years were already filtered by them in 06
    (except FIA-only rows, which aren't LCMAP-masked and so also get a value
    here)
Only the uploaded cells are computed (reduceRegions over each cell's square),
so there is no CONUS-wide step.

Input:
  - cover_cell_years_c04-c06.csv (06_build_cover_datasets.R), uploaded as a
    table asset to <fg.pathAsset>cover/cover_cell_years_c04-c06, with
    longitude and latitude as the X and Y columns. The x, y columns (cell
    centres in the daymet CRS) are used to build each cell's square.
  - <fg.pathAsset>fire/MTBS_burned20yr_100m_2000-2024
    (01_fire_everBurned_windows.js)

Output (Drive): RAP_v3_zeroTree_lt3_cell-years_c04-c06.csv, downloaded with
04_download_gee_output.R

November 3, 2026: changed to cut compute (the earlier version used up the
monthly GEE quota): the burned layer now comes from the 100 m asset instead
of being rebuilt from MTBS at 30 m, and the computation is at `scale` m
rather than 30 m. The earlier c04-c06 output was made with the previous
version of this file (see the git history). Values from this version can
differ slightly from those, because of the coarser resolution.

Author: Martin Holdrege
Started: October 2026
*/

// params -------------------------------------------
var cutoff = 3;               // % tree cover, as in 02_rap_no-tree.js
var zeroTreeFracCutoff = 0.9; // as in 03_rap_frac-zero-tree.js
var scale = 60;               // m; 30 is the 2021 layer's resolution
var tileScale = 4;
var yearStart = 2000;
var yearEnd = 2024;
var vcLabel = 'c04-c06';
var driveFolder = 'PED_vegClimModels2';

// > 0: only this many randomly chosen cell-years, exported to a _test file,
// to check compute (EECU-seconds in the task details) before the full run
var nTest = 2000;

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -------------------------------------
var cellYears = ee.FeatureCollection(fg.pathAsset + 'cover/cover_cell_years_' +
  vcLabel);

// checks the upload and that the year filter below works on its year column
print('cell-years in the upload', cellYears.size());
print('cell-years in 2021', cellYears.filter(ee.Filter.eq('year', 2021)).size());

if (nTest > 0) {
  cellYears = cellYears.randomColumn('r', 1).sort('r').limit(nTest);
}

var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');

// 100 m, band 'year_<Y>': 1 = burned in Y-19 to Y (unmapped counted as
// unburned); created in 01_fire_everBurned_windows.js
var burnedWindow = ee.Image(fg.pathAsset + 'fire/MTBS_burned20yr_100m_2000-2024');

// process ------------------------------------------

// each cell-year's geometry becomes its 1 km snap cell (x, y are the centre)
var proj = ee.Projection(fg.crs);
var cells = cellYears.map(function(f) {
  var x = ee.Number(f.get('x'));
  var y = ee.Number(f.get('y'));
  var square = ee.Geometry.Rectangle(
    [x.subtract(500), y.subtract(500), x.add(500), y.add(500)], proj, false);
  return f.setGeometry(square);
});

// layer for one year: 1 = tree cover < cutoff, masked to natural land
// unburned in the year's window
var notTreeYear = function(year) {
  var tree = ee.Image(rap.filter(ee.Filter.calendarRange(year, year, 'year'))
    .first()).select('TRE');
  var unburned = burnedWindow.select('year_' + year).eq(0);

  return tree.lt(cutoff)
    .updateMask(unburned.and(fg.lcmapMask));
};

// zeroTree from fracZeroTree; left empty where fracZeroTree is (cell with no
// unmasked 30 m pixels)
var addZeroTree = function(f) {
  var frac = f.get('fracZeroTree');
  return f.set('zeroTree', ee.Algorithms.If(
    ee.Algorithms.IsEqual(frac, null),
    null,
    ee.Number(frac).gt(zeroTreeFracCutoff)));
};

// mean over each cell, by year. The mean weights pixels by how much of each
// lies in the cell, as reduceResolution does in 03
var sampled = ee.FeatureCollection([]);
for (var year = yearStart; year <= yearEnd; year++) {
  var s = notTreeYear(year).reduceRegions({
    collection: cells.filter(ee.Filter.eq('year', year)),
    reducer: ee.Reducer.mean().setOutputs(['fracZeroTree']),
    crs: fg.crs,
    scale: scale,
    tileScale: tileScale
  });
  sampled = sampled.merge(s.map(addZeroTree));
}

// export -------------------------------------------
var fileName = 'RAP_v3_zeroTree_lt' + cutoff + '_cell-years_' + vcLabel +
  (nTest > 0 ? '_test' : '');

Export.table.toDrive({
  collection: sampled,
  description: fileName,
  folder: driveFolder,
  fileFormat: 'CSV',
  selectors: ['cell', 'year', 'fracZeroTree', 'zeroTree']
});

if(false) {
  // check: computed 2021 layer vs the stored 02 asset, at 300 2021 cells
  var test = cells.filter(ee.Filter.eq('year', 2021)).limit(300);
  var stored = ee.Image(fg.pathAsset + 'rap/RAP_v3_tree-lt3_masked_2021-2021_30m');
  var both = notTreeYear(2021).rename('computed')
    .addBands(stored.rename('stored'));
  var res = both.reduceRegions({collection: test, reducer: ee.Reducer.mean(),
                                crs: fg.crs, scale: 30});
  print(ui.Chart.feature.byFeature(res, 'stored', ['computed'])
    .setChartType('ScatterChart'));
}