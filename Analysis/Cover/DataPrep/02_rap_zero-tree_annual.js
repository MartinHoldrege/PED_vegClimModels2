/*
Annual zero-tree layers on the daymet 1 km snap grid, yearStart-yearEnd: the
share of a cell's kept 30 m pixels with RAP tree cover < cutoff
(fracZeroTree), and zeroTree = 1 where that share is > zeroTreeFracCutoff.
Same definition as the 2021 layer (02_rap_no-tree.js, 03_rap_frac-zero-tree.js),
for every year. A 30 m pixel is kept if RAP has data there, it is unburned in
the year's 20-year window, and the LCMAP mask keeps it. The burned layer is
the 100 m asset from 01_fire_everBurned_windows.js, so burn edges are at
100 m rather than 30 m.

Built the same way as 01_rap_cover.js, which ran affordably for all years
(the earlier per-cell-year version, 07_rap_zero-tree_cell-years.js, was
far too expensive). Three steps, run in order (set `step`):
  1. For each year, two 30 m 0/1 bands: keep, and notTreeKeep (kept and tree
     cover < cutoff). Averaged to 250 m (reduceResolution) on a grid with the
     snap grid's origin, and saved as a float asset with mean pyramids. 250 m
     nests exactly in the 1 km cells, so each asset's 1000 m pyramid level is
     the mean of each band over the 1 km cell.
  2. (after all step 1 assets have finished) read the assets at 1 km and
     export fracZeroTree = notTreeKeep / keep (the share of kept pixels with
     no trees; NA where nothing is kept) and zeroTree. Two GeoTIFFs, each
     with one band per year ('year_YYYY'), as the RAP cover files.
  3. Delete the 250 m assets (once the step 2 files are downloaded and
     checked).
Cells outside CONUS or failing the 1 km LCMAP mask or that year's 1 km fire
mask are left NA, since the training data (06) never use them; this also
cuts compute. FIA-only field rows in LCMAP-masked cells (not LCMAP-masked in
06) therefore get no zero-tree value; they're only used for the c_needle /
c_broad validation.

To split the run across months (GEE quota), run steps 1-2 for part of the
years at a time (yearStart, yearEnd).

Inputs:
  - RAP v3 cover
  - <fg.pathAsset>fire/MTBS_burned20yr_100m_2000-2024
    (01_fire_everBurned_windows.js)
  - fg.lcmapMask (Functions/gee/general.js)
Outputs (Drive), downloaded with 04_download_gee_output.R:
  - RAP_v3_annualFracZeroTree_lt3_<yearStart>-<yearEnd>_1000m.tif
  - RAP_v3_annualZeroTree_lt3_<yearStart>-<yearEnd>_1000m.tif

Replaces 07_rap_zero-tree_cell-years.js.

Author: Martin Holdrege
Started: November 2026
*/

// params -------------------------------------------
var cutoff = 3;               // % tree cover, as in 02_rap_no-tree.js
var zeroTreeFracCutoff = 0.9; // as in 03_rap_frac-zero-tree.js
var yearStart = 2010;
var yearEnd = 2024;
var step = 2;  // 1: make the 250 m assets; 2: export at 1 km to Drive;
               // 3: delete the 250 m assets

var driveFolder = 'PED_vegClimModels2';

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// 250 m grid with the same origin as the 1 km snap grid
var crsTransform250 = [250, 0, fg.crsTransform[2], 0, -250, fg.crsTransform[5]];

var assetId = function(year) {
  return fg.pathAsset + 'rap/RAP_v3_zeroTreeInputs_lt' + cutoff + '_' + year +
    '_250m';
};

// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');

// 100 m, band 'year_<Y>': 1 = burned in Y-19 to Y (unmapped counted as
// unburned)
var burnedWindow = ee.Image(fg.pathAsset + 'fire/MTBS_burned20yr_100m_2000-2024');

// 30 m bands for one year -------------------------
// keep: 1 = RAP has data, unburned in the year's window, and kept by LCMAP.
// notTreeKeep: 1 = kept and tree cover < cutoff. Within a cell both are 0
// (not masked) where not kept, so their 250 m and 1 km means are shares of
// all the cell's pixels, and notTreeKeep / keep is exactly the share of kept
// pixels with no trees.
// Whole 1 km cells are then masked where the training data never use them
// (outside CONUS, failing the 1 km LCMAP mask or that year's 1 km fire mask,
// as in 06), after the unmask(0), so GEE has less to compute and store.
// setDefaultProjection: reduceResolution averages the input at its default
// projection, so this makes the input RAP's 30 m grid (the other layers are
// sampled onto it) rather than a coarser one
var bands30 = function(year) {
  var tree = ee.Image(rap.filter(ee.Filter.calendarRange(year, year, 'year'))
    .first()).select('TRE');
  var unburned = burnedWindow.select('year_' + year).eq(0);

  var keep = tree.mask().and(unburned).and(fg.lcmapMask).unmask(0);
  var notTreeKeep = keep.and(tree.lt(cutoff).unmask(0));

  var cellKeep = fg.maskConus.gt(0)
    .and(fg.lcmapMaskBinary())
    .and(fg.fireMaskYear(year));

  return ee.Image.cat([keep, notTreeKeep])
    .rename(['keep', 'notTreeKeep'])
    .toByte()
    .updateMask(cellKeep)
    .setDefaultProjection(tree.projection());
};

// step 1: 250 m assets, one per year ---------------
if (step === 1) {
  for (var year = yearStart; year <= yearEnd; year++) {
    // maxPixels: ~70-100 RAP pixels per 250 m cell (more at high latitudes),
    // as in 01_rap_cover.js. If it's too low, GEE silently reads a coarser
    // level instead, so leave headroom
    var img250 = bands30(year)
      .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 512})
      .reproject({crs: fg.crs, crsTransform: crsTransform250})
      .toFloat();  // float, so the pyramid means aren't rounded

    Export.image.toAsset({
      image: img250,
      description: 'RAP_v3_zeroTreeInputs_lt' + cutoff + '_' + year + '_250m',
      assetId: assetId(year),
      crs: fg.crs,
      crsTransform: crsTransform250,
      region: fg.region,
      maxPixels: 1e13,
      pyramidingPolicy: {'.default': 'mean'}
    });
  }
}

// step 2: 1 km exports -----------------------------
// Read at 1 km, each band is its mean over the cell, so the division gives
// the share of the cell's kept pixels with no trees
if (step === 2) {
  var conus = fg.maskConus.gt(0);

  var fracBands = [];
  var zeroBands = [];
  for (var y = yearStart; y <= yearEnd; y++) {
    var a = ee.Image(assetId(y));
    var keepFrac = a.select('keep');
    // 1 km masks applied here too (as in 03), so masked cells are NA even if
    // the step 1 assets were made without them
    var cellKeep = fg.lcmapMaskBinary().and(fg.fireMaskYear(y));
    var frac = a.select('notTreeKeep').divide(keepFrac)
      .updateMask(keepFrac.gt(0))
      .updateMask(cellKeep)
      .rename('year_' + y);
    fracBands.push(frac);
    zeroBands.push(frac.gt(zeroTreeFracCutoff).rename('year_' + y));
  }

  var fracOut = ee.Image.cat(fracBands).updateMask(conus).toFloat();
  var zeroOut = ee.Image.cat(zeroBands).updateMask(conus).toByte();

  // 'annual' keeps these apart from the 2021 layer from 03 (same years in a
  // 2021-only test would otherwise give the same file name)
  var label = '_lt' + cutoff + '_' + yearStart + '-' + yearEnd + fg.resLabel;
  fg.exportDrive(fracOut, 'RAP_v3_annualFracZeroTree' + label, driveFolder);
  fg.exportDrive(zeroOut, 'RAP_v3_annualZeroTree' + label, driveFolder);

  Map.addLayer(zeroOut.select('year_' + yearEnd).selfMask(),
    {palette: ['black']}, 'zero tree ' + yearEnd, false);
}

// step 3: delete the 250 m assets ------------------
// permanent; only run once the step 2 files are downloaded and checked
if (step === 3) {
  for (var d = yearStart; d <= yearEnd; d++) {
    // skip assets that are already gone, but print the error so other
    // failures (e.g. permissions) are still visible
    try {
      ee.data.deleteAsset(assetId(d));
      print('deleted', assetId(d));
    } catch (e) {
      print('not deleted', assetId(d), e.message);
    }
  }
}
