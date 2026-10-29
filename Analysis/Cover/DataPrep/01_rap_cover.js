/*
Annual RAP v3 cover (tree, shrub, herbaceous, bare ground), 2000-2024, as
true means over the cells of the daymet 1 km snap grid.

Two steps, run in order (set `step`):
  1. For each year, RAP 30 m -> 250 m mean (reduceResolution) on a grid with
     the snap grid's origin, saved as a float asset with mean pyramids. 250 m
     nests exactly in the 1 km cells, so each asset's 1000 m pyramid level is
     the mean of RAP over each 1 km cell.
  2. (after all step 1 assets have finished) export each group at 1 km on the
     snap grid (fg.exportDrive), which reads that pyramid level. One GeoTIFF
     per group, one band per year ('year_YYYY').

Not masked by LCMAP or fire. The masks are applied in R: the RAP sample uses
the number of unmasked years per cell, and tree cover is also needed at field
plots in cells the LCMAP mask drops (FIA).

Starts in 2000 so tree cover is available for every year of field data.

The 250 m assets are only intermediates; delete them once the step 2 exports
are in Drive.

Replaces 03_rap_sample.js and 03_rap_cover_2021_250m_test.js.

Author: Martin Holdrege
Started: October 2026
*/

// params -------------------------------------------
var yearStart = 2000;
var yearEnd = 2024;
var step = 1;  // 1: make the 250 m assets; 2: export them at 1 km to Drive;
               // 3: delete the 250 m assets

var driveFolder = 'PED_vegClimModels2';

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// 250 m grid with the same origin as the 1 km snap grid
var crsTransform250 = [250, 0, fg.crsTransform[2], 0, -250, fg.crsTransform[5]];

var pfts = ['tree', 'shrub', 'herbaceous', 'bare_ground'];

var assetId = function(year) {
  return fg.pathAsset + 'rap/RAP_v3_cover_' + year + '_250m';
};

// functional groups; herbaceous = annual + perennial forbs and grasses
var pftBands = function(image) {
  return ee.Image.cat([
    image.select('TRE'),
    image.select('SHR'),
    image.select('AFG').add(image.select('PFG')),
    image.select('BGR')
  ]).rename(pfts);
};

// read in data -------------------------------------
var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');

// check yearEnd is in the collection before running step 1
print('RAP years available', rap.aggregate_array('system:time_start')
  .map(function(t) { return ee.Date(t).get('year'); }));

// step 1: 250 m assets, one per year ---------------
if (step === 1) {
  for (var year = yearStart; year <= yearEnd; year++) {
    var rapYear = ee.Image(
      rap.filter(ee.Filter.calendarRange(year, year, 'year')).first()
    );

    // maxPixels: ~70-100 RAP pixels per 250 m cell (more at high latitudes,
    // as RAP's pixels narrow in metres). If it's too low, GEE silently reads
    // a coarser RAP pyramid level instead, so leave headroom
    var rap250 = pftBands(rapYear)
      .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 512})
      .reproject({crs: fg.crs, crsTransform: crsTransform250})
      .toFloat();  // float, so the pyramid means aren't rounded

    Export.image.toAsset({
      image: rap250,
      description: 'RAP_v3_cover_' + year + '_250m',
      assetId: assetId(year),
      crs: fg.crs,
      crsTransform: crsTransform250,
      region: fg.region,
      maxPixels: 1e13,
      pyramidingPolicy: {'.default': 'mean'}
    });
  }
}

// step 2: 1 km exports, one per group --------------
if (step === 2) {
  var conus = fg.maskConus.gt(0);

  for (var k = 0; k < pfts.length; k++) {
    var pft = pfts[k];

    var bands = [];
    for (var y = yearStart; y <= yearEnd; y++) {
      bands.push(ee.Image(assetId(y)).select(pft).rename('year_' + y));
    }

    var out = ee.Image.cat(bands).updateMask(conus).toFloat();

    var fileName = 'RAP_v3_cover-' + pft + '_' + yearStart + '-' + yearEnd +
      fg.resLabel;
    fg.exportDrive(out, fileName, driveFolder);

    if (pft === 'tree') {
      Map.addLayer(out.select('year_' + yearEnd),
        {min: 0, max: 60, palette: ['white', 'darkgreen']},
        'tree ' + yearEnd + ', 1 km', false);
    }
  }
}

// step 3: delete the 250 m assets ------------------
// permanent; only run once the step 2 files are downloaded and checked
if (step === 3) {
  for (var d = yearStart; d <= yearEnd; d++) {
    ee.data.deleteAsset(assetId(d));
    print('deleted', assetId(d));
  }
}