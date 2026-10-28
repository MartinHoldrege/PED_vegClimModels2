/*
Test: true 1 km means of RAP v3 cover (2021) on the snap grid, to compare with
RAP_v3_cover_2021_1000m from 03_rap_sample.js (which reads RAP's own integer
pyramids, so its 1 km values are rounded and come from ~960 m blocks that
don't line up with the snap grid).

Two steps, run in order (set `step`):
  1. RAP 30 m -> 250 m mean (reduceResolution) on a grid with the snap grid's
     origin, saved as a float asset with mean pyramids. 250 m nests exactly in
     the 1 km cells, so the asset's 1000 m pyramid level is the mean of the 16
     250 m cells = the mean of RAP over each 1 km cell.
  2. (after the step 1 asset has finished) export the asset at 1 km on the snap
     grid (fg.exportDrive), which reads that pyramid level. Same masks as the
     2021 wall-to-wall export in 03_rap_sample.js.

Author: Martin Holdrege
Started: October 2026
*/

// params -------------------------------------------
var year = 2021;
var step = 1;  // 1: make the 250 m asset; 2: export it at 1 km to Drive

var driveFolder = 'PED_vegClimModels2';

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// 250 m grid with the same origin as the 1 km snap grid
var crsTransform250 = [250, 0, fg.crsTransform[2], 0, -250, fg.crsTransform[5]];

var assetName = 'RAP_v3_cover_' + year + '_250m';
var assetId = fg.pathAsset + assetName;

// functional groups, as in 03_rap_sample.js
var pftBands = function(image) {
  return ee.Image.cat([
    image.select('TRE'),
    image.select('SHR'),
    image.select('AFG').add(image.select('PFG')),
    image.select('BGR')
  ]).rename(['tree', 'shrub', 'herbaceous', 'bare_ground']);
};

// step 1: 250 m asset ------------------------------
if (step === 1) {
  var rap = ee.Image(
    ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .first()
  );

  // maxPixels: ~70-100 RAP pixels per 250 m cell (more at high latitudes, as
  // RAP's pixels narrow in metres). If it's too low, GEE silently reads a
  // coarser RAP pyramid level instead, so leave headroom
  var rap250 = pftBands(rap)
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 512})
    .reproject({crs: fg.crs, crsTransform: crsTransform250})
    .toFloat();  // float, so the pyramid means aren't rounded

  Export.image.toAsset({
    image: rap250,
    description: assetName,
    assetId: assetId,
    crs: fg.crs,
    crsTransform: crsTransform250,
    region: fg.region,
    maxPixels: 1e13,
    pyramidingPolicy: {'.default': 'mean'}
  });
}

// step 2: 1 km export from the asset's pyramid -----
if (step === 2) {
  var lcmapMask = fg.lcmapMaskBinary();  // 90% keep threshold
  var fullKeep = fg.maskConus.gt(0).and(lcmapMask).and(fg.fireMaskYear(year));

  var out = ee.Image(assetId).updateMask(fullKeep);

  fg.exportDrive(out, 'RAP_v3_cover_' + year + '_from250m' + fg.resLabel,
                 driveFolder);

  Map.addLayer(out.select('tree'), {min: 0, max: 60, palette: ['white', 'darkgreen']},
               'tree, 1 km from 250 m asset', false);
}