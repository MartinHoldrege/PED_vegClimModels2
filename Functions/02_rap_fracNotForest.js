/*
Fraction of unburned, non-developed 30m pixels per daymet cell that have <1% tree cover (RAP).
Burned pixels (MTBS 2000-2023) and developed/cropland/water pixels (LCMAP 2021)
are excluded from both numerator and denominator.
Used for getting total treeless area for calibrating forest/not forest model

Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStartRap = 2019;
var yearEndRap = 2023;
var yearStartFire = 2000;
var yearEndFire = 2023;

// read in data -------------------------------------

// created in 03_rap_notForest.js
var notForest30 = ee.Image(fg.pathAsset + 'rap/RAP_v3_tree-lt1_' +
  yearStartRap + '-' + yearEndRap + '_30m')
  .unmask(0);

// created in 01_mtbs_everBurned.js
var everBurned = ee.Image(fg.pathAsset + 'fire/MTBS_everBurned_30m_' +
  yearStartFire + '-' + yearEndFire)
  .unmask(0);

// LCMAP 2021 land cover
var lcmap2021 = ee.ImageCollection('projects/sat-io/open-datasets/LCMAP/LCPRI')
  .filterDate('2021-01-01', '2021-12-31')
  .first();

// process ------------------------------------------

// LCMAP mask: 1 = keep (grass/shrub, tree, wetlands, ice/snow, barren)
// 0 = remove (developed, cropland, water)
var lcmapKeep = lcmap2021.remap(
  [1, 2, 3, 4, 5, 6, 7, 8],
  [0, 0, 1, 1, 0, 1, 1, 1]
);

// combined mask: unburned AND not developed/cropland/water
var keepMask = everBurned.not().and(lcmapKeep);

var notForestMasked = notForest30.updateMask(keepMask);

// this layer looks a little weird when zoomed out b/ of the modal
// pyramiding
Map.addLayer(notForestMasked.selfMask(), {}, 'not forest', false)

var fracNotForest = notForestMasked
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 2e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform
  })
  .rename('fracNotForest');

// visualize ----------------------------------------
Map.addLayer(fracNotForest, {min: 0, max: 1, palette: ['white', 'darkgreen']},
  'frac not forest (unburned, natural land only)', false);

// export -------------------------------------------
var fileName = 'RAP_v3_fracNotForest_' +
  yearStartRap + '-' + yearEndRap + fg.resLabel;

Export.image.toAsset({
  image: fracNotForest,
  description: fileName,
  assetId: fg.pathAsset + 'rap/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});