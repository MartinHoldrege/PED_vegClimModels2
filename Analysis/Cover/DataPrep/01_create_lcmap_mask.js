/*
Purpose: calculate % of each daymet gridcell that is water, ag or developed, so 
this can be used as a mask

Author: Martin Holdrege
Script started: April 3, 2026
*/

// dependencies -------------------------------------

// functions general
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -----------------------------
var lcpri = ee.ImageCollection("projects/sat-io/open-datasets/LCMAP/LCPRI");

// process --------------------------------------------------------------

var lcmap2021 = lcpri
  .filterDate('2021-01-01', '2021-12-31')
  .first();


var vizLcpri = {
  min: 1,
  max: 8,
  palette: [
    'e81414',  // 1 - Developed
    'b68c4c',  // 2 - Cropland
    'd2cda4',  // 3 - Grass/Shrub
    '38814e',  // 4 - Tree Cover
    '4bb0d4',  // 5 - Water
    '85c8c8',  // 6 - Wetlands
    'eeeeee',  // 7 - Ice/Snow
    'b5a27e'   // 8 - Barren
  ]
};
Map.addLayer(lcmap2021, vizLcpri, 'lcmap');
// Binary mask: 1 = keep (grass/shrub, tree, wetlands, ice/snow, barren), 0 = remove (developed, cropland, water)
var mask = lcmap2021.remap(
  [1, 2, 3, 4, 5, 6, 7, 8],
  [0, 0, 1, 1, 0, 1, 1, 1]
);

// Fraction of "keep" pixels within each Daymet cell

var fracKeep = mask
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    bestEffort: true,
    maxPixels: 1e3
  })
  .reproject({
    crs: fg.crs,
    crsTransform: fg.crsTransform,
  });

var fracKeep = fracKeep.rename('fracKeep');
// Map.addLayer(fracKeep, {min: 0, max: 1, palette: ['red', 'white', 'green']}, 'fraction keep');

Map.addLayer(fracKeep.geometry(), {}, '')

Export.image.toAsset({
  image: fracKeep,
  description: 'LCMAP_fracKeep_daymet',
  assetId: 'projects/ee-martinholdrege/assets/PED_vegClimModels2/LCMAP_fracKeep_daymet',
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fracKeep.geometry(), // or define a CONUS geometry
  maxPixels: 1e12
});