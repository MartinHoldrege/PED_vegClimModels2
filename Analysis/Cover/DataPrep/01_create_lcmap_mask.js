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

// Binary mask: 1 = keep (grass/shrub, tree, wetlands, ice/snow, barren), 0 = remove (developed, cropland, water)
var mask = fg.lcmapMask; // at native resolution

// process ------------------------------------------

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

Map.addLayer(fracKeep.geometry(), {}, '');

fg.exportAsset(fracKeep, 'LCMAP_fracKeep_daymet' + fg.resLabel, 'masks/');



