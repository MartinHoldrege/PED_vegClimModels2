/*
Export GeoTIFFs of the LCMAP and fire masks, on the snap grid.

Input assets created in:
  - LCMAP fracKeep: 01_create_lcmap_mask.js
  - MTBS fracUnburned: 01_fire_everBurned_windows.js, 02_fire_fracUnburned.js

The fire mask has one band per year ('year_YYYY'), where 1 means less than
10% of the cell burned in the preceding 20 years.

Author: Martin Holdrege
Started: August 2026
*/

// dependencies -------------------------------------

var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------

var yearStart = 2000;
var yearEnd = 2024;
var windowLength = 20;
var maskCutoffLcmap = 0.9;
var maskCutoffFire = 0.9;
var driveFolder = 'PED_vegClimModels2';

// export toggles
var exportLcmapMask = true;
var exportFireMask = true;

// LCMAP mask (fracKeep + binary) -------------------

// created in 01_create_lcmap_mask.js
var fracKeep = ee.Image(fg.pathAsset + 'masks/LCMAP_fracKeep' + fg.resLabel)
  .rename('fracKeep')
  .updateMask(fg.maskConus);

var fracKeepBinary = fracKeep
  .gte(maskCutoffLcmap)
  .rename('lcmap_mask')
  .toByte();

if (exportLcmapMask) {
  var lcmapFileName = 'LCMAP_fracKeep_gte' + maskCutoffLcmap * 100 + fg.resLabel;
  fg.exportDrive(fracKeepBinary, lcmapFileName, driveFolder);
}

// Fire mask (one band per year) --------------------

// created in 02_fire_fracUnburned.js; already binary, 1 = keep
var fireMask = ee.Image(fg.pathAsset + 'fire/MTBS_fracUnburned_gte' +
    maskCutoffFire * 100 + '_' + windowLength + 'yr_' +
    yearStart + '-' + yearEnd + fg.resLabel)
  .updateMask(fg.maskConus);

print('fire mask bands', fireMask.bandNames());

if (exportFireMask) {
  var fireFileName = 'MTBS_fracUnburned_gte' + maskCutoffFire * 100 + '_' +
    windowLength + 'yr_' + yearStart + '-' + yearEnd + fg.resLabel;
  fg.exportDrive(fireMask.toByte(), fireFileName, driveFolder);
}

// visualize ----------------------------------------

Map.addLayer(fracKeepBinary, {min: 0, max: 1, palette: 'red,white'},
  'lcmap keep', false);
Map.addLayer(fireMask.select('year_' + yearEnd),
  {min: 0, max: 1, palette: 'red,white'}, 'fire keep ' + yearEnd, false);