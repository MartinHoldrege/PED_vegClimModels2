/*
Export GeoTIFFs of the LCMAP and fire masks, on the snap grid.

Input assets created in:
  - 01_create_lcmap_mask.js        LCMAP_fracKeep_1000m         (1km fracKeep)
  - 01_fire_everBurned_windows.js  MTBS_burned20yr_100m_*       (100m burned,
                                                                mean pyramiding)
  - 02_fire_fracUnburned.js        MTBS_fracUnburned_gte90_*    (1km per-year
                                                                binary)

Every fire layer here is about the preceding 20 years: for target year Y the
window is Y-19 to Y, clamped to the 1984 start of the MTBS record.

Exports:
  - LCMAP mask: 1 where at least 90% of the cell is not developed, crop or
    water
  - fire mask: one band per year ('year_YYYY'), 1 where less than 10% of the
    cell burned in that year's window
  - fire mean mask: 1 where the mean over 2010-2023 of the burned fraction of
    the cell is at most 10%. The 100m asset is read at the 1km grid, where
    mean pyramiding makes each value that cell's burned fraction; those
    fractions are averaged and thresholded once. The per-year asset is
    already thresholded and cannot be used for this.

All three export as byte.

Started: August 2026
*/

// dependencies -------------------------------------

var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------

var yearStart = 2000;
var yearEnd = 2024;
var windowLength = 20;
var scaleIn = 100;          // m, resolution of the burned asset
var maskCutoffLcmap = 0.9;
var maskCutoffFire = 0.9;   // keep where at least this fraction is unburned
 // years over which to average the burned fraction ~ years with most cover data
var meanYearStart = 2010;  
var meanYearEnd = 2023;
var driveFolder = 'PED_vegClimModels2';

// export toggles
var exportLcmapMask = true;
var exportFireMask = true;
var exportFireMean = true;

// LCMAP mask ---------------------------------------

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

// Fire mask, one band per year ---------------------

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

// Fire mask, mean over the target decade -----------

// created in 01_fire_everBurned_windows.js; read at 1km so mean pyramiding
// gives the burned fraction of each cell
var burnedName = 'MTBS_burned' + windowLength + 'yr_' + scaleIn + 'm_' +
  yearStart + '-' + yearEnd;

var burnedFrac = ee.Image(fg.pathAsset + 'fire/' + burnedName)
  .updateMask(fg.maskConus);

var meanBands = [];
for (var y = meanYearStart; y <= meanYearEnd; y++) {
  meanBands.push('year_' + y);
}

var burnedFracMean = burnedFrac.select(meanBands)
  .reduce(ee.Reducer.mean())
  .rename('burnedFracMean');

var fireMeanBinary = burnedFracMean
  .lte(1 - maskCutoffFire)
  .rename('fire_mask')
  .toByte()
  .set({
    windowLength: windowLength,
    maskCutoffFire: maskCutoffFire,
    meanYears: meanYearStart + '-' + meanYearEnd,
    definition: '1 where the mean over ' + meanYearStart + '-' + meanYearEnd +
      ' of the fraction of the cell burned in the preceding ' + windowLength +
      ' years is at most ' + (1 - maskCutoffFire),
    sourceAsset: burnedName
  });

if (exportFireMean) {
  var fireMeanFileName = 'MTBS_burnedFracMean_gte' + maskCutoffFire * 100 +
    '_' + windowLength + 'yr_' + meanYearStart + '-' + meanYearEnd + fg.resLabel;
  fg.exportDrive(fireMeanBinary, fireMeanFileName, driveFolder);
}

// visualize ----------------------------------------

Map.addLayer(fracKeepBinary, {min: 0, max: 1, palette: 'red,white'},
  'lcmap keep', false);
Map.addLayer(fireMask.select('year_' + yearEnd),
  {min: 0, max: 1, palette: 'red,white'}, 'fire keep ' + yearEnd, false);
Map.addLayer(burnedFracMean, {min: 0, max: 0.5, palette: 'white,red'},
  'mean burned fraction ' + meanYearStart + '-' + meanYearEnd, false);
Map.addLayer(fireMeanBinary, {min: 0, max: 1, palette: 'red,white'},
  'fire keep mean', false);