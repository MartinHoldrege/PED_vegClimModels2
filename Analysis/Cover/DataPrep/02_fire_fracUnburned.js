/*
Fire mask on the daymet 1km grid: 1 where less than 10% of the cell burned in
the preceding 20 years, 0 otherwise. One band per year, 2000-2024.

Reads the 100m burned layers created in 01_mtbs_burned20yr.js. Because those
were exported with mean pyramiding, reading them at the 1km grid returns the
burned proportion, which is then thresholded.

Author: Martin Holdrege
Started: August 2026
*/


// params -------------------------------------------
var yearStart = 2000;
var yearEnd = 2024;
var windowLength = 20;
var scaleIn = 100;        // m, resolution of the input asset
var burnedThreshold = 0.1; // keep cells with less than this burned fraction

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -------------------------------------
// objected created by Cover/DataPrep/01_fire_everBurned_windows.js
var inName = 'MTBS_burned' + windowLength + 'yr_' + scaleIn + 'm_' +
  yearStart + '-' + yearEnd;

var burnedFrac = ee.Image(fg.pathAsset + 'fire/' + inName);

print('input bands', burnedFrac.bandNames());

// process ------------------------------------------
// 1 = keep (burned fraction below threshold), 0 = drop
var mask = burnedFrac.lt(burnedThreshold).toByte()
  .set({
    windowLength: windowLength,
    burnedThreshold: burnedThreshold,
    definition: '1 = less than ' + (burnedThreshold * 100) +
      '% of cell burned in preceding ' + windowLength + ' years',
    sourceAsset: inName
  });

// visualize ----------------------------------------
Map.addLayer(mask.select('year_' + yearEnd),
  {min: 0, max: 1, palette: ['red', 'white']},
  'unburned mask ' + yearEnd, false);

// export -------------------------------------------
var fileName = 'MTBS_fracUnburned_gte' + ((1 - burnedThreshold) * 100) + '_' +
  windowLength + 'yr_' + yearStart + '-' + yearEnd + fg.resLabel;

fg.exportAsset(mask, fileName, 'fire/');