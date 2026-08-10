/*
Annual RAP v3 cover (tree, shrub, herbaceous) on the daymet snap grid, masked to
the thinned (every 5th cell) CONUS snap raster. One GeoTIFF per functional
group, one band per year.

Author: Martin Holdrege
Started: August 2026
*/

// params -------------------------------------------
var yearStart = 2011;
var yearEnd = 2023;

var driveFolder = 'PED_vegClimModels2';
var bandPrefix = '';  // set to 'y' if leading-digit band names cause trouble

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// read in data -------------------------------------
// created in R (create_daymet_conus_snap_raster.R); cell ids, masked elsewhere
var snapMask = ee.Image(fg.pathAsset + 'masks/daymet_conus_snap_1000m_thin5')
  .gt(0);


var rap = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));

// process ------------------------------------------
// pull the three functional groups out of a single year's image
var pftBands = function(image) {
  return ee.Image.cat([
    image.select('TRE'),
    image.select('SHR'),
    image.select('AFG').add(image.select('PFG'))
  ]).rename(['tree', 'shrub', 'herbaceous']);
};

// build one image per pft, with a band per year
var pfts = ['tree', 'shrub', 'herbaceous'];
var stacks = {tree: null, shrub: null, herbaceous: null};

for (var year = yearStart; year <= yearEnd; year++) {
  var yearImage = pftBands(
    ee.Image(rap.filter(ee.Filter.calendarRange(year, year, 'year')).first())
  );

  for (var j = 0; j < pfts.length; j++) {
    var pft = pfts[j];
    var band = yearImage.select(pft).rename(bandPrefix + year);
    stacks[pft] = (stacks[pft] === null) ? band : stacks[pft].addBands(band);
  }
}

// export -------------------------------------------
for (var k = 0; k < pfts.length; k++) {
  var pftName = pfts[k];
  var out = stacks[pftName].updateMask(snapMask).toFloat();

  var fileName = 'RAP_v3_cover-' + pftName + '_' + yearStart + '-' + yearEnd +
    '_thin5' + fg.resLabel;

  fg.exportDrive(out, fileName, driveFolder);
}

// visualize ----------------------------------------
Map.addLayer(snapMask.selfMask(), {palette: 'red'}, 'thin5 mask', false);
Map.addLayer(stacks.shrub.select(bandPrefix + yearEnd).updateMask(snapMask),
  {min: 0, max: 40, palette: ['white', 'black']}, 'shrub ' + yearEnd, false);