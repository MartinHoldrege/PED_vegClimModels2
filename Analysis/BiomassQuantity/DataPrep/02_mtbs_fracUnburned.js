/*
fraction of daymet cell that has burned in given time period
Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------

var yearStart = 2000;
var yearEnd = 2023;

// created in 01_mtbs_everburned.js
var everBurned = ee.Image(fg.pathAsset + 'fire/MTBS_everBurned_30m_' + yearStart + '-' + yearEnd);

// everBurned * pixelArea, then unmask to 0 so unburned/unmapped areas contribute 0 area
var areaBurned = everBurned.multiply(ee.Image.pixelArea()).unmask(0);

var fracBurned = areaBurned.divide(ee.Image.pixelArea());
var fracUnburned = ee.Image(1).subtract(fracBurned).rename('fracUnburned');

Map.addLayer(fracUnburned, {min:0, max: 1, palette: ['black', 'white']});

// export -------------------------------------------

var fileName = 'MTBS_fracUnburned_' + yearStart + '-' + yearEnd + fg.resLabel;

Export.image.toAsset({
  image: fracUnburned,
  description: fileName,
  assetId: fg.pathAsset + 'fire/' + fileName,
  crs: fg.crs,
  crsTransform: fg.crsTransform,
  region: fg.region,
  maxPixels: 1e12
});


