



// grid definition ----------------------------------
// Everything in this repo is exported on the grid of the snap raster.
// These values are read off the asset itself (print(fg.maskConus) -> bands[0]).
// Do not hardcode a grid anywhere else.

exports.pathAsset = 'projects/ee-martinholdrege/assets/PED_vegClimModels2/';

exports.maskConus = ee.Image(exports.pathAsset + 'masks/daymet_conus_snap_1000m')
  .rename('cell_id');  // asset ingests as 'b1'

exports.crs = "PROJCS[\"unnamed\",    GEOGCS[\"unknown\",      DATUM[\"unknown\",        SPHEROID[\"Spheroid\", 6378137.0, 298.257223563]],      PRIMEM[\"Greenwich\", 0.0],      UNIT[\"degree\", 0.017453292519943295],      AXIS[\"Longitude\", EAST],      AXIS[\"Latitude\", NORTH]],    PROJECTION[\"Lambert_Conformal_Conic_2SP\"],    PARAMETER[\"central_meridian\", -100.0],    PARAMETER[\"latitude_of_origin\", 42.5],    PARAMETER[\"standard_parallel_1\", 60.0],    PARAMETER[\"false_easting\", 0.0],    PARAMETER[\"false_northing\", 0.0],    PARAMETER[\"scale_factor\", 1.0],    PARAMETER[\"standard_parallel_2\", 25.0],    UNIT[\"m\", 1.0],    AXIS[\"Easting\", EAST],    AXIS[\"Northing\", NORTH]]";

exports.crsTransform = [1000, 0, -1951750, 0, -1000, 946500];

exports.resolution = 1000;
exports.resLabel = '_1000m';

// region = the snap raster footprint, in the snap raster's own projection
exports.region = ee.Geometry.Rectangle({
  coords: [-1951750, -1786500, 2429250, 946500],
  proj: exports.crs,
  geodesic: false
});

// export helpers -----------------------------------
// Use these instead of calling Export.image.* directly, so that crs,
// crsTransform, and region can never diverge between scripts.

/**
 * Export an image to Drive on the snap grid.
 * @param {ee.Image} image
 * @param {string} name - used for both description and file name
 * @param {string} folder - Drive folder
 */
exports.exportDrive = function(image, name, folder) {
  Export.image.toDrive({
    image: image,
    description: name,
    folder: folder,
    fileNamePrefix: name,
    crs: exports.crs,
    crsTransform: exports.crsTransform,
    region: exports.region,
    maxPixels: 1e12,
    fileFormat: 'GeoTIFF'
  });
};

/**
 * Export an image to an asset on the snap grid.
 * @param {ee.Image} image
 * @param {string} name
 * @param {string} subDir - subdirectory under pathAsset, e.g. 'rap/'
 */
exports.exportAsset = function(image, name, subDir) {
  Export.image.toAsset({
    image: image,
    description: name,
    assetId: exports.pathAsset + subDir + name,
    crs: exports.crs,
    crsTransform: exports.crsTransform,
    region: exports.region,
    maxPixels: 1e12
  });
};
// biomass related functions -------------------------------------

// rap biomass code adapted from: 
// https://code.earthengine.google.com/d172689436c5d4a1bc6bd5e64f52784a

/**
 * Convert a single RAP NPP image to aboveground herbaceous biomass (Mg/ha).
 * @param {ee.Image} image - Two-band image (afgNPP, pfgNPP)
 * @returns {ee.Image} Three bands: afgAGB, pfgAGB, herbaceousAGB (Mg/ha)
 */
var nppToBiomass = function(image) {
  var mat = ee.ImageCollection('projects/rap-data-365417/assets/gridmet-MAT');
  var year = ee.Date(image.get('system:time_start')).format('YYYY');
  var matYear = mat.filterDate(year).first();
  var fANPP = matYear.multiply(0.0129).add(0.171).rename('fANPP');
  
  // 0.0001: NPP scalar -> kg C/m2
  // fANPP: fraction aboveground
  // 2.1276: C to biomass
  // 10: kg/m2 to Mg/ha
  var agb = image.multiply(0.0001)
               .multiply(fANPP)
               .multiply(2.1276)
               .multiply(10)
               .rename(['afgAGB', 'pfgAGB'])
               .copyProperties(image, ['system:time_start'])
               .set('year', year);
  agb = ee.Image(agb);
  var herbaceous = agb.select('afgAGB').add(agb.select('pfgAGB')).rename('herbaceousAGB');
  return agb.addBands(herbaceous);
};

/**
 * Mean RAP herbaceous aboveground biomass (Mg/ha) over a year range.
 * @param {number} yearStart - First year (inclusive)
 * @param {number} yearEnd - Last year (inclusive)
 * @returns {ee.Image} Single band: herbaceousAGB (Mg/ha)
 */
exports.rapHerbBiomass = function(yearStart, yearEnd) {
  var npp = ee.ImageCollection('projects/rap-data-365417/assets/npp-partitioned-v3')
    .select(['afgNPP', 'pfgNPP'])
    .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'));
  
  return npp.map(nppToBiomass).select('herbaceousAGB').mean();
};



/**
 * NBCD aboveground biomass converted to Mg/ha.
 * Raw asset is in total Mg per 240m cell; divide by 5.76 ha/cell.
 * Mean pyramiding preserves this conversion at coarser scales.
 * @returns {ee.Image} Single band: AGB (Mg/ha)
 */
exports.nbcdAGB = function() {
  return ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/biomass/NBCD_countrywide_biomass_mosaic')
    .divide(5.76)
    .rename('AGB');
};

// lcmap mask

var lcmap2021 = ee.ImageCollection('projects/sat-io/open-datasets/LCMAP/LCPRI')
  .filterDate('2021-01-01', '2021-12-31')
  .first();

// process ------------------------------------------

// LCMAP mask: 1 = keep (grass/shrub, tree, wetlands, ice/snow, barren)
// 0 = remove (developed, cropland, water)
exports.lcmapMask = lcmap2021.remap(
  [1, 2, 3, 4, 5, 6, 7, 8],
  [0, 0, 1, 1, 0, 1, 1, 1]
);

var fg = exports
// force a planar rectangle in the daymet CRS
print(fg.maskConus)
print('bounds', fg.maskConus.geometry(1, fg.crs, false).bounds(1, fg.crs));
print('transform', fg.maskConus.projection().transform());