

// CRS for daymet data
// exports.crs="PROJCRS[\"unnamed\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"unknown\",\n            ELLIPSOID[\"Spheroid\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1,\n                    ID[\"EPSG\",9001]]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Lambert Conic Conformal (2SP)\",\n        METHOD[\"Lambert Conic Conformal (2SP)\",\n            ID[\"EPSG\",9802]],\n        PARAMETER[\"Latitude of false origin\",42.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",25,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",60,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]";

exports.resolution = 1000; // resolution of daymet data

var projection = ee.ImageCollection('NASA/ORNL/DAYMET_V4')
                  .first()
                  .projection();

exports.crs = "PROJCS[\"unnamed\",    GEOGCS[\"NAD83\",      DATUM[\"WGS_1984\",        SPHEROID[\"WGS 84\", 6378137.0, 298.257223563, AUTHORITY[\"EPSG\",\"7030\"]],        AUTHORITY[\"EPSG\",\"6326\"]],      PRIMEM[\"Greenwich\", 0.0],      UNIT[\"degree\", 0.017453292519943295],      AXIS[\"Longitude\", EAST],      AXIS[\"Latitude\", NORTH]],    PROJECTION[\"Lambert_Conformal_Conic_2SP\"],    PARAMETER[\"central_meridian\", -100.0],    PARAMETER[\"latitude_of_origin\", 42.5],    PARAMETER[\"standard_parallel_1\", 60.0],    PARAMETER[\"false_easting\", 0.0],    PARAMETER[\"false_northing\", 0.0],    PARAMETER[\"scale_factor\", 1.0],    PARAMETER[\"standard_parallel_2\", 25.0],    UNIT[\"m\", 1.0],    AXIS[\"x\", EAST],    AXIS[\"y\", NORTH]]";
var crsTransform = [1000, 0, -5802750, 0, -1000, 4984500];
exports.crsTransform = crsTransform;
     
     
exports.region = ee.ImageCollection("projects/sat-io/open-datasets/LCMAP/LCPRI")
  .first()
  .geometry();
  
var scale = Math.abs(crsTransform[0]);
exports.resLabel = '_' + scale + 'm';  // '1000m'


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





