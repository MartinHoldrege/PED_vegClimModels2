

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
