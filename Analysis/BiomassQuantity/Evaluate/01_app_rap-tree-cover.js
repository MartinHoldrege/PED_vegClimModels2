/*
App for exploring tree cover at native 30m resolution.
Supports RAP and USFS tree canopy cover datasets.
Optional masking by fire (MTBS) and land cover (LCMAP).
Split panel with independent controls on each side.

Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStart = 2019;
var yearEnd = 2023;

// data ---------------------------------------------

// RAP tree cover (mean 2019-2023)
var rapTree = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .select(['TRE'])
  .mean();

// USFS tree canopy cover (mean 2019-2023)
var usfsTree = ee.ImageCollection('USGS/NLCD_RELEASES/2023_REL/TCC/v2023-5')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .filter('study_area == "CONUS"')
  .select(['Science_Percent_Tree_Canopy_Cover'])
  .mean();

var everBurned = ee.Image(fg.pathAsset + 'fire/MTBS_everBurned_30m_2000-2023')
  .unmask(0);

var datasets = {
  'RAP tree cover': rapTree,
  'USFS tree canopy cover': usfsTree
};
var datasetNames = Object.keys(datasets);

// viz --------------

var paletteCov = [
    'CDA066',
    'D7C29E',
    'C2D096',
    'B7D692',
    'ADDD8E',
    '78C679',
    '5CB86B',
    '41AB5D',
    '39A156',
    '329750',
    '238443',
    '11763D',
    '006837',
    '004529'
  ]

// defaults -----------------------------------------
var defaults = {
  minTree: 0,
  maxTree: 100,
  colorMax: 60,
  maskFire: false,
  maskLcmap: false,
  dataset: datasetNames[0]
};

// functions ----------------------------------------

var buildPanel = function(map, position) {
  var vals = {
    minTree: defaults.minTree,
    maxTree: defaults.maxTree,
    colorMax: defaults.colorMax,
    maskFire: defaults.maskFire,
    maskLcmap: defaults.maskLcmap,
    dataset: defaults.dataset
  };

  var updateMap = function() {
    var img = datasets[vals.dataset];

    var mask = ee.Image(1);
    if (vals.maskFire) {
      mask = mask.and(everBurned.not());
    }
    if (vals.maskLcmap) {
      mask = mask.and(fg.lcmapMask);
    }

    var inRange = img.gte(vals.minTree).and(img.lt(vals.maxTree));
    var treeCov = img.updateMask(mask).updateMask(inRange);
    var outOfRange = ee.Image(0).updateMask(mask).updateMask(img.lt(vals.minTree))
      .blend(ee.Image(1).updateMask(mask).updateMask(img.gte(vals.maxTree)));

    var viz = {min: 0, max: vals.colorMax, palette: paletteCov};
    
    map.layers().reset([
      ui.Map.Layer(outOfRange, {min: 0, max: 1, palette: ['black', 'gray']},
        'outside range (' + vals.minTree + '–' + vals.maxTree + '%)'),
      ui.Map.Layer(treeCov, viz, vals.dataset + ' (%)')
    ]);
  };

  var dropdown = ui.Select({
    items: datasetNames,
    value: defaults.dataset,
    style: {stretch: 'horizontal'},
    onChange: function(value) {
      vals.dataset = value;
      updateMap();
    }
  });

  var makeSlider = function(label, min, max, step, defaultVal, key) {
    var slider = ui.Slider({
      min: min, max: max, value: defaultVal, step: step,
      style: {stretch: 'horizontal', padding: '0px 0px 4px 0px'},
      onChange: function(value) {
        vals[key] = value;
        updateMap();
      }
    });
    return ui.Panel(
      [ui.Label(label, {fontSize: '11px', margin: '4px 0px 0px 0px'}), slider],
      ui.Panel.Layout.flow('vertical')
    );
  };

  var maskFireCheck = ui.Checkbox({
    label: 'Mask burned areas (MTBS 2000–2023)',
    value: defaults.maskFire,
    style: {fontSize: '11px'},
    onChange: function(value) {
      vals.maskFire = value;
      updateMap();
    }
  });

  var maskLcmapCheck = ui.Checkbox({
    label: 'Mask developed/cropland/water (LCMAP)',
    value: defaults.maskLcmap,
    style: {fontSize: '11px'},
    onChange: function(value) {
      vals.maskLcmap = value;
      updateMap();
    }
  });

  var panel = ui.Panel({
    widgets: [
      ui.Label('Tree cover ' + yearStart + '–' + yearEnd,
        {fontWeight: 'bold', fontSize: '12px'}),
      ui.Label('Dataset:', {fontSize: '11px', margin: '4px 0px 0px 0px'}),
      dropdown,
      makeSlider('Min tree cover (%; below is black)', 0, 100, 1, defaults.minTree, 'minTree'),
      makeSlider('Max tree cover (%; above is gray)', 0, 100, 1, defaults.maxTree, 'maxTree'),
      makeSlider('Color saturation (%)', 1, 100, 1, defaults.colorMax, 'colorMax'),
      ui.Label('Masking', {fontWeight: 'bold', fontSize: '11px', margin: '8px 0px 2px 0px'}),
      maskFireCheck,
      maskLcmapCheck
    ],
    style: {position: position, padding: '4px', width: '260px'}
  });

  map.add(panel);
  updateMap();
};

// build app ----------------------------------------
ui.root.clear();

var leftMap = ui.Map();
var rightMap = ui.Map();
leftMap.setControlVisibility(true);
rightMap.setControlVisibility(true);

buildPanel(leftMap, 'top-left');
buildPanel(rightMap, 'top-right');

var splitPanel = ui.SplitPanel({
  firstPanel: leftMap,
  secondPanel: rightMap,
  wipe: true,
  style: {stretch: 'both'}
});

ui.root.widgets().reset([splitPanel]);
var linker = ui.Map.Linker([leftMap, rightMap]);
leftMap.setCenter(-110, 40, 5);

// info
var infoPanel = ui.Panel({
  widgets: [ui.Label(
    'Tree cover (30m), ' + yearStart + '–' + yearEnd + ' mean\n' +
    'RAP: Rangeland Analysis Platform v3\n' +
    'USFS: NLCD Tree Canopy Cover v2023-5',
    {fontSize: '10px', whiteSpace: 'pre'})],
  style: {position: 'bottom-right', padding: '4px',
    backgroundColor: 'rgba(255,255,255,0.8)'}
});
rightMap.add(infoPanel);