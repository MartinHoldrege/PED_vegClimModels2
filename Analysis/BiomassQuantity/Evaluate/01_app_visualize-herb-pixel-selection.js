/*
App for exploring RAP cover thresholds for identifying herb-dominant pixels.
Split panel with independent sliders on each side.

Started: April 6, 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// data ---------------------------------------------
var yearStart = 2019;
var yearEnd = 2023;
var fileName = 'RAP_v3_cover_' + yearStart + '-' + yearEnd + fg.resLabel;
var cover = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + fileName);
var fracKeep = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/LCMAP_fracKeep_daymet');

// viz params ---------------------------------------
var vizCover = {min: 0, max: 60, palette: ['white', 'green']};
var vizMask = {min: 0, max: 1, palette: ['white', 'black']};
var vizBinary = {palette: ['lime']};

// defaults -----------------------------------------
var defaults = {
  treeCutoff: 5,
  shrubCutoff: 5,
  herbThreshold: 5,
  fracKeepThreshold: 0.9
};

// layer names (dropdown options)
var layerNames = [
  'Herb-dominant (all criteria)',
  'Tree cover',
  'Shrub cover',
  'Herb cover',
  'Tree < cutoff',
  'Shrub < cutoff',
  'Herb >= threshold',
  'Mask'
];

// functions ----------------------------------------

var computeLayer = function(name, vals) {
  var mask = fracKeep.gte(vals.fracKeepThreshold);
  var coverMasked = cover.updateMask(mask);
  
  if (name === 'Tree cover') {
    return {image: coverMasked.select('totalTreeCov'), viz: vizCover};
  } else if (name === 'Shrub cover') {
    return {image: coverMasked.select('totalShrubCov'), viz: vizCover};
  } else if (name === 'Herb cover') {
    return {image: coverMasked.select('totalHerbaceousCov'), viz: vizCover};
  } else if (name === 'Tree < cutoff') {
    return {image: coverMasked.select('totalTreeCov').lt(vals.treeCutoff).selfMask(), viz: vizBinary};
  } else if (name === 'Shrub < cutoff') {
    return {image: coverMasked.select('totalShrubCov').lt(vals.shrubCutoff).selfMask(), viz: vizBinary};
  } else if (name === 'Herb >= threshold') {
    return {image: coverMasked.select('totalHerbaceousCov').gte(vals.herbThreshold).selfMask(), viz: vizBinary};
  } else if (name === 'Mask') {
    return {image: mask.selfMask(), viz: vizMask};
  } else {
    var result = coverMasked.select('totalTreeCov').lt(vals.treeCutoff)
      .and(coverMasked.select('totalShrubCov').lt(vals.shrubCutoff))
      .and(coverMasked.select('totalHerbaceousCov').gte(vals.herbThreshold))
      .selfMask();
    return {image: result, viz: vizBinary};
  }
};

var buildPanel = function(map, position) {
  var vals = {
    treeCutoff: defaults.treeCutoff,
    shrubCutoff: defaults.shrubCutoff,
    herbThreshold: defaults.herbThreshold,
    fracKeepThreshold: defaults.fracKeepThreshold
  };
  var selectedLayer = layerNames[0];
  
  var updateMap = function() {
    var layer = computeLayer(selectedLayer, vals);
    map.layers().reset([ui.Map.Layer(layer.image, layer.viz, selectedLayer)]);
  };
  
  var dropdown = ui.Select({
    items: layerNames,
    value: selectedLayer,
    style: {stretch: 'horizontal'},
    onChange: function(value) {
      selectedLayer = value;
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
    return ui.Panel([ui.Label(label, {fontSize: '11px', margin: '4px 0px 0px 0px'}), slider],
      ui.Panel.Layout.flow('vertical'));
  };
  
  var panel = ui.Panel({
    widgets: [
      ui.Label('RAP cover ' + yearStart + '–' + yearEnd, {fontWeight: 'bold', fontSize: '12px'}),
      ui.Label('Display:', {fontSize: '11px', margin: '4px 0px 0px 0px'}),
      dropdown,
      makeSlider('Max tree cover (%)', 0, 20, 1, defaults.treeCutoff, 'treeCutoff'),
      makeSlider('Max shrub cover (%)', 0, 20, 1, defaults.shrubCutoff, 'shrubCutoff'),
      makeSlider('Min herb cover (%)', 0, 20, 1, defaults.herbThreshold, 'herbThreshold'),
      makeSlider('Mask (fraction undeveloped)', 0.5, 1, 0.05, defaults.fracKeepThreshold, 'fracKeepThreshold')
    ],
    style: {position: position, padding: '4px', width: '240px'}
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
leftMap.centerObject(cover, 5);

// info label
var infoLabel = ui.Label(
  'RAP v3 cover data averaged to ' + fg.resLabel + ' for ' + yearStart + '–' + yearEnd,
  {fontSize: '11px', position: 'bottom-center', backgroundColor: 'rgba(255,255,255,0.7)'}
);
leftMap.add(infoLabel);