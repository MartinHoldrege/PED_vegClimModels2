/*
App for comparing biomass and cover datasets with adjustable masking.
Split panel with independent controls on each side.
Author: Martin Holdrege
Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var coverYearStart = 2000;
var coverYearEnd = 2023;
var rapBioYearStart = 2000;
var rapBioYearEnd = 2023;
var coverMax = 60;

// data ---------------------------------------------
var coverFileName = 'RAP_v3_cover_' + coverYearStart + '-' + coverYearEnd + fg.resLabel;
var cover = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + coverFileName);

var bioFileName = 'RAP_v3_herbaceousAGB_' + rapBioYearStart + '-' + rapBioYearEnd + fg.resLabel;
var rapHerbAGB = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/rap/' + bioFileName);

var gedi = ee.Image('LARSE/GEDI/GEDI04_B_002').select('MU');
var nbcd = fg.nbcdAGB();
var fracKeep = ee.Image('projects/ee-martinholdrege/assets/PED_vegClimModels2/LCMAP_fracKeep_daymet');

// layer types --------------------------------------
var woodyLayers = ['GEDI AGB (Mg/ha)', 'NBCD AGB (Mg/ha)'];
var herbLayers = ['RAP herb AGB (Mg/ha)'];
var coverLayers = ['RAP tree cover (%)', 'RAP shrub cover (%)', 'RAP herb cover (%)'];

var layerImages = {
  'GEDI AGB (Mg/ha)':          gedi,
  'NBCD AGB (Mg/ha)':          nbcd,
  'RAP herb AGB (Mg/ha)':      rapHerbAGB,
  'RAP tree cover (%)':        cover.select('totalTreeCov'),
  'RAP shrub cover (%)':       cover.select('totalShrubCov'),
  'RAP herb cover (%)':        cover.select('totalHerbaceousCov')
};
var layerNames = Object.keys(layerImages);

// defaults -----------------------------------------
var defaults = {
  minTree: 0,
  maxTree: 100,
  minShrub: 0,
  maxShrub: 100,
  fracKeepThreshold: 0.9,
  woodyAGBmax: 200,
  herbAGBmax: 2
};

// functions ----------------------------------------

var buildMask = function(vals) {
  var m = fracKeep.gte(vals.fracKeepThreshold);
  m = m.and(cover.select('totalTreeCov').gte(vals.minTree));
  m = m.and(cover.select('totalTreeCov').lte(vals.maxTree));
  m = m.and(cover.select('totalShrubCov').gte(vals.minShrub));
  m = m.and(cover.select('totalShrubCov').lte(vals.maxShrub));
  return m;
};

var getViz = function(layerName, vals) {
  if (woodyLayers.indexOf(layerName) >= 0) {
    return {min: 0, max: vals.woodyAGBmax, palette: ['white', 'green']};
  } else if (herbLayers.indexOf(layerName) >= 0) {
    return {min: 0, max: vals.herbAGBmax, palette: ['white', 'green']};
  } else {
    return {min: 0, max: coverMax, palette: ['white', 'green']};
  }
};

var buildPanel = function(map, position) {
  var vals = {
    minTree: defaults.minTree,
    maxTree: defaults.maxTree,
    minShrub: defaults.minShrub,
    maxShrub: defaults.maxShrub,
    fracKeepThreshold: defaults.fracKeepThreshold,
    woodyAGBmax: defaults.woodyAGBmax,
    herbAGBmax: defaults.herbAGBmax
  };
  var selectedLayer = layerNames[0];

  var updateMap = function() {
    var image = layerImages[selectedLayer];
    var mask = buildMask(vals);
    var masked = image.updateMask(mask);
    var viz = getViz(selectedLayer, vals);
    map.layers().reset([ui.Map.Layer(masked, viz, selectedLayer)]);
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
    return ui.Panel(
      [ui.Label(label, {fontSize: '11px', margin: '4px 0px 0px 0px'}), slider],
      ui.Panel.Layout.flow('vertical')
    );
  };

  var makeHeader = function(text) {
    return ui.Label(text, {fontSize: '11px', fontWeight: 'bold', margin: '8px 0px 2px 0px'});
  };

  var panel = ui.Panel({
    widgets: [
      ui.Label('Layer:', {fontSize: '11px', margin: '4px 0px 0px 0px'}),
      dropdown,
      makeHeader('Apply mask to layers based on cover and land use'),
      makeSlider('Min tree cover (%)', 0, 100, 5, defaults.minTree, 'minTree'),
      makeSlider('Max tree cover (%)', 0, 100, 5, defaults.maxTree, 'maxTree'),
      makeSlider('Min shrub cover (%)', 0, 100, 5, defaults.minShrub, 'minShrub'),
      makeSlider('Max shrub cover (%)', 0, 100, 5, defaults.maxShrub, 'maxShrub'),
      makeSlider('Mask (fraction undeveloped)', 0.5, 1, 0.05, defaults.fracKeepThreshold, 'fracKeepThreshold'),
      makeHeader('Change AGB value at which color saturates'),
      makeSlider('Woody AGB max (Mg/ha)', 1, 300, 1, defaults.woodyAGBmax, 'woodyAGBmax'),
      makeSlider('Herb AGB max (Mg/ha)', 1, 10, 1, defaults.herbAGBmax, 'herbAGBmax')
    ],
    style: {position: position, padding: '4px', width: '240px'}
  });

  map.add(panel);
  updateMap();
};

// info panel ---------------------------------------
var infoText =
  'RAP cover: ' + coverYearStart + '–' + coverYearEnd + ' mean\n' +
  'RAP herb AGB: ' + rapBioYearStart + '–' + rapBioYearEnd + ' mean\n' +
  'GEDI: 2019–2023 | NBCD: 2000\n' +
  'Cover: 0–' + coverMax + '%';

var infoPanel = ui.Panel({
  widgets: [ui.Label(infoText, {fontSize: '10px', whiteSpace: 'pre'})],
  style: {position: 'bottom-right', padding: '4px', backgroundColor: 'rgba(255,255,255,0.8)'}
});

// build app ----------------------------------------
ui.root.clear();

var leftMap = ui.Map();
var rightMap = ui.Map();
leftMap.setControlVisibility(true);
rightMap.setControlVisibility(true);

buildPanel(leftMap, 'top-left');
buildPanel(rightMap, 'top-right');
rightMap.add(infoPanel);

var splitPanel = ui.SplitPanel({
  firstPanel: leftMap,
  secondPanel: rightMap,
  wipe: true,
  style: {stretch: 'both'}
});

ui.root.widgets().reset([splitPanel]);
var linker = ui.Map.Linker([leftMap, rightMap]);
leftMap.centerObject(cover, 5);