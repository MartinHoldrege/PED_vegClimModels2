/*
App for exploring tree cover at native 30m resolution.
Supports RAP and USFS tree canopy cover datasets.
Optional masking by fire (MTBS) and land cover (LCMAP).
Split panel with independent controls on each side.
CDF chart compares left and right panel datasets.

Started: April 2026
*/

// dependencies -------------------------------------
var fg = require('users/MartinHoldrege/PED_vegClimModels2:Functions/gee/general.js');

// params -------------------------------------------
var yearStart = 2019;
var yearEnd = 2023;

// data ---------------------------------------------

var rapTree = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .select(['TRE'])
  .mean();

var usfsTree = ee.ImageCollection('USGS/NLCD_RELEASES/2023_REL/TCC/v2023-5')
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .filter('study_area == "CONUS"')
  .select(['Science_Percent_Tree_Canopy_Cover', 'NLCD_Percent_Tree_Canopy_Cover'])
  .mean();

var everBurned = ee.Image(fg.pathAsset + 'fire/MTBS_everBurned_30m_2000-2023')
  .unmask(0);
  
var percTreed = function(cutoff) {
  // created in 03_rap_fracNotForest.js
  var fileName = 'RAP_v3_fracNotForest_lt' + cutoff + '_' +
    yearStart+ '-' + yearEnd + fg.resLabel;
    
  var fracNotForest = ee.Image(fg.pathAsset + 'rap/' + fileName);
  var mask = fracNotForest.mask()
    .neq(ee.Image(0)); // becomes 0 where should mask. want to avoid partial masking
  
  // to keep plotting easier on the scale of 0-100% like tree cover,
  // converting to percent >= cutoff. 
  var percTreed = ee.Image(1).subtract(fracNotForest).multiply(100)
    .unmask()
    .updateMask(mask);
  return percTreed;
};

var datasets = {
  'RAP tree cover': rapTree,
  'USFS tree canopy cover (Science)': usfsTree.select('Science_Percent_Tree_Canopy_Cover'),
  'USFS tree canopy cover (NLCD)': usfsTree.select('NLCD_Percent_Tree_Canopy_Cover'),
  'RAP % of 1km w/ >= 1% trees': percTreed(1),
  'RAP % of 1km w/ >= 3% trees': percTreed(3),
  'RAP % of 1km w/ >= 5% trees': percTreed(5),
  'RAP % of 1km w/ >= 10% trees': percTreed(10),
};

var datasetNames = Object.keys(datasets);

// viz --------------
var paletteCov = [
  'CDA066', 'D7C29E', 'C2D096', 'B7D692', 'ADDD8E', '78C679', '5CB86B',
  '41AB5D', '39A156', '329750', '238443', '11763D', '006837', '004529'
];

// defaults -----------------------------------------
var defaults = {
  minTree: 0,
  maxTree: 100,
  colorMax: 60,
  maskFire: false,
  maskLcmap: false,
  dataset: datasetNames[0]
};

// shared state -------------------------------------
var leftVals = null;
var rightVals = null;
var leftMap = ui.Map();
var rightMap = ui.Map();
var chartPanel = null;

// chart helpers ------------------------------------

var getMasked = function(vals) {
  var img = datasets[vals.dataset];
  var mask = ee.Image(1);
  if (vals.maskFire) mask = mask.and(everBurned.not());
  if (vals.maskLcmap) mask = mask.and(fg.lcmapMask);
  return img.updateMask(mask);
};

var parseResult = function(result) {
  var values = [];
  var keys = Object.keys(result);
  var prefix = keys[0].split('_p')[0]; // e.g. 'TRE'
  for (var i = 0; i <= 100; i++) {
    var key = prefix + '_p' + i;
    values.push(result[key] !== undefined ? result[key] : null);
  }
  return values;
};

var updateChart = function() {
  if (chartPanel) {
    leftMap.remove(chartPanel);
  }
  chartPanel = ui.Panel({
    style: {position: 'bottom-left', width: '350px', height: '250px',
      backgroundColor: 'rgba(255,255,255,0.9)'}
  });
  leftMap.add(chartPanel);
  chartPanel.add(ui.Label('Computing ...', {fontSize: '10px'}));

  var pctList = ee.List.sequence(0, 100, 1);
  var region = leftMap.getBounds(true);
  var leftMasked = getMasked(leftVals);
  var rightMasked = getMasked(rightVals);

  var leftPct = leftMasked.reduceRegion({
    reducer: ee.Reducer.percentile(pctList),
    geometry: region, scale: 250, maxPixels: 1e8,
    bestEffort: true
  });
  var rightPct = rightMasked.reduceRegion({
    reducer: ee.Reducer.percentile(pctList),
    geometry: region, scale: 250, maxPixels: 1e8,
    bestEffort: true
  });

  ee.Dictionary(leftPct).evaluate(function(leftResult) {
    ee.Dictionary(rightPct).evaluate(function(rightResult) {
      if (!leftResult || !rightResult) {
        chartPanel.clear();
        chartPanel.add(ui.Label('No data in view', {fontSize: '10px'}));
        return;
      }
      var leftY = parseResult(leftResult);
      var rightY = parseResult(rightResult);
      var rows = [];
      for (var i = 0; i <= 100; i++) {
        rows.push([i, leftY[i], rightY[i]]);
      }
      var dataTable = [
        ['Percentile', 'Left: ' + leftVals.dataset, 'Right: ' + rightVals.dataset]
      ].concat(rows);
      var chart = ui.Chart(dataTable)
        .setChartType('LineChart')
        .setOptions({
          title: 'Tree cover CDF (~ >= 250m, current view)',
          hAxis: {title: 'Percentile'},
          vAxis: {title: 'Tree cover (%)'},
          lineWidth: 2,
          colors: ['blue', 'red']
        });
      chartPanel.clear();
      chartPanel.add(chart);
    });
  });
};
// build panel function -----------------------------

var buildPanel = function(map, position) {
  var vals = {
    minTree: defaults.minTree,
    maxTree: defaults.maxTree,
    colorMax: defaults.colorMax,
    maskFire: defaults.maskFire,
    maskLcmap: defaults.maskLcmap,
    dataset: defaults.dataset
  };

  if (position === 'top-left') {
    leftVals = vals;
  } else {
    rightVals = vals;
  }

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
      .blend(ee.Image(1).updateMask(mask).updateMask(img.gt(vals.maxTree)));

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

  var widgets = [
    ui.Label('Tree cover ' + yearStart + '–' + yearEnd,
      {fontWeight: 'bold', fontSize: '12px'}),
    ui.Label('Dataset:', {fontSize: '11px', margin: '4px 0px 0px 0px'}),
    dropdown,
    makeSlider('Min tree cover (%), or min % of 1km >= x% trees (below is black)', 0, 100, 1, defaults.minTree, 'minTree'),
    makeSlider('Max tree cover (%), or max % of 1km >= x% trees (above is black)', 0, 100, 1, defaults.maxTree, 'maxTree'),
    makeSlider('Color saturation (%)', 1, 100, 1, defaults.colorMax, 'colorMax'),
    ui.Label('Masking', {fontWeight: 'bold', fontSize: '11px', margin: '8px 0px 2px 0px'}),
    maskFireCheck,
    maskLcmapCheck
  ];

  if (position === 'top-right') {
    var chartButton = ui.Button({
      label: 'Show ECDF (current view, \nbased on pixels in current view)',
      style: {stretch: 'horizontal'},
      onClick: updateChart
    });
    widgets.push(chartButton);
  }

  var panel = ui.Panel({
    widgets: widgets,
    style: {position: position, padding: '4px', width: '260px'}
  });

  map.add(panel);
  updateMap();
};

// build app ----------------------------------------
ui.root.clear();

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
    '% of 1 km pixel that is >=x% RAP trees is applied only to\n' +
    '30m pixels that do not have the LCMAP mask applied.\n' +
    'USFS: NLCD Tree Canopy Cover v2023-5',
    {fontSize: '10px', whiteSpace: 'pre'})],
  style: {position: 'bottom-right', padding: '4px',
    backgroundColor: 'rgba(255,255,255,0.8)'}
});
rightMap.add(infoPanel);