//All the outputs of this code are stored in assets and used in anomaly prediction
// and vulnerability mapping.
// The required input is - 'region of interest' geometry

Map.addLayer(pa)
Map.centerObject(pa,8)

//Mask to filter clouds and to scale the DN values to surface reflectance
function maskL8sr(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloudnb
  // Bit 2 - Cirrus
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);

  // Replace the original bands with the scaled ones and apply the masks.
  return image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask);
}

// Function to compute spectral indices
var addIndices = function(image) {

  var mirbi = image.expression(
    '(10*swir1) - (9.8*swir2) + 2',{
    'swir1':image.select('SR_B6'),
    'swir2':image.select('SR_B7')
      
    }).rename('mirbi')
    

  var evi = image.expression(
    '2.5 * ((NIR - RED)/(NIR + 6*RED - 7.5*BLUE + 1))', {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4'),
      'BLUE': image.select('SR_B2')
    }).rename('evi');

  var bsi = image.expression(
      '(( X + Y ) - (A + B)) /(( X + Y ) + (A + B)) ', {
        'X': image.select('SR_B6'),
        'Y': image.select('SR_B4'),
        'A': image.select('SR_B5'),
        'B': image.select('SR_B2'),
    }).rename('bsi');
  var ndmi = image.normalizedDifference(['SR_B5', 'SR_B6'])
    .rename('ndmi');
  var nbr = image.normalizedDifference(['SR_B5', 'SR_B7'])
    .rename('nbr');
  var nbr2 = image.normalizedDifference(['SR_B6', 'SR_B7'])
    .rename('nbr2');
  
  var msavi = image.expression(
    '((2 * NIR + 1) - sqrt((2 * NIR + 1) ** 2 - 8 * (NIR - RED))) / 2',
    {
      'NIR': image.select('SR_B5'),
      'RED': image.select('SR_B4')
      
    }
  ).rename('msavi');

  return image
    .addBands(mirbi)
    .addBands(evi)
    .addBands(bsi)
    .addBands(ndmi)
    .addBands(nbr)
    .addBands(nbr2)
    .addBands(msavi);
};




// VIIRS data is available only from 20th Jan
// hence considering landsat also from 20th jan. Landsat 8 data is available from 2013.
// Preparing the trend data for landsat 8-the data from 2013 to Feb 2023 
// 
var outlierDate = ee.Date('2016-06-22');
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                     .filterDate('2013-03-20', '2023-02-28').filterBounds(pa)
                     .filterMetadata('CLOUD_COVER', 'less_than', 10)
                     .map(maskL8sr)
                     .filter(ee.Filter.or(
                    ee.Filter.date('2013-03-20', outlierDate),  // Keep images before outlier
                    ee.Filter.date(outlierDate.advance(1, 'day'), '2023-02-28') // Keep images after outlier
                     ))
                     .map(function(image) {
                    return image.clip(pa);
                  })
                  .map(addIndices);


// Function to calculate SMI for the collection
var addSMI = function(image) {
  // Access globally defined min/max values
  var smi = image.select('SR_B7').subtract(swirMin)
                 .divide(swirMax.subtract(swirMin))
                 .multiply(-1).add(1).rename('smi');
  return image.addBands(smi);
};
// Calculate SWIR min/max for the entire collection
var swirStats = (landsat.mean().select('SR_B7')).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});

var swirMin = ee.Number(swirStats.get('SR_B7_min'));
var swirMax = ee.Number(swirStats.get('SR_B7_max'));

landsat = landsat.map(addSMI)

landsat = landsat.select([ 'evi', 'mirbi','bsi','ndmi','nbr','nbr2','msavi','smi','ST_B10']);
landsat = landsat.sort('system:time_start')
print(landsat,'landsat')


// Add a time band (in years) with explicit casting
var addTimeBand = function(image) {
  var years = ee.Date(image.get('system:time_start')).difference(ee.Date('2013-03-20'), 'year');
  return image.addBands(ee.Image.constant(years).float().rename('time')); // Explicitly cast to float
};

var collectionWithTime = landsat.map(addTimeBand);

// List of indices to process
var indices = ['ndvi', 'evi', 'mirbi','ndfi','bsi','ndmi','nbr','nbr2','msavi','smi','ST_B10'];


// Function to calculate trends for a single index
var calculateTrend = function(index) {
  // Select the index and time bands
  var stacked = collectionWithTime.select(['time', index]);
  // Perform linear regression
  var regression = stacked.reduce(ee.Reducer.linearFit());
  // Extract slope and intercept
  var slope = regression.select('scale').rename(index + '_Slope'); // Trend slope
  var intercept = regression.select('offset').rename(index + '_Intercept'); // Baseline
  return slope.addBands(intercept);
};

// Calculate trends for all indices and combine results
var trendResults = indices.map(calculateTrend);



// Export each trend result separately
trendResults.forEach(function(image, index) {
  var indexName = indices[index]; // Name of the index (e.g., NDVI, EVI)
  Export.image.toAsset({
    image: ee.Image(image),
    description: 'Trend2023_' + indexName,
    assetId: 'path_to_export' + indexName,
    region: pa,
    scale: 30,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
});










///Precipitation trend analysis
// Load CHIRPS daily data for the time period 1982-2022
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filter(ee.Filter.date('1982-01-01', '2022-12-31'))
              .filterBounds(pa);

// Function to compute the yearly mean and include the year as metadata
var createYearlySum = function(year) {
  // Filter for the specific year
  var yearFilter = chirps.filter(ee.Filter.calendarRange(year, year, 'year'));

  // Create the yearly sum (sum precipitation for the year)
  var yearlySum = yearFilter.sum().clip(pa);

  // Add the year as a property
  return yearlySum.set('year', year).rename('precipitation')
                  .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
};

// Create a list of yearly means for each year from 1982 to 2024
var years = ee.List.sequence(1982, 2021);  // Define the range of years
var yearlySums = years.map(createYearlySum);

var prepVis = {bands: ['precipitation'],palette:['#d3ffd2', '#99FF33', 
'#FFFF00', '#FF9900', '#FF3300', '#FF0000'],
 min: 1000, max: 5000}


// Convert the list of yearly means into an ImageCollection
var yearlySumsCollection = ee.ImageCollection(yearlySums);
print(yearlySumsCollection)
//Map.addLayer(yearlyMeansCollection.first(),prepVis,'precipitation')


// Add a 'time' band to the collection for trend calculation
var collectionWithTime = yearlySumsCollection.map(function(image) {
  var year = ee.Date(image.get('system:time_start')).difference(ee.Date('1982-01-01'), 'year');
  return image.addBands(ee.Image.constant(year).rename('time').float());
});

// Function to calculate the trend for precipitation
var calculatePrecipitationTrend = function() {
  var stacked = collectionWithTime.select(['time', 'precipitation']);
  var regression = stacked.reduce(ee.Reducer.linearFit()); // Perform linear regression
  var slope = regression.select('scale').rename('rain_Slope'); // Trend slope
  var intercept = regression.select('offset').rename('rain_Intercept'); // Baseline
  return slope.addBands(intercept); // Combine slope and intercept
};

// Calculate precipitation trend
var precipitationTrend = calculatePrecipitationTrend();
print('Precipitation Trend:', precipitationTrend);

// Visualization parameters
var slopeVis = {min: -6, max: 30, palette: ['blue', 'white', 'red']};
var interceptVis = {min: 0, max: 3000, palette: ['white', 'green']};


// Export the trend results
Export.image.toAsset({
  image: precipitationTrend,
  description: 'Precipitation_Trend_1982_2022',
  assetId: 'path_to_export'+'/file_name',
  region: pa, // Replace with your region of interest
  scale: 30, // Adjust scale to match data resolution
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


// Create single-day average of SMAP L4 version 7, 3-hourly soil moisture
// to view. Plot surface soil moisture and root zone soil moisture over
// user-determined time period.
var soilMoisture = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005').filterBounds(pa)
                        .map(function(image){
                        var x = image.clip(pa)
                        return x;
                        }).filterDate('2015-04-01', '2023-02-28')
                            .select(['soil_moisture_am']);

// Define visualization parameters.
var soilMoistureVis = {
  min: 0.0,
  max: 0.7,
  palette: ['A67C00', 'FFE625', 'C2E5D3', '90DCD0',
            '2FBDBD', '0C9BBD', '068682'],
};


// Function to add a time band to an ImageCollection
var addTimeBand = function(collection) {
  return collection.map(function(image) {
    var time = ee.Date(image.get('system:time_start')).difference(ee.Date('2015-04-01'), 'year');
    return image.addBands(ee.Image.constant(time).rename('time').float());
  });
};

// **1. Soil Moisture Trend**
var soilMoistureTrend = function() {
  // Add time band
  var soilMoistureWithTime = addTimeBand(soilMoisture.select('soil_moisture_am'));
  // Perform linear regression
  var regression = soilMoistureWithTime.reduce(ee.Reducer.linearFit());
  var slope = regression.select('scale').rename('sm_surface_Slope');
  var intercept = regression.select('offset').rename('sm_surface_Intercept');
  return slope.addBands(intercept);
};
var smTrend = soilMoistureTrend();
Map.addLayer(smTrend.select('sm_surface_Slope'), {min: -0.01, max: 0.01, palette: ['blue', 'white', 'red']}, 'Soil Moisture Slope');

// Export the result
Export.image.toAsset({
  image: smTrend,
  description: 'SM_Trend2015_2023',
  assetId:'path_to_export'+'/file_name',
  region:pa,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});





// Load ERA5 dataset
var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")
  .filterDate('1980-01-01', '2023-02-28').filterBounds(pa)
  .map(function(image){
    var x = image.clip(pa)
    return x
  });  // Adjust the date range as necessary

// Get the required bands (dew point temperature and air temperature)
var dewPoint = (era5.select('dewpoint_temperature_2m')).mean();
var temperature = (era5.select('temperature_2m')).mean();
print(temperature,'temp')
// Function to calculate vapor pressure from temperature
var calcVaporPressure = function(temp) {
  //temp = ee.Image(temp)
  var denom = temp.add(243.5)
  var num = temp.multiply(19.67)
  // Convert temperature to vapor pressure using Clausius-Clapeyron equation
  var exponent = (num.divide(denom)).exp();
  return ee.Image(exponent.multiply(6.112))//.rename('relHumidity');  // e(T)
};

// Recalculate RH for each time step in the ERA5 dataset
var rhCollection = era5.map(function(image) {
  // Extract temperature and dew point
  var dewPoint = image.select('dewpoint_temperature_2m');
  var temperature = image.select('temperature_2m');

  // Calculate vapor pressures
  var eT = calcVaporPressure(temperature);
  var eTd = calcVaporPressure(dewPoint);

  // Compute relative humidity
  var RH = eTd.divide(eT).multiply(100).rename('RH');

  // Add time property
  return RH.set('system:time_start', image.get('system:time_start'));
});



// Add a time band to RH collection
var rhWithTime = rhCollection.map(function(image) {
  var time = ee.Date(image.get('system:time_start')).difference(ee.Date('1980-01-01'), 'year');
  return image.addBands(ee.Image.constant(time).rename('time').float());
});


// Perform linear regression to calculate trends
var rhRegression = rhWithTime.reduce(ee.Reducer.linearFit());
var rhSlope = rhRegression.select('scale').rename('rh_Slope');
var rhIntercept = rhRegression.select('offset').rename('rh_Intercept');

// Visualize the slope
Map.addLayer(rhSlope, {min: -1, max: 1, palette: ['blue', 'white', 'red']}, 'RH Slope');

// Export the result
Export.image.toAsset({
  image: rhSlope.addBands(rhIntercept),
  description: 'RH_Trend1980_2023',
  assetId: 'path_to_export'+'file_name',
  region: pa,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

//Trend calculation for NICFI --------------------------------------------------------------------
var nicfi = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/asia')
            .map(function(image){
            return image.clip(pa)
});

// Filter basemaps by date and get the first image from filtered results
var basemap= nicfi.filter(ee.Filter.date('2016-03-01','2023-02-28'))
              .map(function(image){
                var n = image.normalizedDifference(['N','R']).rename('NDVI')
                return image.addBands(n)
              })


var vis = {'bands':['R','G','B'],'min':64,'max':5454,'gamma':1.8};

//Map.addLayer(basemap, vis, '2021-03 mosaic');
//Map.addLayer(basemap.normalizedDifference(['N','R']).rename('NDVI'),{min:-0.55,max:0.8,palette: 
//['8bc4f9', 'c9995c', 'c7d270','8add60','097210'
//    ]}, 'NDVI', false);
    
    
    

// Add a 'time' band to the collection for trend calculation
var planetTime = basemap.map(function(image) {
  var year = ee.Date(image.get('system:time_start')).difference(ee.Date('2016-01-01'),'year');
  return image.addBands(ee.Image.constant(year).rename('time').float());
});

var indices = ['R','G','B','N','NDVI'];


// Function to calculate trends for a single index
var calculateTrend = function(index) {
  // Select the index and time bands
  var stacked = planetTime.select(['time', index]);
  // Perform linear regression
  var regression = stacked.reduce(ee.Reducer.linearFit());
  // Extract slope and intercept
  var slope = regression.select('scale').rename(index + '_Slope'); // Trend slope
  var intercept = regression.select('offset').rename(index + '_Intercept'); // Baseline
  return slope.addBands(intercept);
};

// Calculate trends for all indices and combine results
var planetResults = indices.map(calculateTrend);
print(planetResults,'trends')
planetResults = ee.List(planetResults)
var red = planetResults.get(0)
var green = planetResults.get(1)
var blue= planetResults.get(2)
var nir = planetResults.get(3)
var NDVI = planetResults.get(4)




//Trend Calculation for MODIS LST ------------------------------------------------------------------
// Load MODIS Daily LST (MOD11A1 - Terra)
var modis = ee.ImageCollection('MODIS/061/MOD11A1')
                  .filter(ee.Filter.date('2013-01-01', '2023-02-28'))
                  .filterBounds(pa)
                  .select(['LST_Day_1km', 'LST_Night_1km', 'QC_Day']);

// Function to apply cloud masking based on QC_Day band
var cloudMaskLST = function(image) {
  var qc = image.select('QC_Day');  // Quality control band
  var mask = qc.bitwiseAnd(1).eq(0); // Keep only pixels with best quality
  
  var lstDay = image.select('LST_Day_1km').updateMask(mask).multiply(0.02).subtract(273.15).rename('LST_Day');
  var lstNight = image.select('LST_Night_1km').updateMask(mask).multiply(0.02).subtract(273.15).rename('LST_Night');
  
  return image.addBands(lstDay).addBands(lstNight).select(['LST_Day', 'LST_Night']).clip(pa);
};

// Apply cloud masking and unit conversion
var modisMasked = modis.map(cloudMaskLST);
print(modisMasked, 'Cloud Masked MODIS LST Collection');

// Add a 'time' band for trend analysis
var modisTime = modisMasked.map(function(image) {
  var year = ee.Date(image.get('system:time_start')).difference(ee.Date('2013-01-01'), 'year');
  return image.addBands(ee.Image.constant(year).rename('time').float());
});

// Function to compute trend
var lstTrend = function(index) {
  var stacked = modisTime.select(['time', index]); // Select LST and time bands
  var regression = stacked.reduce(ee.Reducer.linearFit()); // Perform linear regression
  var slope = regression.select('scale').rename(index + '_Slope'); // Trend slope
  var intercept = regression.select('offset').rename(index + '_Intercept'); // Baseline
  return slope.addBands(intercept);
};

// Compute trends for LST_Day and LST_Night
var lstDayTrend = lstTrend('LST_Day');
var lstNightTrend = lstTrend('LST_Night');

// Export the result
Export.image.toAsset({
  image: lstDayTrend,
  description: 'LSTDayMODIS_Trend2013_2023',
  assetId:'path_to_export'+'file_name',
  region: pa,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
// Export the result
Export.image.toAsset({
  image: lstNightTrend,
  description: 'LSTNightMODIS_Trend2013_2023',
  assetId:'path_to_export'+'file_name',
  region: pa,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
// Visualization parameters
var lstVis = { min: -2, max: 2, palette: ['blue', 'white', 'red'] };

// Display results
Map.centerObject(pa, 6);
//Map.addLayer(lstDay.select('LST_Day_Slope'), lstVis, 'LST Daytime Trend');
//Map.addLayer(lstNyt.select('LST_Night_Slope'), lstVis, 'LST Nighttime Trend');



//Import the proximity to roads layer
//var road = ____________________

//Import the DEM slope layer
//var dem_slope = __________________________


//creating input_resampled.tif by combining all the above layers into a single tif.
//Import all the files exported to assets in the above steps


//add all the imported assets and NICFI assets processed locally to "var images" list - an example is given below - modify the names 
var images = [
    rain, rh,stb10, evi,  msavi, mirbi, bsi,sm,ndmi,
  nbr, nbr2,  dem.select('slope'),smi, red,green,blue,nir,NDVI,lstDay,lstNyt,
  road,//waterBuilt15_22,//canopy_ht,roads,ndvi - i have removed landsat ndvi, alosdem.select('b1'), ee.Terrain.slope(alosdem.select('b1'))
];

// Start with the first image
var inputFeat = ee.Image(images[0]);

// Add each subsequent image as bands
for (var i = 1; i < images.length; i++) {
  inputFeat = inputFeat.addBands(ee.Image(images[i]));
}


//Resampling all the layers to 30m 
var proj = inputFeat.select('evi_Slope').projection()
var inputFeatProj = inputFeat.setDefaultProjection(proj)
// Choose the grid size and projection
var gridScale = 30;
var gridProjection = ee.Projection('EPSG:3857')
  .atScale(gridScale);
// Aggregate pixels with 'mean' statistics
var inputResamp = inputFeatProj
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 1024
  })
  .reproject({
    crs: gridProjection
});



Export.image.toAsset({
  image: inputResamp.clip(pa),
  description: 'inputResampled',
  assetId: 'users/sravanthi_mopati/' + 'inputResampled',
  region: pa,
  scale: 30,
  maxPixels: 1e10
});


//Input resampled tif is the output of this code - this tif file has 40 bands.
