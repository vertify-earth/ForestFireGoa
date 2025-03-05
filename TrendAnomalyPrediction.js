//This is the trend anomaly prediction code
//Anomaly prediction - uses the outputs from TrendFire code - 14 layers
//These layers consist of vegetation indices, burn indices and meterological parameters

//true fire events from VIIRS SNPP 375m active fires
var fire13_23 = fire13_19.merge(fire20_23)
var fire13_22 = fire13_23.filter(ee.Filter.stringContains('acq_date', '-2023').not());
var fireMarch2023 = fire13_23.filter(ee.Filter.stringContains('acq_date', '03-2023'))


Map.addLayer(fire13_22,{},'allFirePts',false)
Map.addLayer(fireMarch2023,{},'2023March fire')
Map.centerObject(pa,10)

//Mask to filter cloudy pixels and conversion of DN to surface reflectance
function maskL8sr(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
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
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4'])
    .rename('ndvi');

  var mirbi = image.expression(
    '(10*swir1) - (9.8*swir2) + 2',{
    'swir1':image.select('SR_B6'),
    'swir2':image.select('SR_B7')
      
    }).rename('mirbi')
    

  var ndfi = image.normalizedDifference(['ST_B10', 'SR_B6'])
    .rename('ndfi');

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
    .addBands(ndvi)
    .addBands(mirbi)
    .addBands(ndfi)
    .addBands(evi)
    .addBands(bsi)
    .addBands(ndmi)
    .addBands(nbr)
    .addBands(nbr2)
    .addBands(msavi);
};






//preparing landsat data for time series chart 
// Time series chart generation data is hashed below.
var outlierDate = ee.Date('2016-06-22');
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                     .filterDate('2013-03-20', '2024-02-28').filterBounds(pa)
                     .filterMetadata('CLOUD_COVER', 'less_than', 10)
                     .map(maskL8sr)
                     .map(function(image) {
                    return image.clip(pa);
                  })
                  .map(addIndices);


landsat = landsat
  .filterDate('2013-03-20', '2024-02-28')
  .filterBounds(pa)
  .filter(
    ee.Filter.or(
      ee.Filter.date('2013-03-20', outlierDate),  // Keep images before outlier
      ee.Filter.date(outlierDate.advance(1, 'day'), '2024-02-28') // Keep images after outlier
    )
  );

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

landsat = landsat.select(['ndvi', 'evi', 'mirbi','ndfi','bsi','ndmi','nbr','nbr2','msavi','smi','ST_B10']);
landsat = landsat.sort('system:time_start');

print(landsat,'landsat')

// Function to extract mean values over the polygon
var extractTimeSeries = function(image) {
  var stats = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: pa,  // Use your polygon
    scale: 30,
    bestEffort: true,
    maxPixels: 1e13
  });
  
  // Add the date as a property
  return ee.Feature(null, stats.set('date', image.date().format('YYYY-MM-dd')));
};

// Apply the function to the image collection
var timeSeries = landsat.map(extractTimeSeries);

// Convert to FeatureCollection
var timeSeriesFC = ee.FeatureCollection(timeSeries);



// List of indices to plot separately
var indices = ['ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2', 'msavi', 'smi', 'ST_B10'];


/*
// Loop through each index and create a separate chart - a chart to view the time series of all indices
indices.forEach(function(index) {
  var chart = ui.Chart.feature.byFeature(timeSeriesFC, 'date', [index])
    .setOptions({
      title: index.toUpperCase() + ' Over Time',
      hAxis: {title: 'Date', format: 'YYYY-MM', gridlines: {count: 10}},
      vAxis: {title: index.toUpperCase() + ' Value'},
      lineWidth: 2,
      pointSize: 3
    });

  print(chart);
});
*/


//Trend Anomaly prediction begins here-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//The inputs from TrendFire need to be inserted here - first import them as 
//assets, rename them and put them in the list of indices here
var images = [ smi,stb10, ndvi,evi,  msavi, mirbi, ndmi,ndfi,nbr, nbr2,bsi];

// Start with the first image
var inputFeat = ee.Image(images[0]);

// Add each subsequent image as bands - image = stack of all indices as bands in the image
for (var i = 1; i < images.length; i++) {
  inputFeat = inputFeat.addBands(ee.Image(images[i]));
}


print(inputFeat)






///Hotspot prediction 

// Define March 2023 time in years since 2013-03-20
var const2023 = ee.Date('2023-03-01').difference(ee.Date('2013-03-20'), 'year');
print(const2023)
print(inputFeat.bandNames())
// Select all bands ending with "_Slope" using a proper regex pattern

var slope = inputFeat.bandNames().filter(ee.Filter.stringContains('item', 'Slope'));
var slopeBands = inputFeat.select(slope);
var inter = inputFeat.bandNames().filter(ee.Filter.stringContains('item', 'Intercept'));
var interceptBands = inputFeat.select(inter);

print(slopeBands,'slopeBands')
print(interceptBands,'interceptBands')

// Compute predicted values for all indices using the regression expression
var predictedIndices = slopeBands.multiply(const2023).add(interceptBands)
                        .rename('smi','ST_B10','ndvi','evi','msavi','mirbi'
                        ,'ndmi','ndfi','nbr','nbr2','bsi');

print(predictedIndices,'predIndices')
//the number of years here more than 10 years - 1982 to 2021 ~ 42.4 
var constRain = ee.Date('2022-01-01').difference(ee.Date('1982-01-01'), 'year');
var predRain = (rain.select('rain_Slope').multiply(constRain).add(rain.select('rain_Intercept'))).rename('rain')
//print(predRain,'predRain')
//Map.addLayer(predRain,{},'Pred rain')

var constRh = ee.Date('2023-03-01').difference(ee.Date('1980-01-01'), 'year');
var predRh = (rh.select('rh_Slope').multiply(constRain).add(rh.select('rh_Intercept'))).rename('rh')

var constSm = ee.Date('2023-03-01').difference(ee.Date('2015-04-01'), 'year');
var predSm = (sm.select('sm_surface_Slope').multiply(constRain).add(sm.select('sm_surface_Intercept'))).rename('sm')

predictedIndices = predictedIndices.addBands([predRain,predRh,predSm])

print(predictedIndices,'predictedIndices')


var presentls = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                .filterDate('2023-02-01', '2023-02-28').filterBounds(pa)
                .map(maskL8sr).map(addIndices)
                .map(function(image){
                  return image.clip(pa)
                });



// Function to calculate SMI for the collection
var addSMI = function(image) {
  // Access globally defined min/max values
  var smi = image.select('SR_B7').subtract(swirMin)
                 .divide(swirMax.subtract(swirMin))
                 .multiply(-1).add(1).rename('smi');
  return image.addBands(smi);
};
// Calculate SWIR min/max for the entire collection
var swirStats = (presentls.mean().select('SR_B7')).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});

var swirMin = ee.Number(swirStats.get('SR_B7_min'));
var swirMax = ee.Number(swirStats.get('SR_B7_max'));

presentls = presentls.map(addSMI)
presentls = presentls.mosaic()
Map.addLayer(presentls,{},'March 2023',false)




// Load ERA5 dataset
var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")
  .filterDate('2023-02-01', '2023-02-28').filterBounds(pa)
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
  var RH = eTd.divide(eT).multiply(100).rename('rh');

  // Add time property
  return RH.set('system:time_start', image.get('system:time_start'));
});
print(rhCollection)
var present_rh = rhCollection.mosaic()


// Visualize the relative humidity
//Map.addLayer(present_rh, {min: 0, max: 100, palette: ['blue', 'green', 'yellow', 'red']}, 'Relative Humidity',false);




var soilMoisture = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005').filterBounds(pa)
                        .map(function(image){
                        var x = image.clip(pa)
                        return x;
                        }).filterDate('2023-02-01', '2023-02-28')
                            .select(['soil_moisture_am']);
var present_sm = soilMoisture.mosaic().rename('sm')
// Define visualization parameters.
var soilMoistureVis = {
  min: 0.0,
  max: 0.7,
  palette: ['A67C00', 'FFE625', 'C2E5D3', '90DCD0',
            '2FBDBD', '0C9BBD', '068682'],
};
//Map.addLayer(present_sm,soilMoistureVis,'SM_present')

//var rain2022 =  



///Precipitation trend analysis
// Load CHIRPS daily data for the time period 1982-2024
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
              .filter(ee.Filter.date('2022-01-01', '2022-12-31'))
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

// Create a list of yearly means for each year from 2021 to 2022
var years = ee.List.sequence(2021,2022);  // Define the range of years
var yearlySums = years.map(createYearlySum);

var rain2022 = (yearlySums.filter(ee.Filter.eq('year',2022)))
rain2022 = ee.Image(rain2022.get(0)).rename('rain')
print(rain2022,'rain2022')
var prepVis = {bands: ['precipitation'],palette:['#d3ffd2', '#99FF33', 
'#FFFF00', '#FF9900', '#FF3300', '#FF0000'],
 min: 1000, max: 5000}

// Convert the list of yearly means into an ImageCollection
var yearlySumsCollection = ee.ImageCollection(yearlySums);
print(yearlySumsCollection)








var present_all = presentls.addBands([rain2022,present_sm,present_rh.select('rh')])//merge all presentls+present_smi+present_rh+present_sm
print(present_all,'presentAll')
var past_Pred = predictedIndices.select('smi','ST_B10','ndvi','evi','msavi','mirbi'
                        ,'ndmi','ndfi','nbr','nbr2','bsi','rain','rh','sm')

print(past_Pred,'past predicted')
                      
var anomalyImage = present_all.select('smi','ST_B10','ndvi','evi','msavi','mirbi'
                        ,'ndmi','ndfi','nbr','nbr2','bsi','rain','rh','sm').subtract(past_Pred)
  .rename(past_Pred.bandNames().map(function(b) 
  {
    return ee.String(b).cat('_Anomaly');
  }));

print(anomalyImage,'anomalyImage')

//Rain Anomaly ---------------------------------------------------------------------------------
// Define fire hotspot conditions
var rainAnomaly = anomalyImage.select('rain_Anomaly')
var rainhs = (rainAnomaly.lt(0)).selfMask();//.or(rainAnomaly.lt(50))
Map.addLayer(rainhs,{},'rain Hotspot')


///RH Anamoly ----------------------------------------------------------------------------------
// Calculate RH min/max for the entire collection
var rhAnomaly = anomalyImage.select('rh_Anomaly')
var rhMinMax = (rhAnomaly).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});
print(rhMinMax,'rhMinMax')
var rhhs = (rhAnomaly.lt(60)).selfMask()
Map.addLayer(rhhs,{},'rh hotspot')




///Soil Moisture Anamoly------------------------------------------------------------------------
// Calculate soil Moisture min/max for the entire collection
var smAnomaly = anomalyImage.select('sm_Anomaly')
var smMinMax = (smAnomaly).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});
print(smMinMax,'smMinMax')
var smhs = (smAnomaly.lt(-50)).selfMask()
//Map.addLayer(smhs,{},'sm hotspot')

//Landsat indices anomaly ----------------------------------------------------------------------
var vegIndices = anomalyImage.select('ndvi_Anomaly','evi_Anomaly','msavi_Anomaly','ndmi_Anomaly')
var vegMinMax = (vegIndices).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});
print(vegMinMax,'vegMinMax')
var vegIndices_hs = vegIndices.lt(-0.15)
Map.addLayer(vegIndices_hs,{opacity:0.5},'vegetation indices Hotspots')
//Map.addLayer(vegIndices,{},'vegIndices')

//Landsat burn indices
var burnIndices = anomalyImage.select('mirbi_Anomaly','ST_B10_Anomaly','bsi_Anomaly')
var burnMinMax = (burnIndices).reduceRegion({
  reducer: ee.Reducer.minMax(),
  //selectors: ['SR_B7'],
  geometry: pa,  // Use AOI
  scale: 30  
});
print(burnMinMax,'burnMinMax')
Map.addLayer(anomalyImage.select('ST_B10_Anomaly').gt(3.5),{},'b10Anomaly')


var burnIndices_hs = burnIndices.gt(3.2)//.selfMask()
Map.addLayer(burnIndices_hs,{opacity:0.5},'Burn indices Hotspots')

Export.image.toDrive({
  image: burnIndices.select('ST_B10_Anomaly').gt(3.2),
  description: 'fireHotspotSTB10',
  scale:30,
  //crs: projection.crs,
  //crsTransform: projection.transform,
  region: pa
});



Export.image.toDrive({
  image: vegIndices.select('ndmi_Anomaly').lt(-0.15),
  description: 'fireHotspotNDMI',
  scale:30,
  //crs: projection.crs,
  //crsTransform: projection.transform,
  region: pa
});





//Using LULC to show interc=section between detected hotspots and LULC 


var colFilter = ee.Filter.and(
    ee.Filter.bounds(pa),
    ee.Filter.date('2023-02-01', '2023-02-28'));

var dwCol = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filter(colFilter);
print(dwCol)

dwCol = dwCol.map(function(image){
  var n = image.clip(pa)
  return n
})

Map.centerObject(dwCol,10)

var CLASS_NAMES = [
    'water', 'trees', 'grass', 'flooded_vegetation', 'crops',
    'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];
var VIS_PALETTE = [
    '419bdf', '397d49', '88b053', '7a87c6', 'e49635', 'dfc35a', 'c4281b',
    'a59b8f', 'b39fe1'];

var mos = dwCol.mosaic()
mos = ee.Image(mos)
//Map.addLayer(mos,{min:1,max:1,palette: VIS_PALETTE},'Mosaic')
print(mos,'mosaic')


// Filter for the 'crops' class
var cropsOnly = mos.updateMask(mos.select('label').eq(4)).select('label');

// Visualize the 'crops' class - the pallete is irrelevant here
Map.addLayer(cropsOnly, {min: 0, max: 0, }, 'Crops Only');

var dwRgb = mos
    .select('label')
    .visualize({min: 0, max: 8, palette: VIS_PALETTE})
    .divide(255);
print(dwRgb,'dwrgb')
Map.addLayer(dwRgb,{},'dwrgb')



//Create buffer around the crops layer:
//Then check the intersection between Crops and firehotspots.

// Convert the 'crops' mask to a vector feature collection
var cropVectors = cropsOnly.reduceToVectors({
  geometry:pa,
  geometryType: 'polygon',
  reducer: ee.Reducer.countEvery(),
  scale: 30, // Adjust based on your data resolution
  maxPixels: 1e8
});

// Create a 1 km buffer around the crop features
var cropBuffer = cropVectors.map(function(feature) {
  return (feature.bounds()).buffer(100); // 1 km buffer
});

// Visualize the buffered area
Map.addLayer(cropBuffer, {color: 'blue'}, '1 km crop Buffer');

// Mask the thresholded image using the buffer
var interBurn = burnIndices.clip(cropBuffer);
var interVeg = vegIndices.clip(cropBuffer)
// Visualize the intersection
Map.addLayer(interBurn, {}, 'Burn Hotspots intersecting buffer');
Map.addLayer(interVeg, {}, 'Veg Hotspots intersecting buffer');









