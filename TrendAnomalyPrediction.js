//This is the trend anomaly prediction code
//Anomaly prediction - uses the outputs from TrendFire code - inputResampled.tif with 40 bands
//These layers consist of vegetation indices, burn indices and meterological parameters

//true fire events from VIIRS SNPP 375m active fires
var fire13_23 = fire13_19.merge(fire20_23)
var fire13_22 = fire13_23.filter(ee.Filter.stringContains('acq_date', '-2023').not());
var fireMarch2023 = fire13_23.filter(ee.Filter.stringContains('acq_date', '03-2023'))


Map.addLayer(fire13_22,{},'allFirePts',false)
Map.addLayer(fireMarch2023,{},'2023March fire')
Map.centerObject(pa,10)



//Trend Anomaly prediction begins here-------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//The inputs from TrendFire need to be inserted here - first import them as 
//assets, rename them and put them in the list of indices here - here we are using landsat indices to predict forest fire prone hotspot areas
//the below code is a sample evaluation for this area of interest.

var inputFeat = inputResampled.select(['smi_Slope','smi_Intercept','ST_B10_Slope','ST_B10_Intercept','evi_Slope','evi_Intercept',
'msavi_Slope','msavi_Intercept','mirbi_Slope','mirbi_Intercept','nbr_Slope','nbr_Intercept','nbr2_Slope','nbr2_Intercept','bsi_Slope','bsi_Intercept'])



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









