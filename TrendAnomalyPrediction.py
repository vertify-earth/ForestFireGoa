#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TrendAnomalyPrediction.py: Python implementation of TrendAnomalyPrediction.js using Earth Engine Python API.
This script calculates anomalies between current conditions and predicted values from historical trends
to identify potential fire hotspots in Goa, India.

The script performs the following steps:
1. Processes current Landsat imagery to calculate various spectral indices
2. Uses trend coefficients to predict expected values based on historical data
3. Calculates anomalies between current and predicted values
4. Identifies hotspots based on anomalies in rainfall, RH, soil moisture, and vegetation indices
5. Integrates land use/land cover data to identify vulnerable areas
6. Exports hotspot maps to Google Drive and/or Earth Engine assets

Equivalent to the TrendAnomalyPrediction.js script in the original implementation.
"""

import ee
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

# Initialize the Earth Engine API
try:
    ee.Initialize()
    print("Earth Engine API initialized successfully.")
except Exception as e:
    print(f"Error initializing Earth Engine API: {e}")
    print("Make sure you have authenticated with Earth Engine using ee.Authenticate()")

# Define the Goa study area boundary
def get_study_boundary():
    """
    Get the study area boundary for Goa, India.
    Returns an ee.Geometry.
    """
    # Method 1: Load from a shapefile if available locally
    try:
        boundary_path = os.path.join('data', 'pa_boundary.shp')
        if os.path.exists(boundary_path):
            boundary_gdf = gpd.read_file(boundary_path)
            # Convert to GeoJSON
            boundary_geojson = boundary_gdf.geometry.__geo_interface__
            # Create an EE geometry
            return ee.Geometry.Polygon(boundary_geojson['features'][0]['geometry']['coordinates'])
    except Exception as e:
        print(f"Could not load boundary from shapefile: {e}")
    
    # Method 2: Use predefined feature from Earth Engine
    try:
        # Get the Goa district boundary from Earth Engine's administrative boundaries
        goa = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(ee.Filter.eq('ADM1_NAME', 'Goa'))
        return goa.geometry()
    except Exception as e:
        print(f"Could not load boundary from Earth Engine: {e}")
        
    # Method 3: Define manually as fallback
    # These coordinates would need to be adjusted to match the study area
    goa_coords = [
        [73.6765, 15.7560],
        [74.3161, 15.7560],
        [74.3161, 14.8922],
        [73.6765, 14.8922],
        [73.6765, 15.7560]
    ]
    return ee.Geometry.Polygon(goa_coords)

# Function to load fire events data
def load_fire_events(pa):
    """
    Load fire events data from 2013-2023 and merge into a single feature collection.
    
    Args:
        pa: ee.Geometry representing the study area
    
    Returns:
        ee.FeatureCollection of fire events
    """
    print("Loading fire events data...")
    
    # Load fire events data (2013-2019)
    try:
        # If the data is available in Earth Engine, you can load it directly
        fire_13_19 = ee.FeatureCollection("users/your_username/fire13_19")
    except:
        # Otherwise, try to find and upload it, or use a placeholder
        print("Fire data 2013-2019 not available in Earth Engine. Using placeholder.")
        fire_13_19 = ee.FeatureCollection([])
    
    # Load fire events data (2020-2023)
    try:
        # If the data is available in Earth Engine, you can load it directly
        fire_20_23 = ee.FeatureCollection("users/your_username/fire20_23")
    except:
        # Otherwise, try to find and upload it, or use a placeholder
        print("Fire data 2020-2023 not available in Earth Engine. Using placeholder.")
        fire_20_23 = ee.FeatureCollection([])
    
    # Merge the two collections
    fire_all = fire_13_19.merge(fire_20_23)
    
    # Filter by the study area boundary
    fire_all = fire_all.filterBounds(pa)
    
    # Add a 'year' property based on the acquisition date
    def add_year(feature):
        acq_date = ee.String(feature.get('ACQ_DATE'))
        year = ee.Number.parse(acq_date.slice(0, 4))
        return feature.set('year', year)
    
    fire_all = fire_all.map(add_year)
    
    # Count the number of fire points
    fire_count = fire_all.size().getInfo()
    print(f"Number of fire events in the study area: {fire_count}")
    
    return fire_all

# Function to mask clouds in Landsat 8/9 imagery
def maskL8sr(image):
    """
    Mask clouds and scale pixel values for Landsat 8/9 SR imagery.
    """
    # Get QA band
    qa = image.select('QA_PIXEL')
    
    # Bits 3 and 5 are cloud shadow and cloud, respectively
    cloud_shadow_bit_mask = 1 << 3
    clouds_bit_mask = 1 << 5
    
    # Both flags should be set to zero, indicating clear conditions
    mask = qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0).And(
           qa.bitwiseAnd(clouds_bit_mask).eq(0))
    
    # Scale the optical bands
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    
    # Scale the thermal bands
    thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    
    # Return the masked and scaled image
    return image.select(['QA_.*']).addBands(optical_bands).addBands(thermal_bands).updateMask(mask)

# Function to calculate spectral indices and add them as bands
def addIndices(image):
    """
    Calculate various spectral indices and add them as bands to the image.
    """
    # NDVI (Normalized Difference Vegetation Index)
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('ndvi')
    
    # EVI (Enhanced Vegetation Index)
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select('SR_B5'),
            'RED': image.select('SR_B4'),
            'BLUE': image.select('SR_B2')
        }
    ).rename('evi')
    
    # MSAVI (Modified Soil Adjusted Vegetation Index)
    msavi = image.expression(
        '(2 * NIR + 1 - sqrt((2 * NIR + 1) * (2 * NIR + 1) - 8 * (NIR - RED))) / 2',
        {
            'NIR': image.select('SR_B5'),
            'RED': image.select('SR_B4')
        }
    ).rename('msavi')
    
    # MIRBI (Mid-Infrared Burn Index)
    mirbi = image.expression(
        '10 * SWIR2 - 9.8 * SWIR1 + 2',
        {
            'SWIR1': image.select('SR_B6'),
            'SWIR2': image.select('SR_B7')
        }
    ).rename('mirbi')
    
    # NDMI (Normalized Difference Moisture Index)
    ndmi = image.normalizedDifference(['SR_B5', 'SR_B6']).rename('ndmi')
    
    # NDFI (Normalized Difference Fraction Index)
    ndfi = image.normalizedDifference(['SR_B7', 'SR_B5']).rename('ndfi')
    
    # NBR (Normalized Burn Ratio)
    nbr = image.normalizedDifference(['SR_B5', 'SR_B7']).rename('nbr')
    
    # NBR2 (Normalized Burn Ratio 2)
    nbr2 = image.normalizedDifference(['SR_B6', 'SR_B7']).rename('nbr2')
    
    # BSI (Bare Soil Index)
    bsi = image.expression(
        '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))',
        {
            'SWIR1': image.select('SR_B6'),
            'RED': image.select('SR_B4'),
            'NIR': image.select('SR_B5'),
            'BLUE': image.select('SR_B2')
        }
    ).rename('bsi')
    
    # Add all indices as bands
    return image.addBands([ndvi, evi, msavi, mirbi, ndmi, ndfi, nbr, nbr2, bsi])

# Function to calculate SMI (Soil Moisture Index)
def calculateSMI(image):
    """
    Calculate Soil Moisture Index (SMI) from Landsat data.
    """
    Ts = image.select('ST_B10')
    ndvi = image.select('ndvi')
    
    # Apply the SMI calculation
    smi = image.expression(
        '(Ts_max - Ts) / (Ts_max - Ts_min)',
        {
            'Ts': Ts,
            'Ts_max': Ts.reduceRegion(
                reducer=ee.Reducer.max(),
                geometry=image.geometry(),
                scale=30,
                maxPixels=1e9
            ).get('ST_B10'),
            'Ts_min': Ts.reduceRegion(
                reducer=ee.Reducer.min(),
                geometry=image.geometry(),
                scale=30,
                maxPixels=1e9
            ).get('ST_B10')
        }
    ).rename('smi')
    
    return image.addBands(smi)

# Function to get current Landsat data
def get_current_landsat(pa, date_range):
    """
    Get the most recent Landsat 8/9 imagery for a specified date range.
    
    Args:
        pa: ee.Geometry representing the study area
        date_range: Dictionary with 'start' and 'end' dates
    
    Returns:
        ee.Image with all spectral indices and bands
    """
    print(f"Getting current Landsat data for {date_range['start']} to {date_range['end']}...")
    
    # Load Landsat collection and filter by date and location
    ls_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')) \
        .filterDate(date_range['start'], date_range['end']) \
        .filterBounds(pa) \
        .filter(ee.Filter.lt('CLOUD_COVER', 20))
    
    # Apply cloud masking
    ls_masked = ls_collection.map(maskL8sr)
    
    # Add indices to all images in collection
    ls_indices = ls_masked.map(addIndices)
    
    # Calculate SMI for each image
    ls_with_smi = ls_indices.map(calculateSMI)
    
    # Create mosaic of the collection
    current_ls = ls_with_smi.mosaic().clip(pa)
    
    return current_ls

# Function to calculate vapor pressure from temperature
def calcVaporPressure(temp):
    """
    Calculate vapor pressure from temperature using the Clausius-Clapeyron equation.
    
    Args:
        temp: ee.Image with temperature in Celsius
    
    Returns:
        ee.Image with vapor pressure (hPa)
    """
    denom = temp.add(243.5)
    num = temp.multiply(17.67)
    exponent = num.divide(denom).exp()
    return exponent.multiply(6.112)

# Function to get current ERA5 relative humidity data
def get_current_relative_humidity(pa, date_range):
    """
    Get the most recent ERA5 relative humidity data.
    
    Args:
        pa: ee.Geometry representing the study area
        date_range: Dictionary with 'start' and 'end' dates
    
    Returns:
        ee.Image with relative humidity
    """
    print(f"Getting current relative humidity data for {date_range['start']} to {date_range['end']}...")
    
    # Load ERA5 data
    era5 = ee.ImageCollection('ECMWF/ERA5/DAILY') \
        .filterDate(date_range['start'], date_range['end']) \
        .filterBounds(pa) \
        .select(['dewpoint_temperature_2m', 'temperature_2m'])
    
    # Calculate relative humidity
    def calculateRH(image):
        dewPoint = image.select('dewpoint_temperature_2m')
        temperature = image.select('temperature_2m')
        
        # Calculate vapor pressures
        eT = calcVaporPressure(temperature)
        eTd = calcVaporPressure(dewPoint)
        
        # Compute relative humidity (%)
        RH = eTd.divide(eT).multiply(100).rename('rh')
        
        # Add time property
        return RH.set('system:time_start', image.get('system:time_start'))
    
    # Apply RH calculation to the collection
    rh_collection = era5.map(calculateRH)
    
    # Create mosaic of the collection
    current_rh = rh_collection.mosaic().clip(pa)
    
    return current_rh

# Function to get current SMAP soil moisture data
def get_current_soil_moisture(pa, date_range):
    """
    Get the most recent SMAP soil moisture data.
    
    Args:
        pa: ee.Geometry representing the study area
        date_range: Dictionary with 'start' and 'end' dates
    
    Returns:
        ee.Image with soil moisture
    """
    print(f"Getting current soil moisture data for {date_range['start']} to {date_range['end']}...")
    
    # Load SMAP data
    smap = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005') \
        .filterDate(date_range['start'], date_range['end']) \
        .filterBounds(pa) \
        .select(['soil_moisture_am'])
    
    # Create mosaic of the collection
    current_sm = smap.mosaic().clip(pa).rename('sm')
    
    return current_sm

# Function to get current CHIRPS precipitation data
def get_current_precipitation(pa, date_range):
    """
    Get the most recent CHIRPS precipitation data.
    
    Args:
        pa: ee.Geometry representing the study area
        date_range: Dictionary with 'start' and 'end' dates
    
    Returns:
        ee.Image with precipitation
    """
    print(f"Getting current precipitation data for {date_range['start']} to {date_range['end']}...")
    
    # Load CHIRPS data
    chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
        .filterDate(date_range['start'], date_range['end']) \
        .filterBounds(pa)
    
    # Create sum of the collection
    current_rain = chirps.sum().clip(pa).rename('rain')
    
    return current_rain

# Function to load trend coefficients from the TrendFire output
def load_trend_coefficients(pa):
    """
    Load trend coefficients previously calculated using TrendFire.py.
    
    Args:
        pa: ee.Geometry representing the study area
    
    Returns:
        ee.Image with all trend coefficients
    """
    print("Loading trend coefficients...")
    
    try:
        # Try to load the trend layers from Earth Engine assets
        trends = ee.Image("users/your_username/goa_fire_trends").clip(pa)
        print("Successfully loaded trend coefficients from Earth Engine asset.")
        return trends
    except:
        print("Trend coefficients not found in Earth Engine assets. Using placeholder.")
        # Create a placeholder image with zeros for all bands
        # This should be replaced with actual trend calculations
        placeholder = ee.Image(0).rename('ndvi_trend')
        return placeholder

# Function to predict index values based on trend coefficients
def predict_from_trends(trend_coefficients, current_time, reference_time):
    """
    Predict index values for the current time based on historical trends.
    
    Args:
        trend_coefficients: ee.Image with trend coefficients
        current_time: Current time in decimal years since start
        reference_time: Reference time in decimal years for the intercept
    
    Returns:
        ee.Image with predicted values
    """
    print("Predicting index values from trend coefficients...")
    
    # Construct regression expressions
    # Predicted value = intercept + slope * (current_time - reference_time)
    time_diff = ee.Number(current_time).subtract(reference_time)
    
    # Create prediction expressions for each index
    indices = [
        'ndvi', 'evi', 'msavi', 'mirbi', 'ndmi', 'ndfi', 
        'nbr', 'nbr2', 'bsi', 'smi', 'ST_B10', 'rain', 'rh', 'sm'
    ]
    
    predicted_indices = ee.Image.constant(0)
    
    for index in indices:
        # Try to get trend coefficient for this index
        try:
            trend = trend_coefficients.select(f'{index}_trend')
            
            # Calculate predicted value
            # Here we're assuming the intercept is the current average value
            # This is a simplification - in a full implementation we'd have both slope and intercept
            predicted = trend.multiply(time_diff).rename(index)
            
            # Add to the predicted indices image
            predicted_indices = predicted_indices.addBands(predicted)
        except:
            print(f"Warning: Could not find trend coefficient for {index}")
    
    return predicted_indices

# Function to calculate anomalies between current and predicted values
def calculate_anomalies(current, predicted):
    """
    Calculate anomalies between current observations and predicted values.
    
    Args:
        current: ee.Image with current observations
        predicted: ee.Image with predicted values
    
    Returns:
        ee.Image with anomalies for each band
    """
    print("Calculating anomalies...")
    
    # Calculate anomalies (current - predicted)
    anomaly_image = current.subtract(predicted)
    
    # Rename bands to indicate they are anomalies
    band_names = anomaly_image.bandNames()
    new_names = band_names.map(lambda name: ee.String(name).cat('_Anomaly'))
    
    return anomaly_image.rename(new_names)

# Function to identify hotspots based on anomalies
def identify_hotspots(anomaly_image, thresholds):
    """
    Identify potential fire hotspots based on anomalies.
    
    Args:
        anomaly_image: ee.Image with anomalies for each band
        thresholds: Dictionary with threshold values for each anomaly
    
    Returns:
        Dictionary with hotspot layers for different indices
    """
    print("Identifying potential fire hotspots...")
    
    # Rainfall hotspot (deficit)
    rain_anomaly = anomaly_image.select('rain_Anomaly')
    rain_hotspot = rain_anomaly.lt(thresholds['rain']).selfMask().rename('rain_hotspot')
    
    # Relative humidity hotspot (deficit)
    rh_anomaly = anomaly_image.select('rh_Anomaly')
    rh_hotspot = rh_anomaly.lt(thresholds['rh']).selfMask().rename('rh_hotspot')
    
    # Soil moisture hotspot (deficit)
    sm_anomaly = anomaly_image.select('sm_Anomaly')
    sm_hotspot = sm_anomaly.lt(thresholds['sm']).selfMask().rename('sm_hotspot')
    
    # Vegetation indices hotspots (deficit indicating stress)
    veg_indices = anomaly_image.select(['ndvi_Anomaly', 'evi_Anomaly', 'msavi_Anomaly', 'ndmi_Anomaly'])
    veg_hotspot = veg_indices.lt(thresholds['veg']).selfMask().rename('veg_hotspot')
    
    # Burn indices hotspots (excess indicating fire potential)
    burn_indices = anomaly_image.select(['mirbi_Anomaly', 'ST_B10_Anomaly', 'bsi_Anomaly'])
    burn_hotspot = burn_indices.gt(thresholds['burn']).selfMask().rename('burn_hotspot')
    
    # Combine all hotspots
    all_hotspots = ee.Image.cat([
        rain_hotspot, rh_hotspot, sm_hotspot, veg_hotspot, burn_hotspot
    ])
    
    return {
        'rain': rain_hotspot,
        'rh': rh_hotspot,
        'sm': sm_hotspot,
        'veg': veg_hotspot,
        'burn': burn_hotspot,
        'all': all_hotspots
    }

# Function to get land use/land cover data
def get_land_cover(pa, date_range):
    """
    Get land use/land cover data for the study area.
    
    Args:
        pa: ee.Geometry representing the study area
        date_range: Dictionary with 'start' and 'end' dates
    
    Returns:
        ee.Image with land cover classification
    """
    print("Getting land use/land cover data...")
    
    # Load Dynamic World land cover data
    dw_collection = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1') \
        .filterDate(date_range['start'], date_range['end']) \
        .filterBounds(pa)
    
    # Create a cloud-masked composite
    dw_composite = dw_collection.map(lambda img: img.clip(pa)).mosaic()
    
    # Get the label band (classification)
    lulc = dw_composite.select('label')
    
    return lulc

# Function to analyze hotspots by land cover type
def analyze_hotspots_by_land_cover(hotspots, lulc, pa):
    """
    Analyze hotspots by land cover type.
    
    Args:
        hotspots: Dictionary with hotspot layers
        lulc: ee.Image with land cover classification
        pa: ee.Geometry representing the study area
    
    Returns:
        Dictionary with hotspot statistics by land cover type
    """
    print("Analyzing hotspots by land cover type...")
    
    # Define land cover classes
    land_cover_classes = {
        0: 'water',
        1: 'trees',
        2: 'grass',
        3: 'flooded_vegetation',
        4: 'crops',
        5: 'shrub_and_scrub',
        6: 'built',
        7: 'bare',
        8: 'snow_and_ice'
    }
    
    # Analyze each hotspot type by land cover
    results = {}
    
    for hotspot_name, hotspot_image in hotspots.items():
        if hotspot_name == 'all':
            continue
            
        # Create a binary mask of the hotspot
        hotspot_mask = hotspot_image.gt(0)
        
        # Calculate area of hotspot by land cover class
        hotspot_by_lulc = hotspot_mask.multiply(ee.Image.pixelArea()).addBands(lulc) \
            .reduceRegion(
                reducer=ee.Reducer.sum().group(1),
                geometry=pa,
                scale=30,
                maxPixels=1e13
            )
        
        # Extract results
        groups = ee.List(hotspot_by_lulc.get('groups'))
        
        # Store results
        results[hotspot_name] = groups
    
    return results

# Function to export hotspot maps to Google Drive
def export_hotspots(hotspots, pa, output_folder='ForestFireGoa'):
    """
    Export hotspot maps to Google Drive.
    
    Args:
        hotspots: Dictionary with hotspot layers
        pa: ee.Geometry representing the study area
        output_folder: Google Drive folder name
    """
    print(f"Exporting hotspot maps to Google Drive folder '{output_folder}'...")
    
    # Export each hotspot type
    for hotspot_name, hotspot_image in hotspots.items():
        task = ee.batch.Export.image.toDrive({
            'image': hotspot_image,
            'description': f'fire_hotspot_{hotspot_name}',
            'folder': output_folder,
            'fileNamePrefix': f'fire_hotspot_{hotspot_name}',
            'region': pa,
            'scale': 30,
            'maxPixels': 1e13
        })
        
        task.start()
        print(f"Started export task for {hotspot_name} hotspot with ID: {task.id}")

# Main function to run the entire workflow
def main():
    """
    Main function to run the trend anomaly prediction workflow.
    """
    print("Starting forest fire trend anomaly prediction for Goa, India...")
    
    # Get the study area boundary
    pa = get_study_boundary()
    
    # Define date range for current conditions (e.g., latest month)
    # These dates should be adjusted based on your specific analysis period
    current_date_range = {
        'start': '2023-02-01',
        'end': '2023-02-28'
    }
    
    # Define reference time for trend calculations (e.g., start of historical period)
    reference_time = 2013.0  # Start of analysis in decimal years
    
    # Calculate current time in decimal years
    current_year = 2023
    current_month = 2
    current_time = current_year + (current_month - 1) / 12
    
    # Load fire events data
    fire_events = load_fire_events(pa)
    
    # Get current observational data
    current_ls = get_current_landsat(pa, current_date_range)
    current_rh = get_current_relative_humidity(pa, current_date_range)
    current_sm = get_current_soil_moisture(pa, current_date_range)
    current_rain = get_current_precipitation(pa, current_date_range)
    
    # Combine all current data
    current_all = current_ls.addBands([current_rain, current_sm, current_rh])
    
    # Load trend coefficients from TrendFire output
    trend_coefficients = load_trend_coefficients(pa)
    
    # Predict values based on trend coefficients
    predicted_values = predict_from_trends(trend_coefficients, current_time, reference_time)
    
    # Calculate anomalies
    anomaly_image = calculate_anomalies(current_all, predicted_values)
    
    # Define thresholds for hotspot identification
    # These thresholds would need to be calibrated for your specific study area
    thresholds = {
        'rain': 0,       # Rainfall below predicted value
        'rh': 60,        # RH below 60% of predicted value
        'sm': -50,       # Soil moisture 50% below predicted value
        'veg': -0.15,    # Vegetation indices 0.15 below predicted value
        'burn': 3.2      # Burn indices 3.2 above predicted value
    }
    
    # Identify hotspots
    hotspots = identify_hotspots(anomaly_image, thresholds)
    
    # Get land cover data
    lulc = get_land_cover(pa, current_date_range)
    
    # Analyze hotspots by land cover
    hotspot_analysis = analyze_hotspots_by_land_cover(hotspots, lulc, pa)
    
    # Export hotspot maps to Google Drive (uncomment to use)
    # export_hotspots(hotspots, pa)
    
    print("Forest fire trend anomaly prediction complete!")
    print("To visualize the results, check your Google Drive for the exported hotspot maps.")

if __name__ == "__main__":
    main() 