#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TrendFire.py: Python implementation of TrendFire.js using Earth Engine Python API.
This script calculates various vegetation, burn, and climate trends for a study area in Goa, India.

The script performs the following steps:
1. Processes Landsat 8 imagery for various spectral indices (NDVI, EVI, MSAVI, etc.)
2. Calculates long-term trends using linear regression for each index
3. Processes CHIRPS precipitation data and calculates rainfall trends
4. Processes SMAP soil moisture data and calculates soil moisture trends
5. Processes ERA5 relative humidity data and calculates humidity trends
6. Exports trend layers as GEE assets and local GeoTIFF files

Equivalent to the TrendFire.js script in the original implementation.
"""

import ee
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import Polygon
import ee.data # Import the data module

# Initialize the Earth Engine API
def initialize_ee(project_id='ee-crop-health-telangana'):
    try:
        # Attempt authentication. Might prompt user if credentials are not found.
        ee.Authenticate()
        
        # Initialize with the specified project ID.
        ee.Initialize(project=project_id)
        
        print(f"Earth Engine API initialized successfully for project: {project_id}")
    except ee.EEException as e:
        print(f"Error initializing Earth Engine API for project {project_id}: {e}")
        print("Please ensure:")
        print("  1. You have authenticated with `ee.Authenticate()`.")
        print(f"  2. The project ID '{project_id}' is correct.")
        print("  3. The Earth Engine API is enabled for this project in Google Cloud Console.")
        print("  4. The account used for authentication has access to this project.")
        raise # Re-raise the exception to stop execution if initialization fails
    except Exception as e:
        # Catch other potential errors during initialization
        print(f"An unexpected error occurred during Earth Engine initialization: {e}")
        raise

# Function to debug the shapefile
def debug_shapefile(boundary_path):
    """
    Debug function to check shapefile contents and geometry validity.
    """
    try:
        print("\nDebugging shapefile...")
        gdf = gpd.read_file(boundary_path)
        print(f"\nShapefile contains {len(gdf)} features")
        print("\nFirst few rows of the data:")
        print(gdf.head())
        print("\nGeometry type:", gdf.geometry.type.unique())
        print("\nCRS:", gdf.crs)
        
        # Check if geometries are valid
        invalid_geoms = gdf[~gdf.geometry.is_valid]
        if len(invalid_geoms) > 0:
            print("\nWarning: Found invalid geometries!")
            print("Number of invalid geometries:", len(invalid_geoms))
            print("\nInvalid geometries:")
            print(invalid_geoms)
        
        return gdf
    except Exception as e:
        print(f"Error debugging shapefile: {e}")
        return None

# Define the Goa study area boundary
def get_goa_boundary():
    """
    Get the boundary of the study area directly from the GEE asset.
    Returns an ee.Geometry object.
    """
    try:
        asset_id = 'users/jonasnothnagel/pa_boundary'
        print(f"Loading study area boundary from GEE asset: {asset_id}...")
        goa_fc = ee.FeatureCollection(asset_id)
        
        # Ensure the collection is not empty
        if goa_fc.size().getInfo() == 0:
            print(f"Error: Asset {asset_id} is empty or inaccessible.")
            return None
        
        # Get the geometry of the first feature (assuming it's a single polygon feature)
        # Use .geometry() to simplify the FeatureCollection to a single Geometry
        # Use .dissolve() first to merge geometries if it might be a multi-part feature collection
        goa_geometry = goa_fc.geometry().dissolve()
        
        # Optional: Check geometry validity (can be computationally intensive)
        # is_valid = goa_geometry.isValid().getInfo()
        # if not is_valid:
        #     print(f"Warning: Geometry from asset {asset_id} may not be valid according to GEE.")
        
        print("Successfully loaded boundary from GEE asset.")
        return goa_geometry
        
    except ee.EEException as e:
        print(f"Error loading boundary from GEE asset {asset_id}: {e}")
        print("Please ensure the asset exists and you have read permissions.")
        return None
    except Exception as e:
        # Catch other potential errors 
        print(f"An unexpected error occurred loading the boundary asset: {e}")
        return None

# Function to mask clouds in Landsat 8 imagery
def maskL8sr(image):
    """
    Mask clouds and scale pixel values for Landsat 8 SR imagery.
    Matches the JavaScript implementation.
    """
    # Bits 0-4 are Fill, Dilated Cloud, Cirrus, Cloud, Cloud Shadow.
    qaMask = image.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
    saturationMask = image.select('QA_RADSAT').eq(0)

    # Scale the optical bands
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)

    # Scale the thermal bands
    thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)

    # Return the masked and scaled image, applying both masks
    return image.addBands(optical_bands, None, True) \
              .addBands(thermal_bands, None, True) \
              .updateMask(qaMask) \
              .updateMask(saturationMask)

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
    
    # NDFI (Normalized Difference Fraction Index) - Corrected formula
    ndfi = image.normalizedDifference(['ST_B10', 'SR_B6']).rename('ndfi')
    
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

# Main function to process Landsat 8 data and calculate trends
def process_landsat_trends(goa, start_date='2013-03-20', end_date='2023-02-28'):
    """
    Process Landsat 8 imagery, calculate vegetation indices and trends.
    """
    print("Processing Landsat 8 data and calculating trends...")
    
    # Define the outlier date
    outlier_date = ee.Date('2016-06-22')
    
    # Load Landsat 8 collection and filter by date and location
    ls8_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterDate(start_date, end_date) \
        .filterBounds(goa) \
        .filterMetadata('CLOUD_COVER', 'less_than', 10) \
        .map(maskL8sr) \
        .filter(ee.Filter.Or(  # Add outlier filter
            ee.Filter.date(start_date, outlier_date),
            ee.Filter.date(outlier_date.advance(1, 'day'), end_date)
        )) \
        .map(lambda image: image.clip(goa)) \
        .map(addIndices)
    
    # Calculate SWIR min/max for SMI calculation
    swir_stats = ls8_collection.mean().select('SR_B7').reduceRegion(
        reducer=ee.Reducer.minMax(),
        geometry=goa,
        scale=30
    )
    
    swir_min = ee.Number(swir_stats.get('SR_B7_min'))
    swir_max = ee.Number(swir_stats.get('SR_B7_max'))
    
    # Function to calculate SMI using SWIR bands
    def addSMI(image):
        smi = image.select('SR_B7').subtract(swir_min) \
            .divide(swir_max.subtract(swir_min)) \
            .multiply(-1).add(1).rename('smi')
        return image.addBands(smi)
    
    # Add SMI to collection
    ls8_with_smi = ls8_collection.map(addSMI)
    
    # Select bands for trend analysis
    trend_bands = ['ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2', 'msavi', 'smi', 'ST_B10']
    ls8_with_smi = ls8_with_smi.select(trend_bands).sort('system:time_start')
    
    # Add a time band (in years)
    def addTimeBand(image):
        years = ee.Date(image.get('system:time_start')).difference(ee.Date(start_date), 'year')
        return image.addBands(ee.Image.constant(years).float().rename('time'))
    
    collection_with_time = ls8_with_smi.map(addTimeBand)
    
    # Calculate trends for each band and combine into a single image
    all_trends = None
    
    for band in trend_bands:
        print(f"Calculating trend for {band}...")
        
        # Select the index and time bands
        stacked = collection_with_time.select(['time', band])
        
        # Perform linear regression
        regression = stacked.reduce(ee.Reducer.linearFit())
        
        # Extract slope and intercept
        slope = regression.select('scale').rename(f'{band}_Slope')
        intercept = regression.select('offset').rename(f'{band}_Intercept')
        
        # Combine with previous trends
        if all_trends is None:
            all_trends = slope.addBands(intercept)
        else:
            all_trends = all_trends.addBands(slope).addBands(intercept)
    
    # Export the combined trends
    export_to_asset(
        all_trends,
        'users/jonasnothnagel/Trend2023_landsat_py',
        goa,
        'Landsat_Trends_2023_py',
        scale=30
    )
    
    return all_trends, ls8_with_smi

# Function to process CHIRPS precipitation data
def process_chirps_trends(goa, start_date='1982-01-01', end_date='2022-12-31'):
    """
    Process CHIRPS precipitation data and calculate rainfall trends.
    """
    print("Processing CHIRPS precipitation data and calculating trends...")
    
    # Load CHIRPS data
    chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
        .filterDate(start_date, end_date) \
        .filterBounds(goa)
    
    # Function to compute the yearly sum
    def createYearlySum(year):
        # Filter for the specific year
        year_filter = chirps.filter(ee.Filter.calendarRange(year, year, 'year'))
        
        # Create the yearly sum
        yearly_sum = year_filter.sum().clip(goa)
        
        # Add the year as a property
        return yearly_sum.set('year', year) \
                         .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis()) \
                         .rename('precipitation')
    
    # Create a list of yearly sums
    years = ee.List.sequence(1982, 2021)  # Match JS implementation
    yearly_sums = years.map(createYearlySum)
    
    # Convert to ImageCollection
    yearly_sums_collection = ee.ImageCollection.fromImages(yearly_sums)
    
    # Add time band
    def addTimeBand(image):
        time = ee.Date(image.get('system:time_start')).difference(ee.Date('1982-01-01'), 'year')
        return image.addBands(ee.Image.constant(time).rename('time').float())
    
    collection_with_time = yearly_sums_collection.map(addTimeBand)
    
    # Calculate trend
    stacked = collection_with_time.select(['time', 'precipitation'])
    regression = stacked.reduce(ee.Reducer.linearFit())
    slope = regression.select('scale').rename('rain_Slope')
    intercept = regression.select('offset').rename('rain_Intercept')
    trend = slope.addBands(intercept)
    
    # Export the trend
    export_to_asset(
        trend,
        'users/jonasnothnagel/Trend2022_rain_new_py',
        goa,
        'Precipitation_Trend_1982_2022_py',
        scale=30
    )
    
    return trend, yearly_sums_collection

# Function to process SMAP soil moisture data
def process_smap_trends(goa, start_date='2015-04-01', end_date='2023-02-28'):
    """
    Process SMAP soil moisture data and calculate trends.
    Note: SMAP data starts from April 2015.
    """
    print("Processing SMAP soil moisture data and calculating trends...")
    
    # Load SMAP data
    smap = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005') \
        .filterDate(start_date, end_date) \
        .filterBounds(goa) \
        .select(['soil_moisture_am']) \
        .map(lambda image: image.clip(goa))
    
    # Add time band
    def addTimeBand(image):
        time = ee.Date(image.get('system:time_start')).difference(ee.Date('2015-04-01'), 'year')
        return image.addBands(ee.Image.constant(time).rename('time').float())
    
    # Apply time band to each image in the collection
    soil_moisture_with_time = smap.map(addTimeBand)
    
    # Calculate trend
    regression = soil_moisture_with_time.reduce(ee.Reducer.linearFit())
    slope = regression.select('scale').rename('sm_surface_Slope')
    intercept = regression.select('offset').rename('sm_surface_Intercept')
    trend = slope.addBands(intercept)
    
    # Export the trend
    export_to_asset(
        trend,
        'users/jonasnothnagel/Trend2023_SM_new_py',
        goa,
        'SM_Trend2015_2023_py',
        scale=30
    )
    
    return trend, smap

# Function to process ERA5 relative humidity data
def process_era5_trends(goa, start_date='1980-01-01', end_date='2023-02-28'):
    """
    Process ERA5 data and calculate relative humidity trends.
    """
    print("Processing ERA5 relative humidity data and calculating trends...")
    
    # Load ERA5 data
    era5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR") \
        .filterDate(start_date, end_date) \
        .filterBounds(goa) \
        .map(lambda image: image.clip(goa))
    
    # Function to calculate vapor pressure
    def calcVaporPressure(temp):
        denom = temp.add(243.5)
        num = temp.multiply(19.67)  # Note: JS uses 19.67 instead of 17.67
        exponent = num.divide(denom).exp()
        return exponent.multiply(6.112)
    
    # Calculate RH for each time step
    def calculateRH(image):
        dewPoint = image.select('dewpoint_temperature_2m')
        temperature = image.select('temperature_2m')
        
        # Calculate vapor pressures
        eT = calcVaporPressure(temperature)
        eTd = calcVaporPressure(dewPoint)
        
        # Compute relative humidity (%)
        RH = eTd.divide(eT).multiply(100).rename('RH')
        
        # Add time property
        return RH.set('system:time_start', image.get('system:time_start'))
    
    # Apply RH calculation to the collection
    rh_collection = era5.map(calculateRH)
    
    # Add time band
    def addTimeBand(image):
        time = ee.Date(image.get('system:time_start')).difference(ee.Date('1980-01-01'), 'year')
        return image.addBands(ee.Image.constant(time).rename('time').float())
    
    rh_with_time = rh_collection.map(addTimeBand)
    
    # Calculate trend
    regression = rh_with_time.reduce(ee.Reducer.linearFit())
    slope = regression.select('scale').rename('rh_Slope')
    intercept = regression.select('offset').rename('rh_Intercept')
    trend = slope.addBands(intercept)
    
    # Export the trend
    export_to_asset(
        trend,
        'users/jonasnothnagel/Trend2024_RH_new_py',
        goa,
        'RH_Trend1980_2023_py',
        scale=30
    )
    
    return trend, rh_collection

# Function to merge all trend layers
def merge_trend_layers(ls_trends, rain_trend, sm_trend, rh_trend):
    """
    Merge all trend layers into a single multi-band image.
    """
    # Combine all trend layers
    all_trends = ls_trends.addBands(rain_trend).addBands(sm_trend).addBands(rh_trend)
    
    return all_trends

# Function to export the trends to Google Earth Engine assets
def export_to_asset(image, asset_name, region, description=None, scale=30, export_bounds=None):
    """
    Export an image to a Google Earth Engine asset.
    Deletes the asset first if it already exists.
    Allows specifying explicit export bounds.
    """
    if description is None:
        description = asset_name.split('/')[-1]

    # --- Check and delete existing asset --- 
    try:
        asset_info = ee.data.getInfo(asset_name)
        if asset_info:
            print(f"Asset {asset_name} already exists. Deleting...")
            ee.data.deleteAsset(asset_name)
            print(f"Asset {asset_name} deleted.")
    except ee.EEException as e:
        # Handle case where asset doesn't exist (getInfo throws an EEException)
        if 'not found' in str(e):
            print(f"Asset {asset_name} not found. Proceeding with export.")
        else:
            print(f"Error checking/deleting asset {asset_name}: {e}. Proceeding with export attempt...")
    except Exception as e:
        # Catch other potential errors during check/delete
        print(f"Unexpected error checking/deleting asset {asset_name}: {e}. Proceeding with export attempt...")
    # --------------------------------------

    # Use provided bounds if available, otherwise use the region geometry's bounds
    export_region = export_bounds if export_bounds else region

    print(f"Starting export task to Asset: {description}")
    task = ee.batch.Export.image.toAsset(
        image=image, # Image should be pre-projected and cast
        description=description,
        assetId=asset_name,
        region=export_region,
        scale=scale,
        crs='EPSG:3857',
        maxPixels=1e13
    )

    task.start()
    print(f"Task started (id: {task.id}). Check GEE Tasks or use monitoring tools.")
    return task

# Function to export the trends to Google Drive
def export_to_drive(image, filename, region, description=None, folder='GEE_Exports', scale=30, export_bounds=None):
    """
    Export an image to Google Drive.
    Allows specifying explicit export bounds.
    """
    if description is None:
        description = filename
        
    # Use provided bounds if available, otherwise use the region geometry's bounds
    export_region = export_bounds if export_bounds else region
    
    task = ee.batch.Export.image.toDrive(
        image=image, # Image should be pre-projected and cast
        description=description,
        folder=folder,
        fileNamePrefix=filename,
        region=export_region,
        scale=scale,
        crs='EPSG:3857',       
        maxPixels=1e13
    )
    
    task.start()
    print(f"Started export task to Drive: {description}")
    return task

def main():
    """
    Main function to process all trends and export results.
    Exports two versions: one matching inputResampled extent, one covering the full study area.
    """
    # Initialize Earth Engine
    initialize_ee()
    
    # Get Goa boundary
    print("Retrieving Goa boundary...")
    goa = get_goa_boundary()
    if goa is None:
        print("Error: Could not retrieve Goa boundary")
        return
    
    print("Successfully retrieved Goa boundary")
    
    # Process Landsat trends
    print("\nProcessing Landsat trends...")
    landsat_trends, landsat_collection = process_landsat_trends(goa)
    print("Landsat trends processed successfully")
    
    # Process CHIRPS precipitation trends
    print("\nProcessing CHIRPS precipitation trends...")
    rain_trends, rain_collection = process_chirps_trends(goa)
    print("CHIRPS precipitation trends processed successfully")
    
    # Process SMAP soil moisture trends
    print("\nProcessing SMAP soil moisture trends...")
    sm_trends, sm_collection = process_smap_trends(goa)
    print("SMAP soil moisture trends processed successfully")
    
    # Process ERA5 relative humidity trends
    print("\nProcessing ERA5 relative humidity trends...")
    rh_trends, rh_collection = process_era5_trends(goa)
    print("ERA5 relative humidity trends processed successfully")
    
    # Merge all trend layers
    print("\nMerging all trend layers...")
    all_trends = landsat_trends.addBands(rain_trends).addBands(sm_trends).addBands(rh_trends)
    
    # --- Define Target Export Bounds (from inputResampled.tif in EPSG:3857) ---
    target_bounds_coords = [
        [8246820.0, 1680600.0], # bottom-left (x, y)
        [8275140.0, 1680600.0], # bottom-right (x, y)
        [8275140.0, 1767960.0], # top-right (x, y)
        [8246820.0, 1767960.0], # top-left (x, y)
        [8246820.0, 1680600.0]  # close the loop
    ]
    target_export_region_limited = ee.Geometry.Polygon(target_bounds_coords, proj='EPSG:3857', evenOdd=False)

    # --- Explicitly set data type and CRS before export ---
    print("Casting final trend image to float32 and setting projection to EPSG:3857...")
    # Apply casting and projection ONCE to the merged image
    # Exports will use this projected image but different bounding regions
    all_trends_export_ready = all_trends.toFloat() \
                                       .reproject(crs='EPSG:3857', scale=30)

    # --- EXPORT 1: Limited Extent (Matching inputResampled.tif) --- 
    print("\n--- Starting Export for Limited Extent (Match JS Output Grid) ---")
    asset_id_limited = 'users/jonasnothnagel/TrendFirePy_LimitedExtent_py'
    drive_filename_limited = 'TrendFirePy_LimitedExtent'

    # Export Asset (Limited)
    export_to_asset(
        all_trends_export_ready, 
        asset_id_limited,
        goa, # Original region for context, bounds override
        description=f'{asset_id_limited.split("/")[-1]}_Asset', 
        scale=30,
        export_bounds=target_export_region_limited # Use the limited bounds
    )

    # Export Drive (Limited)
    export_to_drive(
        all_trends_export_ready, 
        drive_filename_limited, 
        goa, # Original region for context, bounds override
        description=f'{drive_filename_limited}_Drive', 
        folder='GEE_Exports',      
        scale=30,
        export_bounds=target_export_region_limited # Use the limited bounds
    )

    # --- EXPORT 2: Full Extent (Covering pa_boundary.shp) --- 
    print("\n--- Starting Export for Full Study Area Extent ---")
    asset_id_full = 'users/jonasnothnagel/TrendFirePy_FullExtent_py'
    drive_filename_full = 'TrendFirePy_FullExtent'

    # Export Asset (Full) - Use goa geometry for bounds implicitly
    export_to_asset(
        all_trends_export_ready, 
        asset_id_full,
        goa, # Let GEE determine bounds from this geometry
        description=f'{asset_id_full.split("/")[-1]}_Asset', 
        scale=30
        # No export_bounds specified, defaults to region (goa)
    )

    # Export Drive (Full) - Use goa geometry for bounds implicitly
    export_to_drive(
        all_trends_export_ready, 
        drive_filename_full, 
        goa, # Let GEE determine bounds from this geometry
        description=f'{drive_filename_full}_Drive', 
        folder='GEE_Exports',      
        scale=30
        # No export_bounds specified, defaults to region (goa)
    )

    print("\nAnalysis completed successfully!")
    print("Multiple export tasks started. Monitor GEE Tasks.")
    
    # Return the original merged trends (before final reprojection for export)
    return {
        'trends': all_trends,
        'collections': {
            'landsat': landsat_collection,
            'rain': rain_collection,
            'sm': sm_collection,
            'rh': rh_collection
        }
    }

if __name__ == "__main__":
    main() 