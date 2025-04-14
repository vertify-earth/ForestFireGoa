#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TrendAnomalyPrediction.py: Python implementation of TrendAnomalyPrediction.js using Earth Engine Python API.

This script calculates trend anomalies for various indices based on pre-computed trends 
and recent observations.

Steps:
1. Load pre-computed trend layers (Slope and Intercept for various indices).
2. Load the study area boundary.
3. Define prediction date and reference start dates for trends.
4. Calculate predicted index values for the prediction date.
5. Load actual observed data ('present data') for the period preceding the prediction date.
6. Calculate the anomaly (Present - Predicted).
7. Export the anomaly results.
"""

import ee
import os
import datetime
import geopandas as gpd
from shapely.geometry import Polygon

# ==============================================================================
# Helper Functions (Adapted from TrendFire.py)
# ==============================================================================

def initialize_ee(project_id='ee-crop-health-telangana'): # Use the verified project ID
    """Initializes the Earth Engine API."""
    try:
        # Attempt authentication. Might prompt user if credentials are not found.
        ee.Authenticate() # Often needed only once per session/machine
        
        # Initialize with the specified project ID.
        ee.Initialize(project=project_id, opt_url='https://earthengine-highvolume.googleapis.com') # Use high volume endpoint
        
        print(f"Earth Engine API initialized successfully for project: {project_id}")
    except ee.EEException as e:
        print(f"Error initializing Earth Engine API for project {project_id}: {e}")
        print("Please ensure:")
        print("  1. You have authenticated (e.g., via `gcloud auth application-default login`).")
        print(f"  2. The project ID '{project_id}' is correct.")
        print("  3. The Earth Engine API is enabled for this project in Google Cloud Console.")
        print("  4. The account used for authentication has access to this project.")
        raise # Re-raise the exception to stop execution if initialization fails
    except Exception as e:
        # Catch other potential errors during initialization
        print(f"An unexpected error occurred during Earth Engine initialization: {e}")
        raise

def get_goa_boundary(asset_id='users/jonasnothnagel/pa_boundary'):
    """
    Get the boundary of the study area directly from the GEE asset.
    Returns an ee.Geometry object.
    """
    try:
        print(f"Loading study area boundary from GEE asset: {asset_id}...")
        goa_fc = ee.FeatureCollection(asset_id)
        
        # Ensure the collection is not empty
        if goa_fc.size().getInfo() == 0:
            print(f"Error: Asset {asset_id} is empty or inaccessible.")
            return None
        
        # Get the geometry, dissolve potentially multi-part features
        goa_geometry = goa_fc.geometry().dissolve()
        
        print("Successfully loaded boundary from GEE asset.")
        return goa_geometry
        
    except ee.EEException as e:
        print(f"Error loading boundary from GEE asset {asset_id}: {e}")
        print("Please ensure the asset exists and you have read permissions.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred loading the boundary asset: {e}")
        return None

def export_to_asset(image, asset_id, region, description=None, scale=30, crs='EPSG:3857'):
    """Exports an image to a Google Earth Engine asset, deleting if exists."""
    if description is None:
        description = asset_id.split('/')[-1]

    # Check and delete existing asset
    try:
        asset_info = ee.data.getInfo(asset_id)
        if asset_info:
            print(f"Asset {asset_id} already exists. Deleting...")
            ee.data.deleteAsset(asset_id)
            print(f"Asset {asset_id} deleted.")
    except ee.EEException as e:
        if 'not found' in str(e).lower():
            print(f"Asset {asset_id} not found. Proceeding with export.")
        else:
            print(f"Error checking/deleting asset {asset_id}: {e}. Proceeding...")
    except Exception as e:
        print(f"Unexpected error checking/deleting asset {asset_id}: {e}. Proceeding...")

    print(f"Starting export task to Asset: {description} (ID: {asset_id})")
    task = ee.batch.Export.image.toAsset(
        image=image.toFloat(), # Ensure float type for export
        description=description,
        assetId=asset_id,
        region=region,
        scale=scale,
        crs=crs,
        maxPixels=1e13
    )
    task.start()
    print(f"Task started (id: {task.id}). Check GEE Tasks.")
    return task

# Function to debug the shapefile (Optional but good practice)
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
    
    # NDFI (Normalized Difference Fraction Index) - Uses ST_B10 and SR_B6 like JS
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


# Function to calculate SMI using SWIR bands (requires pre-calculated min/max)
def addSMI_local(image, swir_min, swir_max):
    """Calculates SMI, requires swir min/max for the specific collection."""
    # Ensure swir_min and swir_max are ee.Number objects
    swir_min = ee.Number(swir_min)
    swir_max = ee.Number(swir_max)
    # Avoid division by zero or non-numeric results if swir_max <= swir_min
    smi = image.select('SR_B7').subtract(swir_min) \
        .divide(swir_max.subtract(swir_min).max(ee.Number(1e-9))) \
        .multiply(-1).add(1).rename('smi')
    return image.addBands(smi)

# Function to calculate vapor pressure (for RH calculation)
def calcVaporPressure(temp_k):
    """ Calculates vapor pressure from temperature in Kelvin using Clausius-Clapeyron approx. """
    # Convert Kelvin to Celsius for the formula
    temp_c = temp_k.subtract(273.15)
    # Formula parameters from the JS script (note 19.67, standard is often 17.67)
    # Ensure calculations are done per pixel
    denom = temp_c.add(243.5)
    num = temp_c.multiply(19.67)
    exponent = num.divide(denom).exp() # e^(...)
    return exponent.multiply(6.112) #Multiply by saturation vapor pressure at 0°C (in hPa or mbar)

# ==============================================================================
# Main Script Logic
# ==============================================================================

def main():
    """Main script execution."""
    # --- Initialization (Moved to the top) ---
    PROJECT_ID = 'ee-crop-health-telangana' # Or get from environment variable
    initialize_ee(PROJECT_ID)
    
    # --- Configuration ---
    BOUNDARY_ASSET_ID = 'users/jonasnothnagel/pa_boundary'
    OUTPUT_ASSET_BASE = 'users/jonasnothnagel/' # Base path for output assets
    
    # Now create GEE objects after initialization
    PREDICTION_DATE = ee.Date('2023-03-01')
    PRESENT_START_DATE = '2023-02-01' # For LS, SM, RH
    PRESENT_END_DATE = '2023-02-28'   # For LS, SM, RH
    RAIN_YEAR_START = '2022-01-01'    # CHIRPS sum year start
    RAIN_YEAR_END = '2022-12-31'      # CHIRPS sum year end
    
    # Trend reference start dates (match TrendFire.js/py logic)
    TREND_REF_DATES = {
        'landsat': ee.Date('2013-03-20'),
        'rain': ee.Date('1982-01-01'),
        'sm': ee.Date('2015-04-01'),
        'rh': ee.Date('1980-01-01')
    }
    
    # --- Load Boundary (after init) ---
    boundary = get_goa_boundary(BOUNDARY_ASSET_ID)
    if boundary is None:
        print("Exiting: Could not load boundary.")
        return
        
    # --- Processing Steps ---
    # 1. Load Original Trend Assets 
    original_trends = load_original_trend_assets() 
    if original_trends is None:
        print("Exiting: Could not load original trend assets.")
        return
    # print("Placeholder: Skipping load_original_trend_assets()")
    # original_trends = None # Temporary

    # 2. Calculate Predicted Values 
    predicted_image = calculate_predicted_values(original_trends, PREDICTION_DATE, TREND_REF_DATES)
    if predicted_image is None:
        print("Exiting: Could not calculate predicted values.")
        return
    # print("Placeholder: Skipping calculate_predicted_values()")
    # predicted_image = None # Temporary
        
    # 3. Get Present Data
    present_image = get_present_data(boundary, PRESENT_START_DATE, PRESENT_END_DATE, RAIN_YEAR_START, RAIN_YEAR_END)
    if present_image is None:
        print("Exiting: Could not get present data.")
        return
    # print("Placeholder: Skipping get_present_data()")
    # present_image = None # Temporary

    # 4. Calculate Anomaly (Requires present_image and predicted_image)
    if present_image and predicted_image:
        anomaly_image = calculate_anomaly(present_image, predicted_image)
        print("Anomaly image calculated.")
        # print(anomaly_image.bandNames().getInfo()) # Verify band names

        # --- Export Full Anomaly Results ---
        anomaly_asset_id = OUTPUT_ASSET_BASE + 'TrendAnomaly_py' 
        export_to_asset(
            image=anomaly_image,
            asset_id=anomaly_asset_id,
            region=boundary, 
            description='Trend_Anomaly_Python',
            scale=30 
        )
        
        # --- Calculate and Export Specific Hotspots (Matching JS Drive Exports) ---
        print("Calculating and exporting specific hotspots...")
        try:
            # ST_B10 Hotspot (Thermal Anomaly)
            stb10_hotspot = anomaly_image.select('ST_B10_Anomaly').gt(3.2) # Pixels significantly hotter than predicted
            # Masking non-hotspot areas (optional, makes visualization cleaner)
            # stb10_hotspot = stb10_hotspot.selfMask() 
            stb10_asset_id = OUTPUT_ASSET_BASE + 'TrendAnomaly_STB10_Hotspot_py'
            export_to_asset(
                image=stb10_hotspot.unmask(0).toFloat(), # Export boolean as 0/1 float
                asset_id=stb10_asset_id,
                region=boundary,
                description='Hotspot_STB10_py',
                scale=30
            )

            # NDMI Hotspot (Moisture Anomaly)
            ndmi_hotspot = anomaly_image.select('ndmi_Anomaly').lt(-0.15) # Pixels significantly drier than predicted
            # ndmi_hotspot = ndmi_hotspot.selfMask()
            ndmi_asset_id = OUTPUT_ASSET_BASE + 'TrendAnomaly_NDMI_Hotspot_py'
            export_to_asset(
                image=ndmi_hotspot.unmask(0).toFloat(), # Export boolean as 0/1 float
                asset_id=ndmi_asset_id,
                region=boundary,
                description='Hotspot_NDMI_py',
                scale=30
            )
            print("Hotspot export tasks started.")

        except ee.EEException as e:
             print(f"Could not calculate or export hotspots: {e}")
             # Check if specific anomaly bands exist in anomaly_image
             print("Available anomaly bands:", anomaly_image.bandNames().getInfo())
        except Exception as e:
             print(f"Unexpected error during hotspot processing: {e}")

    else:
        print("Skipping anomaly calculation and export due to missing inputs (likely placeholder functions)." ) # Updated message

    print("\nScript finished. Monitor GEE tasks for export status.")

# --- Trend Anomaly Specific Functions ---

def load_original_trend_assets(base_path='users/jonasnothnagel/'):
    """
    Loads the individual trend assets generated by the original TrendFire.js script 
    (assumed to be under the user's GEE path) and merges them.
    Returns a single ee.Image with all trend bands (Slope and Intercept).
    """
    print(f"Loading original trend assets from base path: {base_path}...")
    
    # Define asset suffixes and index names based on TrendFire.js output
    landsat_indices = ['ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2', 'msavi', 'smi', 'ST_B10']
    landsat_asset_prefix = 'Trend2023_' # Corrected from suffix to prefix
    rain_asset_id_part = 'Trend2022_rain_new'
    sm_asset_id_part = 'Trend2023_SM_new'
    rh_asset_id_part = 'Trend2024_RH_new'
    
    all_trend_bands = ee.Image().toFloat() # Start with an empty image, cast to float early
    asset_id = '' # Initialize asset_id for error message context
    
    try:
        # Load Landsat index trends
        for index in landsat_indices:
            asset_id = f"{base_path}{landsat_asset_prefix}{index}"
            print(f"  Loading Landsat trend: {asset_id}")
            trend_img = ee.Image(asset_id)
            # Expecting bands: index_Slope, index_Intercept
            all_trend_bands = all_trend_bands.addBands(trend_img)
                
        # Load Rain trend
        asset_id = f"{base_path}{rain_asset_id_part}"
        print(f"  Loading Rain trend: {asset_id}")
        rain_trend = ee.Image(asset_id) # Bands: rain_Slope, rain_Intercept
        all_trend_bands = all_trend_bands.addBands(rain_trend)

        # Load SM trend
        asset_id = f"{base_path}{sm_asset_id_part}"
        print(f"  Loading SM trend: {asset_id}")
        sm_trend = ee.Image(asset_id) # Bands: sm_surface_Slope, sm_surface_Intercept
        # Ensure band name consistency if different from TrendFire.py
        # sm_trend = sm_trend.rename(['sm_Slope', 'sm_Intercept']) # Example rename if needed
        all_trend_bands = all_trend_bands.addBands(sm_trend)
        
        # Load RH trend
        asset_id = f"{base_path}{rh_asset_id_part}"
        print(f"  Loading RH trend: {asset_id}")
        rh_trend = ee.Image(asset_id) # Bands: rh_Slope, rh_Intercept
        all_trend_bands = all_trend_bands.addBands(rh_trend)
        
        # Remove the empty initial band if necessary (usually not needed with addBands)
        # band_names = all_trend_bands.bandNames()
        # if band_names.length().gt(0) and band_names.get(0).match('^_*'): # Check if first band is empty/default
        #     all_trend_bands = all_trend_bands.select(band_names.slice(1))

        print("Successfully loaded and merged all original trend assets.")
        # Verify band names (optional)
        # print("Loaded trend bands:", all_trend_bands.bandNames().getInfo())
        return all_trend_bands
        
    except ee.EEException as e:
        print(f"Error loading GEE asset during trend loading: {e}")
        print(f"Please ensure asset '{asset_id}' exists and you have permissions.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred loading trend assets: {e}")
        return None

def calculate_predicted_values(trend_image, prediction_date, ref_dates):
    """
    Calculates the predicted values for each index based on trends.
    
    Args:
        trend_image: ee.Image with bands like 'index_Slope', 'index_Intercept'.
        prediction_date: ee.Date object for the prediction target date.
        ref_dates: Dictionary mapping dataset type ('landsat', 'rain', 'sm', 'rh') 
                   to its trend reference start ee.Date.
    
    Returns:
        ee.Image with predicted values for each index (bands named 'index').
    """
    print(f"Calculating predicted values for {prediction_date.format('YYYY-MM-dd').getInfo()}...")

    try:
        # Define the base index names and their corresponding dataset types
        # Note: SM index from JS/Python TrendFire is 'sm_surface', RH is 'RH'
        index_mapping = {
            'ndvi': 'landsat', 'evi': 'landsat', 'mirbi': 'landsat', 'ndfi': 'landsat',
            'bsi': 'landsat', 'ndmi': 'landsat', 'nbr': 'landsat', 'nbr2': 'landsat',
            'msavi': 'landsat', 'smi': 'landsat', 'ST_B10': 'landsat',
            'rain': 'rain',
            'sm_surface': 'sm', # Match SM band name from TrendFire
            'rh': 'rh' # Match RH band name from TrendFire
        }
        
        predicted_bands_list = []
        
        for index, dataset_type in index_mapping.items():
            print(f"  Predicting for index: {index} (type: {dataset_type})")
            
            # Get the specific reference date for this dataset type
            ref_date = ref_dates.get(dataset_type)
            if not ref_date:
                print(f"    Error: Reference date not found for dataset type '{dataset_type}'. Skipping index '{index}'.")
                continue
                
            # Calculate time difference in years
            time_diff = prediction_date.difference(ref_date, 'year')
            
            # Construct band names for slope and intercept
            slope_band = f'{index}_Slope'
            intercept_band = f'{index}_Intercept'
            
            # Check if bands exist in the input image
            if slope_band not in trend_image.bandNames().getInfo() or \
               intercept_band not in trend_image.bandNames().getInfo():
                 print(f"    Error: Bands '{slope_band}' or '{intercept_band}' not found in trend image. Skipping index '{index}'.")
                 continue

            # Select bands
            slope = trend_image.select(slope_band)
            intercept = trend_image.select(intercept_band)
            
            # Calculate predicted value: predicted = slope * time_diff + intercept
            predicted = slope.multiply(time_diff).add(intercept)
            
            # Rename the predicted band to the simple index name
            predicted = predicted.rename(index)
            
            predicted_bands_list.append(predicted)

        # Check if any bands were processed
        if not predicted_bands_list:
            print("Error: No predicted bands were generated. Check input trend image and band names.")
            return None
            
        # Combine all predicted bands into a single image
        predicted_image = ee.ImageCollection(predicted_bands_list).toBands()
        
        # The toBands() function creates names like '0_ndvi', '1_evi'. We need to rename them back.
        original_band_names = [img.bandNames().get(0) for img in predicted_bands_list]
        predicted_image = predicted_image.rename(original_band_names)
        
        print("Predicted values calculated successfully.")
        # print("Predicted bands:", predicted_image.bandNames().getInfo()) # Optional verification
        return predicted_image

    except ee.EEException as e:
        print(f"GEE error during predicted value calculation: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error during predicted value calculation: {e}")
        return None

def get_present_data(boundary, present_start, present_end, rain_year_start, rain_year_end):
    """
    Loads and prepares the 'present' observational data for the period before prediction.
    
    Args:
        boundary: ee.Geometry for clipping.
        present_start: Start date (YYYY-MM-DD) for Landsat, SMAP, ERA5.
        present_end: End date (YYYY-MM-DD) for Landsat, SMAP, ERA5.
        rain_year_start: Start date (YYYY-MM-DD) for the annual rainfall sum year.
        rain_year_end: End date (YYYY-MM-DD) for the annual rainfall sum year.
    
    Returns:
        ee.Image with observed values (bands named 'index').
    """
    print(f"Loading present data ({present_start} to {present_end}, Rain year: {rain_year_start} to {rain_year_end})...")
    
    try:
        # 1. Landsat Present (Feb 2023 Mosaic)
        print("  Processing present-day Landsat...")
        ls_present_col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
            .filterDate(present_start, present_end) \
            .filterBounds(boundary) \
            .map(maskL8sr) \
            .map(lambda img: img.clip(boundary)) \
            .map(addIndices)
        
        # Calculate SWIR min/max for *this specific collection* for SMI
        # Handle potential empty collection after filtering/masking
        collection_size = ls_present_col.size().getInfo()
        present_ls_mosaic = None
        if collection_size > 0:
            swir_stats_present = ls_present_col.mean().select('SR_B7').reduceRegion(
                reducer=ee.Reducer.minMax(),
                geometry=boundary,
                scale=30,
                maxPixels=1e9 # Use lower maxPixels for reduceRegion if needed
            )
            # Handle potential nulls if band doesn't exist or all pixels masked
            swir_min_present = ee.Number(ee.Algorithms.If(swir_stats_present.contains('SR_B7_min'), swir_stats_present.get('SR_B7_min'), 0))
            swir_max_present = ee.Number(ee.Algorithms.If(swir_stats_present.contains('SR_B7_max'), swir_stats_present.get('SR_B7_max'), 1)) # Avoid divide by zero
            
            ls_present_col_smi = ls_present_col.map(lambda img: addSMI_local(img, swir_min_present, swir_max_present))
            
            # Select relevant bands and create mosaic
            present_ls_bands = ['smi', 'ST_B10', 'ndvi', 'evi', 'msavi', 'mirbi', 'ndmi', 'ndfi', 'nbr', 'nbr2', 'bsi']
            present_ls_mosaic = ls_present_col_smi.select(present_ls_bands).mosaic()
            print("  Landsat mosaic created.")
        else:
            print("  Warning: No valid Landsat images found for the present period. Landsat bands will be missing.")
            # Create a placeholder image if needed, or handle downstream

        # 2. ERA5 RH Present (Feb 2023 Mean)
        print("  Processing present-day ERA5 RH...")
        era5_present_col = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR") \
            .filterDate(present_start, present_end) \
            .filterBounds(boundary) \
            .map(lambda img: img.clip(boundary))
            
        def calculateRH(image):
            dewPoint_k = image.select('dewpoint_temperature_2m')
            temp_k = image.select('temperature_2m')
            eT = calcVaporPressure(temp_k)
            eTd = calcVaporPressure(dewPoint_k)
            # Avoid division by zero if eT is 0 and clamp RH
            rh = eTd.divide(eT.max(1e-9)).multiply(100).max(0).min(100) 
            return rh.rename('rh').set('system:time_start', image.get('system:time_start'))

        rh_collection_present = era5_present_col.map(calculateRH)
        # Check if collection is empty before calculating mean
        present_rh = None
        if rh_collection_present.size().getInfo() > 0:
            present_rh = rh_collection_present.mean() # Use mean for monthly agg
            print("  ERA5 RH calculated (mean).")
        else:
             print("  Warning: No ERA5 data found for the present period. RH band will be missing.")
             
        # 3. SMAP Present (Feb 2023 Mosaic)
        print("  Processing present-day SMAP...")
        smap_present_col = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005') \
            .filterDate(present_start, present_end) \
            .filterBounds(boundary) \
            .map(lambda img: img.clip(boundary)) \
            .select(['soil_moisture_am'])
        
        present_sm = None
        if smap_present_col.size().getInfo() > 0:
             present_sm = smap_present_col.mosaic().rename('sm_surface') # Rename to match trend band name
             print("  SMAP mosaic created.")
        else:
            print("  Warning: No SMAP data found for the present period. SM band will be missing.")

        # 4. CHIRPS Rain Present (Sum for 2022)
        print("  Processing present-day CHIRPS Rain (2022 Sum)...")
        chirps_present_col = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
                  .filterDate(rain_year_start, rain_year_end) \
            .filterBounds(boundary)
        
        present_rain = None
        if chirps_present_col.size().getInfo() > 0:
            # Calculate sum for the year 2022
            yearly_sum_2022 = chirps_present_col.sum().clip(boundary).rename('rain')
            present_rain = yearly_sum_2022
            print("  CHIRPS Rain sum calculated.")
        else:
             print("  Warning: No CHIRPS data found for the rain year. Rain band will be missing.")

        # Combine all present-day layers
        print("  Combining present data layers...")
        # Start with an empty image and add bands conditionally
        present_all = ee.Image().toFloat()
        if present_ls_mosaic is not None:
            present_all = present_all.addBands(present_ls_mosaic)
        if present_rain is not None:
            present_all = present_all.addBands(present_rain)
        if present_rh is not None:
            present_all = present_all.addBands(present_rh)
        if present_sm is not None:
            present_all = present_all.addBands(present_sm)
            
        # Check if any bands were added
        if not present_all.bandNames().getInfo():
            print("Error: No present data bands could be generated.")
            return None
            
        print("Present-day data loaded and processed.")
        # Select bands in the same order as predicted_image expects for calculate_anomaly
        # This order comes from the `calculate_predicted_values` function's index_mapping
        final_band_order = [
            'ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2',
            'msavi', 'smi', 'ST_B10', 'rain', 'sm_surface', 'rh'
        ]
        
        # Select only the bands that were successfully created
        available_bands = present_all.bandNames()
        bands_to_select = [b for b in final_band_order if b in available_bands.getInfo()]
        
        if not bands_to_select:
             print("Error: No matching bands found between predicted and present images.")
             return None
             
        present_all = present_all.select(bands_to_select)
        # print("Final present bands:", present_all.bandNames().getInfo()) # Optional verification
        return present_all

    except ee.EEException as e:
        print(f"GEE error during present data retrieval: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error during present data retrieval: {e}")
        return None

def calculate_anomaly(present_image, predicted_image):
    """
    Calculates the anomaly image (Present - Predicted).
    
    Args:
        present_image: ee.Image with observed values ('index' bands).
        predicted_image: ee.Image with predicted values ('index' bands).
        
    Returns:
        ee.Image with anomaly values ('index_Anomaly' bands).
    """
    print("Calculating anomaly image...")
    # Ensure band names match for subtraction
    present_bands = present_image.bandNames()
    predicted_bands = predicted_image.bandNames()
    
    # Assuming predicted_image has the target band names ('ndvi', 'evi', etc.)
    # Select corresponding bands from present_image if needed
    present_image_aligned = present_image.select(predicted_bands)
    
    anomaly = present_image_aligned.subtract(predicted_image)
    
    # Rename bands to include _Anomaly suffix
    new_band_names = predicted_bands.map(lambda b: ee.String(b).cat('_Anomaly'))
    anomaly = anomaly.rename(new_band_names)
    
    return anomaly

# --- Main Execution ---

if __name__ == "__main__":
    # Add dependency check if desired
    try:
        import geopandas
        import shapely
    except ImportError as e:
        print(f"Missing dependency: {e}")
        print("Please install required libraries: pip install geopandas shapely")
    else:
        main() 