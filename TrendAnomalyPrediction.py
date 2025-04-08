#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TrendAnomalyPrediction.py: Calculates trend anomalies based on long-term trends and recent observations.

This script replicates the core logic of TrendAnomalyPrediction.js:
1. Loads pre-calculated trend layers (output from TrendFire.py).
2. Calculates predicted values for a target period (e.g., Feb 2023) based on trends.
3. Calculates actual observed values for the target period using Landsat, ERA5, SMAP, CHIRPS.
4. Computes the anomaly by subtracting predicted values from observed values.
5. Exports the resulting anomaly image.
"""

import ee
import os
import datetime
import geopandas as gpd
from shapely.geometry import Polygon

# ==============================================================================
# Helper Functions (Adapted from TrendFire.py)
# ==============================================================================

def initialize_ee():
    """Initializes the Earth Engine API."""
    try:
        # Check if already initialized
        try:
            ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005').limit(1).size().getInfo()
            print("Earth Engine API already initialized.")
            return True
        except ee.EEException as e:
             # If not initialized, attempt authentication and initialization
            if "cannot be used before ee.Initialize()" in str(e) or \
               "Earth Engine account" in str(e) or \
               "Credentials are not valid" in str(e):
                try:
                   ee.Authenticate(quiet=True) # Use quiet=True for non-interactive
                   ee.Initialize()
                   print("Earth Engine API initialized successfully.")
                   return True
                except Exception as init_e:
                    print(f"Error initializing Earth Engine API: {init_e}")
                    print("Please ensure you have authenticated with GEE CLI: 'earthengine authenticate'")
                    return False
            else:
                # Handle other EE exceptions if needed
                 print(f"An Earth Engine error occurred: {e}")
                 return False

    except Exception as e:
        print(f"An unexpected error occurred during EE initialization: {e}")
        return False


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

# Define the Goa study area boundary
def get_goa_boundary(shp_path=os.path.join('data', 'pa_boundary.shp')):
    """
    Get the boundary of the study area from a shapefile.
    Returns an ee.Geometry object.
    """
    try:
        if os.path.exists(shp_path):
            print(f"Loading study area boundary from {shp_path}...")
            gdf = gpd.read_file(shp_path) # Removed debug call for cleaner flow

            if len(gdf) == 0:
                print("Error: Shapefile contains no features")
                return None

            # Reproject to EPSG:4326 if necessary for GEE Geometry conversion
            if gdf.crs and gdf.crs.to_epsg() != 4326:
                 print(f"Reprojecting boundary from {gdf.crs} to EPSG:4326...")
                 gdf = gdf.to_crs(epsg=4326)

            geometry = gdf.geometry.iloc[0]

            if not geometry.is_valid:
                 print("Warning: Attempting to fix invalid geometry...")
                 geometry = geometry.buffer(0) # Common trick to fix minor issues
                 if not geometry.is_valid:
                     print("Error: Geometry is invalid and could not be fixed.")
                     return None

            # Convert 3D polygon to 2D if needed
            if geometry.has_z:
                print("Converting 3D polygon to 2D...")
                coords = list(geometry.exterior.coords)
                coords_2d = [(x, y) for x, y, z in coords]
                geometry = Polygon(coords_2d)

            # Convert to GeoJSON format suitable for EE
            gdf.geometry = gpd.GeoSeries([geometry]) # Ensure only the valid geometry is used
            boundary_geojson = gdf.__geo_interface__

            # Extract coordinates for EE Geometry Polygon
            ee_coords = boundary_geojson['features'][0]['geometry']['coordinates']
            return ee.Geometry.Polygon(ee_coords)
        else:
            print(f"Error: Could not find boundary shapefile at {shp_path}")
            return None
    except Exception as e:
        print(f"Error loading boundary from shapefile: {e}")
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
    smi = image.select('SR_B7').subtract(swir_min) \
        .divide(swir_max.subtract(swir_min)) \
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

# Function to export trends to Google Drive
def export_to_drive(image, filename, region, description=None, folder='GEE_Exports_Anomaly', scale=30, export_bounds=None):
    """
    Export an image to Google Drive.
    Assumes image is already cast and projected if needed.
    Uses explicit export bounds if provided.
    """
    if description is None:
        description = filename

    export_region = export_bounds if export_bounds else region

    print(f"Starting export task to Drive: {description}")
    task = ee.batch.Export.image.toDrive(
        image=image, # Assume image is pre-processed (float32, projected)
        description=description,
        folder=folder,
        fileNamePrefix=filename,
        region=export_region.getInfo()['coordinates'], # Pass coordinates explicitly
        scale=scale,
        crs='EPSG:3857', # Ensure CRS matches previous step
        maxPixels=1e13
    )
    task.start()
    print(f"Task started (id: {task.id}). Check GEE Tasks or use monitoring tools.")
    return task

# ==============================================================================
# Main Script Logic
# ==============================================================================

def main():
    """Main function to calculate and export trend anomalies."""
    if not initialize_ee():
        return # Stop if EE initialization fails

    pa = get_goa_boundary()
    if pa is None:
        print("Error: Could not load study area boundary. Exiting.")
        return

    print("Successfully loaded study area boundary.")

    # --- Define Time Periods ---
    trend_start_date = '2013-03-20'
    trend_end_date = '2023-02-28' # Corresponds to TrendFire.py output
    present_start_date = '2023-02-01'
    present_end_date = '2023-02-28'
    rain_present_year_start = '2022-01-01'
    rain_present_year_end = '2022-12-31'
    rh_sm_start_date = '1980-01-01' # Start date for RH trend baseline
    sm_start_date = '2015-04-01' # Start date for SM trend baseline
    rain_start_date = '1982-01-01' # Start date for Rain trend baseline


    # --- Load Trend Image Asset ---
    trend_asset_id = 'users/jonasnothnagel/Trend2024_all_new' # Asset from TrendFire.py
    print(f"Loading trend image asset: {trend_asset_id}")
    try:
        input_trends = ee.Image(trend_asset_id)
        # Verify bands (optional)
        # print("Trend bands:", input_trends.bandNames().getInfo())
    except Exception as e:
        print(f"Error loading trend asset {trend_asset_id}: {e}")
        print("Please ensure the asset exists and you have access.")
        return

    # --- (Optional) Load Fire Points ---
    # Assuming Asset IDs - replace with actual IDs if available/needed
    # print("Loading fire points...")
    # try:
    #     fire13_19_asset = 'YOUR_ASSET_ID/fire13_19' # Replace
    #     fire20_23_asset = 'YOUR_ASSET_ID/fire20_23' # Replace
    #     fire13_19 = ee.FeatureCollection(fire13_19_asset)
    #     fire20_23 = ee.FeatureCollection(fire20_23_asset)
    #     fire13_23 = fire13_19.merge(fire20_23)
    #     # Filter as in JS if needed for visualization/context
    #     fireMarch2023 = fire13_23.filter(ee.Filter.stringContains('acq_date', '2023')) # Simplified filter
    #     print(f"Loaded {fire13_23.size().getInfo()} total fire points.")
    # except Exception as e:
    #     print(f"Warning: Could not load fire point assets: {e}. Continuing without them.")
    #     fireMarch2023 = None # Set to None if loading fails


    print("\nCalculating predicted values based on trends...")
    # --- Calculate Predicted Values (Past_Pred) ---
    target_date_landsat = ee.Date(present_start_date) # Use start of month for prediction point
    target_date_rain = ee.Date(rain_present_year_start)
    target_date_rh = target_date_landsat
    target_date_sm = target_date_landsat

    # Calculate time differences in years
    const_landsat_years = target_date_landsat.difference(ee.Date(trend_start_date), 'year')
    const_rain_years = target_date_rain.difference(ee.Date(rain_start_date), 'year')
    const_rh_years = target_date_rh.difference(ee.Date(rh_sm_start_date), 'year') # JS used constRain, seems incorrect, use RH baseline
    const_sm_years = target_date_sm.difference(ee.Date(sm_start_date), 'year') # JS used constRain, seems incorrect, use SM baseline

    # Select slope and intercept bands - Ensure names match TrendFire.py output
    landsat_indices = ['ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2', 'msavi', 'smi', 'ST_B10']
    slope_bands_ls = [f'{idx}_Slope' for idx in landsat_indices]
    intercept_bands_ls = [f'{idx}_Intercept' for idx in landsat_indices]
    other_indices = ['rain', 'sm_surface', 'rh']
    slope_bands_other = [f'{idx}_Slope' for idx in other_indices]
    intercept_bands_other = [f'{idx}_Intercept' for idx in other_indices]

    # Select bands from the input trend image
    slopes_ls = input_trends.select(slope_bands_ls)
    intercepts_ls = input_trends.select(intercept_bands_ls)
    slopes_other = input_trends.select(slope_bands_other)
    intercepts_other = input_trends.select(intercept_bands_other)

    # Calculate predicted values for Landsat-based indices
    predicted_ls = slopes_ls.multiply(const_landsat_years).add(intercepts_ls).rename(landsat_indices)

    # Calculate predicted values for other indices
    pred_rain = slopes_other.select('rain_Slope').multiply(const_rain_years).add(intercepts_other.select('rain_Intercept')).rename('rain')
    pred_sm = slopes_other.select('sm_surface_Slope').multiply(const_sm_years).add(intercepts_other.select('sm_surface_Intercept')).rename('sm')
    pred_rh = slopes_other.select('rh_Slope').multiply(const_rh_years).add(intercepts_other.select('rh_Intercept')).rename('rh')

    # Combine all predicted values
    past_Pred = predicted_ls.addBands([pred_rain, pred_sm, pred_rh])
    # Reorder bands to match expected order for subtraction later
    ordered_pred_bands = ['smi', 'ST_B10', 'ndvi', 'evi', 'msavi', 'mirbi', 'ndmi', 'ndfi', 'nbr', 'nbr2', 'bsi', 'rain', 'rh', 'sm']
    past_Pred = past_Pred.select(ordered_pred_bands)
    print("Predicted values calculated.")
    # print("Predicted bands:", past_Pred.bandNames().getInfo())


    print("\nCalculating present-day observed values...")
    # --- Calculate Present-Day Values (Present_All) ---

    # 1. Landsat Present (Feb 2023 Mosaic)
    print("Processing present-day Landsat (Feb 2023)...")
    present_ls_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterDate(present_start_date, present_end_date) \
        .filterBounds(pa) \
        .map(maskL8sr) \
        .map(lambda img: img.clip(pa)) \
        .map(addIndices)

    # Calculate SWIR min/max for *this specific collection* for SMI
    swir_stats_present = present_ls_collection.mean().select('SR_B7').reduceRegion(
        reducer=ee.Reducer.minMax(),
        geometry=pa,
        scale=30,
        maxPixels=1e9 # Use lower maxPixels for reduceRegion if needed
    )
    # Handle potential nulls if collection is empty
    swir_min_present = ee.Number(ee.Algorithms.If(swir_stats_present.contains('SR_B7_min'), swir_stats_present.get('SR_B7_min'), 0))
    swir_max_present = ee.Number(ee.Algorithms.If(swir_stats_present.contains('SR_B7_max'), swir_stats_present.get('SR_B7_max'), 1)) # Avoid divide by zero

    present_ls_collection_smi = present_ls_collection.map(lambda img: addSMI_local(img, swir_min_present, swir_max_present))

    # Select relevant bands and create mosaic
    present_ls_bands = ['smi', 'ST_B10', 'ndvi', 'evi', 'msavi', 'mirbi', 'ndmi', 'ndfi', 'nbr', 'nbr2', 'bsi']
    present_ls = present_ls_collection_smi.select(present_ls_bands).mosaic()
    print("Landsat mosaic created.")

    # 2. ERA5 RH Present (Feb 2023 Mean)
    print("Processing present-day ERA5 RH (Feb 2023)...")
    era5_present = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR") \
        .filterDate(present_start_date, present_end_date) \
        .filterBounds(pa) \
        .map(lambda img: img.clip(pa))

    def calculateRH(image):
        dewPoint_k = image.select('dewpoint_temperature_2m')
        temp_k = image.select('temperature_2m')
        eT = calcVaporPressure(temp_k)
        eTd = calcVaporPressure(dewPoint_k)
        # Avoid division by zero if eT is 0
        rh = eTd.divide(eT).multiply(100).max(0).min(100) # Clamp RH between 0 and 100
        return rh.rename('rh').set('system:time_start', image.get('system:time_start'))

    rh_collection_present = era5_present.map(calculateRH)
    present_rh = rh_collection_present.mean() # Use mean for monthly agg
    print("ERA5 RH calculated.")


    # 3. SMAP Present (Feb 2023 Mosaic)
    print("Processing present-day SMAP (Feb 2023)...")
    smap_present_col = ee.ImageCollection('NASA/SMAP/SPL3SMP_E/005') \
        .filterDate(present_start_date, present_end_date) \
        .filterBounds(pa) \
        .map(lambda img: img.clip(pa)) \
        .select(['soil_moisture_am'])
    present_sm = smap_present_col.mosaic().rename('sm')
    print("SMAP mosaic created.")


    # 4. CHIRPS Rain Present (Sum for 2022)
    print("Processing present-day CHIRPS Rain (2022 Sum)...")
    chirps_present_col = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
              .filterDate(rain_present_year_start, rain_present_year_end) \
              .filterBounds(pa)

    # Calculate sum for the year 2022
    yearly_sum_2022 = chirps_present_col.sum().clip(pa).rename('rain')
    present_rain = yearly_sum_2022
    print("CHIRPS Rain sum calculated.")

    # Combine all present-day layers
    present_all = present_ls.addBands([present_rain, present_rh, present_sm])
    # Select bands in the same order as past_Pred
    present_all = present_all.select(ordered_pred_bands)
    print("Present-day layers combined.")
    # print("Present bands:", present_all.bandNames().getInfo())


    print("\nCalculating anomaly image...")
    # --- Calculate Anomaly Image ---
    anomalyImage = present_all.subtract(past_Pred)

    # Rename bands to include '_Anomaly'
    anomaly_band_names = [f'{b}_Anomaly' for b in ordered_pred_bands]
    anomalyImage = anomalyImage.rename(anomaly_band_names)
    print("Anomaly image calculated.")
    # print("Anomaly bands:", anomalyImage.bandNames().getInfo())


    # --- Define Target Export Bounds (matching TrendFire.py output) ---
    # Ensures consistent grid with previous step's output
    target_bounds_coords = [
        [8246820.0, 1680600.0], # bottom-left (x, y)
        [8275140.0, 1680600.0], # bottom-right (x, y)
        [8275140.0, 1767960.0], # top-right (x, y)
        [8246820.0, 1767960.0], # top-left (x, y)
        [8246820.0, 1680600.0]  # close the loop
    ]
    target_export_region = ee.Geometry.Polygon(target_bounds_coords, proj='EPSG:3857', evenOdd=False)


    # --- Export Anomaly Image ---
    print("\nExporting anomaly image...")
    export_to_drive(
        image=anomalyImage.toFloat(), # Ensure float32 for export consistency
        filename='TrendAnomalyPy_output',
        region=pa, # Use original geometry for context, bounds dictate output extent
        description='Trend_Anomaly_Python_Output',
        folder='GEE_Exports_Anomaly',
        scale=30,
        export_bounds=target_export_region
    )

    print("\nScript finished. Monitor GEE Tasks for export completion.")


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