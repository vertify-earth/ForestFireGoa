# ForestFireGoa

**Forest fire vulnerability and prediction using vegetation indices anomaly, surface temperature anomaly, and meteorological parameters' anomaly in Google Earth Engine (GEE).**

## Overview
This project analyzes forest fire vulnerability and predicts potential fire-prone areas using long-term trends in vegetation indices, burn indices, and meteorological parameters. The methodology leverages Google Earth Engine (GEE) for data processing and analysis, with implementations available in both JavaScript (original) and Python.

## Methodology Summary
1.  **Calculate Decadal Trends:** Generate long-term linear trends (Slope and Intercept) for various environmental indicators (vegetation indices, LST, precipitation, soil moisture, humidity) using historical satellite and climate data.
2.  **Map Fire Vulnerability:** Develop a spatial model using machine learning (Random Forest) trained on historical fire locations, environmental trends, and static factors (topography, roads) to map the inherent susceptibility of different areas to fire.
3.  **Predict Near-Term Fire Risk using Anomalies:** Identify areas where recent environmental conditions deviate significantly from their long-term trends, potentially indicating heightened short-term fire risk.

## Execution Sequence & Scripts

### 1. Trend Calculation (`Trendfire.js` / `TrendFire.py`)
   - **Goal:** Calculate long-term (typically decadal) linear trends (Slope and Intercept) for key environmental indicators.
   - **Input:** Study area boundary (e.g., `pa_boundary` GEE asset).
   - **Datasets Used:** Landsat 8 (Vegetation/Burn Indices, LST), CHIRPS Daily (Precipitation), SMAP (Soil Moisture), ERA5-Land Monthly (Relative Humidity).
   - **Method:** Filters datasets by date and boundary, applies masking (e.g., clouds), calculates relevant indices/parameters, and fits a linear trend over time for each pixel using `ee.Reducer.linearFit()`.
   - **Output:** GEE assets containing the calculated Slope and Intercept bands for each indicator (e.g., `users/.../Trend2023_ndvi`, `users/.../Trend2022_rain_new`). The Python version can optionally export a single merged asset containing all trend bands.

### 2. Fire Vulnerability Mapping (`FireVulnerability.js` / `FireVulnerability.py`)
   - **Goal:** Map the inherent, long-term susceptibility of areas to fire based on environmental factors and historical fire patterns.
   - **Inputs:**
     - *Individual Trend Layers:* Slope and Intercept assets generated in Step 1.
     - *Additional Trend Layers:* Calculated within this script, such as Planet-NICFI visual/NIR/NDVI trends and MODIS LST Day/Night trends (Note: NICFI data access may be restricted).
     - *Static/Quasi-Static Layers:* DEM-derived slope, Roads layer (rasterized).
     - *Historical Fire Points:* GEE assets containing locations of past fires (e.g., `fire13_19`, `fire20_23`).
   - **Method:**
     1.  Loads all input layers (trends, static data).
     2.  Combines these layers into a single multi-band predictor image.
     3.  Resamples the predictor image to a consistent grid (e.g., 30m EPSG:3857). An asset containing this intermediate resampled image (equivalent to the original project's `inputResampled.tif`) is exported by the Python script as `FireVulnerability_InputsResampled_py`.
     4.  Generates random points (labelled non-fire) and processes historical fire points (potentially categorizing them based on attributes like 'Delta T').
     5.  Samples the predictor values from the *resampled* image at the fire/non-fire point locations to create training/validation datasets.
     6.  Trains a Random Forest classifier using the sampled data to learn the relationship between predictor values and fire occurrence/category.
     7.  Applies the trained classifier to the entire *resampled* predictor image.
   - **Output:** A fire risk classification map (e.g., `FireVulnerability_py` asset) showing predicted vulnerability levels across the study area.

### 3. Fire Prediction using Trend Anomalies (`TrendAnomalyPrediction.js` / `TrendAnomalyPrediction.py`)
   - **Goal:** Identify areas behaving abnormally compared to their long-term trends *just before* a period of interest (e.g., the start of fire season), potentially indicating heightened *near-term* fire risk.
   - **Inputs:**
     - Individual Trend layers (Slope/Intercept bands from Step 1).
     - Study area boundary (`pa_boundary` asset).
   - **Methodology:**
     1.  **Load Trends:** Loads the Slope and Intercept assets generated in Step 1.
     2.  **Predict Expected Conditions:** Uses the linear trend formula (`Predicted = Slope * Time_Difference + Intercept`) to estimate the expected value for each indicator at a specific near-future date (e.g., March 1st), based on the trend's reference start date.
     3.  **Load Observed Conditions:** Loads and processes the most recent available observational data *before* the prediction date (e.g., Landsat, SMAP, ERA5 for the preceding month; CHIRPS sum for the preceding year).
     4.  **Calculate Anomaly:** Computes the difference: `Anomaly = Observed - Predicted`.
     5.  **Identify Hotspots:** Applies thresholds to specific anomaly bands (e.g., `ST_B10_Anomaly > 3.2`, `ndmi_Anomaly < -0.15`) to flag areas where conditions deviate significantly from the trend in a way suggesting increased fire risk (e.g., hotter or drier than expected).
   - **Outputs:**
     - *Full Anomaly Image:* GEE asset containing all calculated anomaly bands (e.g., `TrendAnomaly_py`). Pixel values represent the difference between observed and predicted conditions.
     - *Hotspot Images:* Derived boolean/masked images showing areas exceeding specific anomaly thresholds (e.g., `TrendAnomaly_STB10_Hotspot_py`).
   - **Interpretation:** Hotspots indicate areas where recent conditions are unusually conducive to fire compared to that location's historical trend for that time of year. This provides a dynamic, short-term risk indicator, ideally interpreted alongside the longer-term vulnerability map from Step 2.

## Potential Improvements & Next Steps
*   **Trend Model:** Explore non-linear trend models (e.g., harmonic regression) if significant seasonality or non-linear changes are present, potentially providing a more accurate baseline for anomaly calculation.
*   **Data Gaps:** Investigate methods to fill gaps in input data (especially cloud-related gaps in Landsat) using techniques like temporal interpolation or data fusion to allow for trend calculation in more areas.
*   **Vulnerability Model Inputs:** Evaluate the impact of different input features (e.g., including/excluding NICFI/MODIS LST based on availability/quality) on the vulnerability model performance.
*   **Training Data Strategy:** Refine the generation of non-fire points (e.g., using buffers, stratification) and explore different fire point categorization schemes for training the vulnerability model.
*   **Anomaly Thresholds:** Validate and potentially calibrate the anomaly thresholds used for hotspot detection against historical fire data or field observations to improve their local accuracy and relevance.
*   **Temporal Resolution:** Incorporate more frequent or near-real-time data sources if available to reduce the lag between observation and prediction in the anomaly step.
*   **Model Integration:** Combine the trend anomaly results with the vulnerability map to create a more comprehensive risk assessment (e.g., prioritizing areas that are both historically vulnerable *and* currently showing high-risk anomalies).

This workflow enables proactive fire risk management by identifying vulnerable areas based on historical patterns and predictive modeling using trend anomalies.


