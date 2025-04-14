# ForestFireGoa

**Forest fire vulnerability and prediction using vegetation indices anomaly, surface temperature anomaly, and meteorological parameters' anomaly in Google Earth Engine (GEE).**

## Overview
This project analyzes forest fire vulnerability and predicts potential fire-prone areas using long-term trends in vegetation indices, burn indices, and meteorological parameters. The methodology leverages Google Earth Engine (GEE) for data processing and analysis.

## Methodology
1. **Calculate Decadal Trends**  
   - Generate long-term trends for vegetation and burn indices using Landsat imagery.
   - Compute trends for meteorological parameters such as rainfall and relative humidity.
2. **Fire Vulnerability Mapping**  
   - Use historical fire events and trend layers to assess vulnerability.
   - Identify areas most susceptible to fires.
3. **Fire Prediction**  
   - Utilize trend layers to predict high-risk fire zones in the near future.

## Execution Sequence

**Note on `inputResampled.tif`:** The original README linked to a file named `inputResampled.tif` for multiple steps. Analysis of the code reveals this file is actually an **intermediate output** created *within* the Fire Vulnerability Mapping step (Step 2). It represents a combination and resampling of *all* features used in that step (including TrendFire outputs, NICFI, MODIS LST, DEM slope, roads). It is **not** a direct output of Step 1 (`TrendFire`) and **not** an input to Step 3 (`TrendAnomalyPrediction`). The descriptions below reflect the actual inputs/outputs based on the code logic.

### 1. Trend Calculation (`Trendfire.js` / `TrendFire.py`)
   - **Goal:** Calculate long-term (decadal) linear trends (Slope and Intercept) for various environmental indicators.
   - **Input:** Study area boundary (`pa_boundary` asset).
   - **Datasets Used:** Landsat 8 (Vegetation/Burn Indices, LST), CHIRPS Daily (Precipitation), SMAP (Soil Moisture), ERA5-Land Monthly (Relative Humidity).
   - **Method:** Filters datasets by date and boundary, applies masking, calculates indices, fits a linear trend (`ee.Reducer.linearFit()`).
   - **Outputs:** Individual GEE assets for each calculated trend (e.g., `users/.../Trend2023_ndvi` (Slope+Intercept), `users/.../Trend2022_rain_new` (Slope+Intercept), etc.). The Python version can optionally export a single merged asset containing all trend bands.

### 2. Fire Vulnerability Mapping (`FireVulnerability.js` / `FireVulnerability.py`)
   - **Goal:** Map the inherent, long-term susceptibility of areas to fire based on environmental factors and historical fire patterns.
   - **Inputs:**
     - Individual Trend layers (from Step 1).
     - Additional Trend Layers (calculated within this script): NICFI (R/G/B/N/NDVI trends), MODIS LST (Day/Night trends).
     - Static/Quasi-Static Layers: DEM slope, Roads layer.
     - Historical Fire Points (e.g., `fire13_19`, `fire20_23` assets).
   - **Method:**
     1.  Loads all input layers.
     2.  Combines all input layers/bands into a single multi-band image.
     3.  Resamples this combined image to a consistent grid (e.g., 30m EPSG:3857) - *This resampled image is the equivalent of the data in `inputResampled.tif` and is exported by the Python script as `FireVulnerability_InputsResampled_py`.*
     4.  Generates random points (labelled non-fire) and categorizes historical fire points (e.g., based on 'Delta T' property).
     5.  Samples the *resampled* multi-band image at the fire/non-fire point locations to create training/validation data.
     6.  Trains a classifier (e.g., Random Forest) using the sampled data.
     7.  Applies the trained classifier to the *resampled* multi-band image covering the entire study area.
   - **Output:** Fire risk classification map (e.g., `FireVulnerability_py` asset) showing predicted vulnerability levels (e.g., Low/High risk) across the study area.

### 3. Fire Prediction using Trend Anomalies (`TrendAnomalyPrediction.js` / `TrendAnomalyPrediction.py`)
   - **Goal:** Identify areas behaving abnormally compared to their long-term trends *just before* a period of interest, potentially indicating heightened *near-term* fire risk.
   - **Inputs:**
     - Individual Trend layers (Slope/Intercept bands from Step 1).
     - Study area boundary (`pa_boundary` asset).
   - **Methodology:**
     1.  **Load Trends:** Loads the individual trend assets (Slopes and Intercepts) generated in Step 1.
     2.  **Predict Expected Conditions:** Uses the linear trend formula (`Predicted = Slope * Time_Difference + Intercept`) to predict the expected value for each indicator for a near-future date (e.g., March 1st, 2023), based on the reference start date of each trend.
     3.  **Load Observed Conditions:** Loads and processes the most recent available satellite/climate data before the prediction date (e.g., Landsat, SMAP, ERA5 for Feb 2023; CHIRPS sum for 2022).
     4.  **Calculate Anomaly:** Computes the difference between the observed conditions and the predicted conditions (`Anomaly = Observed - Predicted`).
     5.  **Identify Hotspots:** Applies thresholds to specific anomaly bands (e.g., `ST_B10_Anomaly > 3.2`, `ndmi_Anomaly < -0.15`) to flag areas where conditions deviate significantly from the trend in a way that increases fire risk (e.g., significantly hotter or drier than expected).
   - **Outputs:**
     - **Full Anomaly Image:** GEE asset containing anomaly bands (e.g., `TrendAnomaly_py`, `TrendAnomaly_js`). Pixel values represent the difference between observed and predicted conditions.
     - **Hotspot Images:** Derived boolean/masked images showing areas exceeding specific anomaly thresholds (e.g., `TrendAnomaly_STB10_Hotspot_py`).
   - **Interpretation:** Hotspots indicate areas where recent conditions are unusually conducive to fire compared to the established long-term trend for that location and time of year. Best interpreted alongside the vulnerability map from Step 2.

## Potential Improvements & Next Steps
Based on the anomaly prediction methodology:
*   **Trend Model:** Explore non-linear trend models (e.g., harmonic regression) if significant seasonality or non-linear changes are present in the data, potentially providing a more accurate baseline for anomaly calculation.
*   **Data Gaps:** Investigate methods to fill gaps in input data (especially cloud-related gaps in Landsat) using techniques like temporal interpolation or data fusion (e.g., using MODIS) to allow for trend calculation in more areas.
*   **Threshold Validation:** Validate and potentially calibrate the anomaly thresholds used for hotspot detection against historical fire data or field observations to improve their local accuracy.
*   **Temporal Resolution:** Incorporate more frequent or near-real-time data sources if available to reduce the lag between observation and prediction.
*   **Input Variables:** Evaluate the inclusion of additional static or dynamic variables known to influence fire risk (e.g., detailed fuel type maps, wind data, human activity indicators) in either the trend analysis or the final vulnerability/prediction models.
*   **Model Integration:** Combine the trend anomaly results with the vulnerability map from Step 2 to create a more comprehensive risk assessment (e.g., prioritizing areas that are both historically vulnerable *and* currently showing high-risk anomalies).

This workflow enables proactive fire risk management by identifying vulnerable areas based on historical patterns and predictive modeling using trend anomalies.


