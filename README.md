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
### 1. Trend Calculation (`Trendfire.js` / `TrendFire.py`)
   - Calculates long-term (decadal) linear trends (Slope and Intercept) for various environmental indicators within the study area (`pa_boundary`).
   - **Datasets Used:** Landsat 8 (Vegetation/Burn Indices, LST), CHIRPS Daily (Precipitation), SMAP (Soil Moisture), ERA5-Land Monthly (Relative Humidity).
   - **Method:** Filters datasets by date and boundary, applies relevant masking (e.g., cloud masking for Landsat), calculates indices, fits a linear trend over time (`ee.Reducer.linearFit()`), and exports the resulting slope/intercept bands.
   - **Outputs (JS Version):** Individual GEE assets for each index trend (e.g., `users/.../Trend2023_ndvi`, `users/.../Trend2022_rain_new`).
   - **Outputs (Python Version):** Merged GEE assets containing all trend bands (e.g., `users/.../TrendFirePy_FullExtent_py`) and optionally individual trend assets.

### 2. Fire Vulnerability Mapping (`FireVulnerability.js` / `FireVulnerability.py` - *Python version TBD*)
   - **Inputs:**
     - Trend layers (often combined into a single multi-band image like `inputResampled.tif`)
     - Roads layer
     - DEM (Digital Elevation Model) derived slope
     - Historical fire event datasets (`fire13_19`, `fire20_23`)
     - Other potential factors (e.g., Planet NICFI, MODIS LST)
   - **Method:** Samples input features at locations of historical fires and random non-fire points. Trains a classifier (e.g., Random Forest) to predict fire risk categories based on the input features. Applies the trained classifier to the entire study area.
   - **Output:** Fire risk classification map (e.g., Low/Moderate/High risk) at 30m resolution.

### 3. Fire Prediction using Trend Anomalies (`TrendAnomalyPrediction.js` / `TrendAnomalyPrediction.py`)
   - **Goal:** Identify areas behaving abnormally compared to their long-term trends, potentially indicating heightened near-term fire risk.
   - **Inputs:**
     - Individual Trend layers (Slope/Intercept bands from Step 1, e.g., `users/.../Trend2023_ndvi`)
     - Study area boundary (`pa_boundary`)
   - **Methodology:**
     1.  **Load Trends:** Loads the individual trend assets (Slopes and Intercepts) generated in Step 1.
     2.  **Predict Expected Conditions:** Uses the linear trend formula (`Predicted = Slope * Time_Difference + Intercept`) to predict the expected value for each indicator for a near-future date (e.g., March 1st, 2023), based on the reference start date of each trend.
     3.  **Load Observed Conditions:** Loads and processes the most recent available satellite/climate data before the prediction date (e.g., Landsat, SMAP, ERA5 for Feb 2023; CHIRPS sum for 2022).
     4.  **Calculate Anomaly:** Computes the difference between the observed conditions and the predicted conditions (`Anomaly = Observed - Predicted`).
     5.  **Identify Hotspots:** Applies thresholds to specific anomaly bands (e.g., `ST_B10_Anomaly > 3.2`, `ndmi_Anomaly < -0.15`) to flag areas where conditions deviate significantly from the trend in a way that increases fire risk (e.g., significantly hotter or drier than expected).
   - **Outputs:**
     - **Full Anomaly Image (Python & modified JS):** GEE asset containing anomaly bands (e.g., `ndvi_Anomaly`, `ST_B10_Anomaly`). Pixel values represent the difference between observed and predicted conditions.
     - **Hotspot Images (Python & JS):** Derived boolean/masked images showing areas exceeding specific anomaly thresholds (e.g., `TrendAnomaly_STB10_Hotspot_py`, `TrendAnomaly_NDMI_Hotspot_py`).
   - **Interpretation:** Hotspots indicate areas where recent conditions are unusually conducive to fire compared to the established long-term trend for that location and time of year. For example:
     - *Thermal Hotspot (`ST_B10_Anomaly > threshold`):* Area is significantly hotter than expected.
     - *Moisture Hotspot (`ndmi_Anomaly < threshold`):* Area is significantly drier than expected.

## Potential Improvements & Next Steps
Based on the anomaly prediction methodology:
*   **Trend Model:** Explore non-linear trend models (e.g., harmonic regression) if significant seasonality or non-linear changes are present in the data, potentially providing a more accurate baseline for anomaly calculation.
*   **Data Gaps:** Investigate methods to fill gaps in input data (especially cloud-related gaps in Landsat) using techniques like temporal interpolation or data fusion (e.g., using MODIS) to allow for trend calculation in more areas.
*   **Threshold Validation:** Validate and potentially calibrate the anomaly thresholds used for hotspot detection against historical fire data or field observations to improve their local accuracy.
*   **Temporal Resolution:** Incorporate more frequent or near-real-time data sources if available to reduce the lag between observation and prediction.
*   **Input Variables:** Evaluate the inclusion of additional static or dynamic variables known to influence fire risk (e.g., detailed fuel type maps, wind data, human activity indicators) in either the trend analysis or the final vulnerability/prediction models.
*   **Model Integration:** Combine the trend anomaly results with the vulnerability map from Step 2 to create a more comprehensive risk assessment (e.g., prioritizing areas that are both historically vulnerable *and* currently showing high-risk anomalies).

This workflow enables proactive fire risk management by identifying vulnerable areas based on historical patterns and predictive modeling using trend anomalies.


