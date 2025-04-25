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
### 1. Trend Calculation (`Trendfire.js`)
   - **Input:** `pa_boundary`
   - **Output:** The output consists of 19 trend layers with 2 bands in each layer representing long-term changes in vegetation, burn indices, and meteorological parameters along with 2 constant layers of DEM and road (Proximity to road in tif format) combined into one tif file with 40 bands. This inputResampled.tif file can be generated using the trendFire.js or can be downloaded from here - https://www.dropbox.com/scl/fi/ena1qfeqv5ppshluop8kt/inputResampled.tif?rlkey=arefwl3my7nx3z5qpknu88lsx&st=vj8cn3bd&dl=0

### 2. Fire Vulnerability Mapping (`FireVulnerability.js`)
   - **Inputs:**
     - Trend layers - https://www.dropbox.com/scl/fi/ena1qfeqv5ppshluop8kt/inputResampled.tif?rlkey=arefwl3my7nx3z5qpknu88lsx&st=vj8cn3bd&dl=0
     - Roads
     - DEM (Digital Elevation Model)
     - Fire event datasets (`fire13_19` and `fire20_23`)
   - **Output:**
     - Fire risk classification map (Low/High fire risk) at 30m resolution
     - Kappa coefficient for accuracy assessment

### 3. Fire Prediction (`TrendAnomalyPrediction.js`)
   - **Inputs:**
     - Trend layers - https://www.dropbox.com/scl/fi/ena1qfeqv5ppshluop8kt/inputResampled.tif?rlkey=arefwl3my7nx3z5qpknu88lsx&st=vj8cn3bd&dl=0
     - Fire event datasets (`fire13_19`, `fire20_23`)
     - `pa_boundary`
   - **Output:**
     - Anomaly map comparing trend layers with current conditions

This workflow enables proactive fire risk management by identifying vulnerable areas based on historical patterns and predictive modeling.


