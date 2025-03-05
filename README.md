# ForestFireGoa
Forest fire vulnerability and forest fire prediction using vegetation indices anomaly, surface temperature anomaly, meteorological parameters' anomaly in GEE.

Summary of the methodology:



Sequence of execution:

1. Trendfire.js - ROI is the only required input. 14 trend layers are generaed as output.
2. FireVulnerability.js - Prepare the true fire events file for 2013 to 2023.
3. 2013 - 2021 - true fire events are available from https://gee-community-catalog.org/projects/firms_vector/#notes
4. 2022 and 2023 true fire events from - https://firms.modaps.eosdis.nasa.gov/active_fire/
5. Eg: var fire2019 = ee.FeatureCollection("projects/sat-io/open-datasets/VIIRS/VNP14IMGTDL_NRT_2019")
6. 14 layers from the previous are the inputs for this file.
7. TrendAnomalyPrediction.js - 14 layers as inputs.
Preparing roads layer: Use QuickOSM - to download the roads data for Goa region Use the nearestroad.py - to calculate the distance of the nearest road for every grid cell. "n_dista" column is used to rasterize 9. the resulting vector grid layer = proximity to roads layer. Use IoU.py to check the predicted vs true fire events ratio.
