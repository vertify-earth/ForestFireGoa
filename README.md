# ForestFireGoa
Forest fire vulnerability and forest fire prediction using vegetation indices anomaly, surface temperature anomaly, meteorological parameters' anomaly in GEE.

Summary of the methodology:



Sequence of execution:

1. Trendfire.js - pa_boundary is the only required input. 14 trend layers are generaed as output.
2. FireVulnerability.js - Add trend layers (generated in the previous step), road, dem, fire13_19 and fire20_23 as inputs in the fireVulnerability.js
3. TrendAnomalyPrediction.js - Add trend layers, fire13_19, fire20_23 and pa_boundary as inputs.


