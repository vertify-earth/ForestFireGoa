# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 20:31:33 2025

@author: srava
"""

import geopandas as gpd

# Load the grid and roads layers
joinedGrid = gpd.read_file("D:/EarthAnalytics/ForestFire/grid_with_nearest_buffer_optimized2.shp")  # Replace with actual path
paRoadReproj = gpd.read_file("D:/EarthAnalytics/ForestFire/paRoadReproj.shp")  # Replace with actual path

# Ensure both layers have the same CRS
if joinedGrid.crs != paRoadReproj.crs:
    paRoadReproj = paRoadReproj.to_crs(joinedGrid.crs)

# Perform spatial join (intersection)
joined = gpd.sjoin(joinedGrid, paRoadReproj, how="left", predicate="intersects")

road_counts = joined.groupby(joined.index).size()

# Rename index column to match joinedGrid index
road_counts = road_counts.rename("road_count")

# Ensure the grid has a proper index column for merging
joinedGrid["grid_id"] = joinedGrid.index

# Merge the count back to the original grid
joinedGrid = joinedGrid.merge(road_counts.rename_axis("grid_id").reset_index(), 
                              left_index=True, right_on="grid_id", how="left").fillna({"road_count": 0})

joinedGrid.to_file("D:/EarthAnalytics/ForestFire/delete3.shp")