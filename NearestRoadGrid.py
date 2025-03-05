# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 19:38:16 2025

@author: srava
"""

# -*- coding: utf-8 -*-
"""
Optimized Spatial Join: Assigns Nearest Buffer Category to Grid using cKDTree.
"""

import geopandas as gpd
import numpy as np
from scipy.spatial import cKDTree

# Load the grid and buffer layers
grid = gpd.read_file("D:/EarthAnalytics/ForestFire/grid30m.shp")  
buffers = gpd.read_file("D:/EarthAnalytics/ForestFire/multiringBuffer200.shp")  

# Ensure both layers have the same CRS
grid = grid.to_crs("EPSG:32643")
buffers = buffers.to_crs(grid.crs)

# Get centroids of both layers for nearest neighbor search
grid["centroid"] = grid.geometry.centroid
buffers["centroid"] = buffers.geometry.centroid

# Convert centroids to NumPy arrays for fast querying
grid_points = np.array(list(zip(grid["centroid"].x, grid["centroid"].y)))
buffer_points = np.array(list(zip(buffers["centroid"].x, buffers["centroid"].y)))

# Build the KDTree for buffer centroids
tree = cKDTree(buffer_points)

# Query the nearest buffer for each grid cell
distances, indices = tree.query(grid_points, k=1)  # k=1 → Find only the nearest neighbor

# Assign nearest buffer category & distance to the grid
nearest_buffers = buffers.loc[indices, ["ringId", "distance"]].to_numpy()
grid["nearest_buffer"], grid["near_distance"] = nearest_buffers[:, 0], nearest_buffers[:, 1]

grid["n_dist"] = distances  # Assign actual nearest computed distance


# Drop temporary centroid columns
grid = grid.drop(columns=["centroid"])

# Save output
grid.to_file("D:/EarthAnalytics/ForestFire/grid_with_nearest_buffer_optimized2.shp")

print("Optimized nearest buffer assignment completed successfully.")
