import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Load shapefiles
predicted = gpd.read_file("D:/EarthAnalytics/ForestFire/FeltOutputsPresentation/fireHotspot/fireHSPoly.shp")
true_fire = gpd.read_file("D:/EarthAnalytics/ForestFire/FeltOutputsPresentation/2023fire.shp")

# Ensure both layers use a Projected CRS (change EPSG as per your region)
projected_crs = "EPSG:32643"  
predicted = predicted.to_crs(projected_crs)
true_fire = true_fire.to_crs(projected_crs)

# Buffer both shapefiles by 375m
predicted_buffered = predicted.copy()
true_fire_buffered = true_fire.copy()
predicted_buffered["geometry"] = predicted.geometry.buffer(500)
true_fire_buffered["geometry"] = true_fire.geometry.buffer(500)

# Debugging: Plot to check overlap
fig, ax = plt.subplots(figsize=(6, 8))
predicted_buffered.plot(ax=ax, color="red", alpha=0.5, label="Predicted (Buffered)")
true_fire_buffered.plot(ax=ax, color="blue", alpha=0.5, label="True Fire (Buffered)")
# Add labels
ax.set_xlabel("Easting (UTM)")  # Adjust if using lat/lon
ax.set_ylabel("Northing (UTM)")  # Adjust if using lat/lon
ax.set_title("Predicted Fire vs True Fire")

# Create legend manually
red_patch = mpatches.Patch(color='red', label="Buffered Predicted Fire")
blue_patch = mpatches.Patch(color='blue', label="Buffered True Fire")
ax.legend(handles=[red_patch, blue_patch], loc="upper right")
plt.legend()
# Add legend
ax.legend()
#plt.title("Predicted Fire vs True Fire")
plt.show()

# Compute Intersection and Union
intersection = gpd.overlay(predicted_buffered, true_fire_buffered, how="intersection", keep_geom_type=False)
union = gpd.overlay(predicted_buffered, true_fire_buffered, how="union", keep_geom_type=False)

# Debugging: Check if areas are non-zero
print(f"Intersection Area: {intersection.area.sum()}")
print(f"Union Area: {union.area.sum()}")

# Compute IoU
if union.area.sum() > 0:
    iou = intersection.area.sum() / union.area.sum()
    print(f"IoU: {iou:.4f}")
else:
    print("Warning: Union area is zero. Check CRS and data overlap.")

# Save buffered shapefiles (optional)
predicted_buffered.to_file("D:/EarthAnalytics/ForestFire/FeltOutputsPresentation/fireHS_vegBurn_375mBuffer.shp")
true_fire_buffered.to_file("D:/EarthAnalytics/ForestFire/FeltOutputsPresentation/2023fire_375mBuffer.shp")
