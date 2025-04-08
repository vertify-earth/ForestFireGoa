import rasterio
import os
import numpy as np

def compare_tifs(file1_path, file2_path):
    """
    Compares metadata and basic statistics of two GeoTIFF files.
    """
    print(f"Comparing:\n  (1) {os.path.basename(file1_path)}\n  (2) {os.path.basename(file2_path)}\n")

    if not os.path.exists(file1_path):
        print(f"Error: File not found - {file1_path}")
        return
    if not os.path.exists(file2_path):
        print(f"Error: File not found - {file2_path}")
        return

    try:
        with rasterio.open(file1_path) as ds1, rasterio.open(file2_path) as ds2:
            
            # --- Metadata Comparison ---
            print("--- Metadata Comparison ---")
            
            meta1 = ds1.meta
            meta2 = ds2.meta

            print(f"Number of Bands:  (1) {ds1.count:<5} (2) {ds2.count}")
            print(f"Data Types:       (1) {ds1.dtypes[0]:<10} (2) {ds2.dtypes[0]}") # Show first band's dtype
            print(f"CRS:              (1) {ds1.crs.to_string() if ds1.crs else 'N/A':<10} (2) {ds2.crs.to_string() if ds2.crs else 'N/A'}")
            print(f"Resolution (X):   (1) {ds1.res[0]:<10.8f} (2) {ds2.res[0]:.8f}")
            print(f"Resolution (Y):   (1) {ds1.res[1]:<10.8f} (2) {ds2.res[1]:.8f}")
            print(f"Width (pixels):   (1) {ds1.width:<5} (2) {ds2.width}")
            print(f"Height (pixels):  (1) {ds1.height:<5} (2) {ds2.height}")
            print(f"Bounds:           (1) {ds1.bounds}\n                  (2) {ds2.bounds}")
            # Safely format nodata and compression
            nodata1_str = str(ds1.nodata) if ds1.nodata is not None else 'None'
            nodata2_str = str(ds2.nodata) if ds2.nodata is not None else 'None'
            comp1_str = ds1.compression.value if ds1.compression else 'None'
            comp2_str = ds2.compression.value if ds2.compression else 'None'
            print(f"NoData Value:     (1) {nodata1_str:<10} (2) {nodata2_str}")
            print(f"Compression:      (1) {comp1_str:<10} (2) {comp2_str}")

            if ds1.count != ds2.count:
                 print("\n*** Band count differs! ***")
                 # List expected bands from TrendFire.py calculation
                 ls_indices = ['ndvi', 'evi', 'mirbi', 'ndfi', 'bsi', 'ndmi', 'nbr', 'nbr2', 'msavi', 'smi', 'ST_B10']
                 other_indices = ['rain', 'sm_surface', 'rh']
                 expected_bands = []
                 for index in ls_indices + other_indices:
                     expected_bands.extend([f"{index}_Slope", f"{index}_Intercept"])
                 
                 print(f"Expected bands based on TrendFire.py ({len(expected_bands)} bands total):")
                 # Print first few and last few expected bands for brevity
                 print(f"  {', '.join(expected_bands[:5])} ... {', '.join(expected_bands[-5:])}")
                 
                 # Try reading band descriptions if they exist
                 print("\nAttempting to read band descriptions...")
                 try:
                     print("\nBand descriptions (1) - inputResampled.tif:")
                     if ds1.descriptions:
                         for i, desc in enumerate(ds1.descriptions):
                            print(f"  Band {i+1}: {desc}")
                     else:
                         print("  No descriptions found.")

                     print("\nBand descriptions (2) - TrendFirePy_output.tif:")
                     if ds2.descriptions:
                         for i, desc in enumerate(ds2.descriptions):
                             print(f"  Band {i+1}: {desc}")
                     else:
                         print("  No descriptions found.")
                 except Exception as e:
                     print(f"Could not read band descriptions: {e}")

            # --- Data Comparison (Basic Stats for First Band) ---
            print("\n--- Basic Statistics (First Band) ---")
            try:
                data1 = ds1.read(1, masked=True)
                data2 = ds2.read(1, masked=True)
                
                # Check if arrays are comparable in shape
                if data1.shape == data2.shape:
                    print(f"Min Value:        (1) {np.ma.min(data1):<10.4f} (2) {np.ma.min(data2):.4f}")
                    print(f"Max Value:        (1) {np.ma.max(data1):<10.4f} (2) {np.ma.max(data2):.4f}")
                    print(f"Mean Value:       (1) {np.ma.mean(data1):<10.4f} (2) {np.ma.mean(data2):.4f}")
                    print(f"Std Dev:          (1) {np.ma.std(data1):<10.4f} (2) {np.ma.std(data2):.4f}")
                    
                    # Calculate difference
                    diff = data1 - data2
                    print(f"\nDifference (Band 1) Stats:")
                    print(f"  Min Diff:   {np.ma.min(diff):.4f}")
                    print(f"  Max Diff:   {np.ma.max(diff):.4f}")
                    print(f"  Mean Diff:  {np.ma.mean(diff):.4f}")
                    print(f"  Std Diff:   {np.ma.std(diff):.4f}")
                else:
                    print("Arrays have different shapes, cannot compare data directly.")

            except Exception as e:
                print(f"Error reading or comparing band data: {e}")

    except rasterio.RasterioIOError as e:
        print(f"Error opening GeoTIFF file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Define file paths relative to the script location or use absolute paths
    data_dir = "data"
    file1 = os.path.join(data_dir, "inputResampled.tif")  # Ground truth (JS output)
    file2 = os.path.join(data_dir, "TrendFirePy_output.tif") # Python output

    compare_tifs(file1, file2) 