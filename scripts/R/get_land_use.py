import ee
import pandas as pd

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize(project="ee-nikolasbasler")

# Load CSV with coordinates
csv_file = "data/sample_coordinates.csv"
df = pd.read_csv(csv_file)

# Load Copernicus Global Land Cover dataset (2019)
cropland_fraction = (
    ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")
    .filterDate("2019-01-01", "2019-12-31")
    .first()
    .select("crops-coverfraction")
)


# Function to get land cover classification for a point
def get_cropland_fraction(lat, lon):
    try:
        # Create a point geometry and buffer it to a 2 km radius.
        center = ee.Geometry.Point([lon, lat])
        region = center.buffer(2000)
        
        # Compute the mean cropland fraction over the circular region.
        # Update the band name below if necessary.
        result = cropland_fraction.reduceRegion(
            reducer=ee.Reducer.mean(), 
            geometry=region, 
            scale=100,
            maxPixels=1e9
        )
        
        # For debugging: print the entire result dictionary.
        result_dict = result.getInfo()
        print(f"Region reduction result at ({lat}, {lon}):", result_dict)
        
        # Attempt to get the cropland fraction value.
        # Update the key if the correct band name is different.
        fraction = result.get("crops-coverfraction")
        return fraction.getInfo() if fraction is not None else "No Data"
    except Exception as e:
        print(f"Error at ({lat}, {lon}): {e}")
        return "Error"

# Apply the function to each row
df["cropland_fraction"] = df.apply(
    lambda row: get_cropland_fraction(row["Latitude"], row["Longitude"]), axis=1
)

print(df)

# Save the results to a new CSV
output_file = "data/land_cover_results.csv"
df.to_csv(output_file, index=False)

# print(f"Results saved to {output_file}")
