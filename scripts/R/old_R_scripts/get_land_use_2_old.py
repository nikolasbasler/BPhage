import ee
import pandas as pd

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize(project="ee-nikolasbasler")

# Load CSV with coordinates
csv_file = "data/sample_coordinates.csv"
df = pd.read_csv(csv_file)

# Load Copernicus Global Land Cover dataset (2019) and select the discrete classification band.
land_cover = (
    ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")
    .filterDate("2019-01-01", "2019-12-31")
    .first()
    .select("discrete_classification")
)

# Define the cropland class code.
cropland_code = 40

# Create a binary mask: pixels classified as cropland become 1, others 0.
cropland_mask = land_cover.eq(cropland_code)

# Multiply the mask by the pixel area (in m²) to get an image where each pixel's value is its cropland area.
# Rename the band so it’s easier to retrieve later.
cropland_area_image = cropland_mask.multiply(ee.Image.pixelArea()).rename("cropland_area")

# Function to get total cropland area for a point
def get_cropland_area(lat, lon):
    try:
        # Create a point geometry and buffer it to a 2 km radius.
        center = ee.Geometry.Point([lon, lat])
        region = center.buffer(2000)
        
        # Sum the cropland area (in m²) over the region.
        result = cropland_area_image.reduceRegion(
            reducer=ee.Reducer.sum(),
            geometry=region,
            scale=100,    # Adjust if needed to match the dataset's resolution
            maxPixels=1e9
        )
        
        # For debugging: print the entire result dictionary.
        result_dict = result.getInfo()
        print(f"Region reduction result at ({lat}, {lon}):", result_dict)
        
        # Retrieve the computed area using the renamed band.
        area = result.get("cropland_area")
        return area.getInfo() if area is not None else "No Data"
    except Exception as e:
        print(f"Error at ({lat}, {lon}): {e}")
        return "Error"

# Apply the function to each row of the dataframe
df["cropland_area"] = df.apply(
    lambda row: get_cropland_area(row["Latitude"], row["Longitude"]), axis=1
)

print(df)
