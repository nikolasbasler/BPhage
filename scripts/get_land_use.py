import ee
import pandas as pd

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize()

# Load CSV with coordinates
csv_file = "data/metadata/coordinates.csv"
df = pd.read_csv(csv_file)

# Load Copernicus Global Land Cover dataset (2019)
land_cover = (
    ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")
    .filterDate("2019-01-01", "2019-12-31")
    .first()
    .select("discrete_classification")
)


# Function to get land cover classification for a point
def get_land_cover(lat, lon):
    try:
        point = ee.Geometry.Point([lon, lat])
        land_class = land_cover.reduceRegion(
            reducer=ee.Reducer.first(), geometry=point, scale=100
        ).get("discrete_classification")

        # Fetch the result from Earth Engine
        return land_class.getInfo() if land_class else "No Data"
    except Exception as e:
        print(f"Error at ({lat}, {lon}): {e}")
        return "Error"


# Apply the function to each row
df["land_cover_class"] = df.apply(
    lambda row: get_land_cover(row["latitude"], row["longitude"]), axis=1
)

# Save the results to a new CSV
output_file = "data/metadata/land_cover_results.csv"
df.to_csv(output_file, index=False)

print(f"Results saved to {output_file}")
