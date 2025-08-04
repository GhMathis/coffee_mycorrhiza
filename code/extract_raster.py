import os
import gzip
import tempfile
import numpy as np
import shutil
import urllib
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import matplotlib.colors as colors
from datetime import datetime

import netCDF4
from azure.storage.blob import ContainerClient
import cartopy.crs as ccrs

# Set up temp directory
temp_dir = os.path.join(tempfile.gettempdir(), 'ghe')
os.makedirs(temp_dir, exist_ok=True)

# Azure storage constants
storage_account_name = 'ghe'
container_name = 'noaa-ghe'
storage_account_url = f'https://{storage_account_name}.blob.core.windows.net'
ghe_blob_root = f'{storage_account_url}/{container_name}/'

ghe_container_client = ContainerClient(
    account_url=storage_account_url,
    container_name=container_name,
    credential=None
)

# Download URL helper
def download_url(url):
    url_as_filename = url.replace('://', '_').replace('/', '_')
    destination_filename = os.path.join(temp_dir, url_as_filename)
    if not os.path.exists(destination_filename):
        urllib.request.urlretrieve(url, destination_filename)
    return destination_filename

# Load grid
grid_file_url = 'https://ghe.blob.core.windows.net/noaa-ghe/NPR.GEO.GHE.v1.Navigation.netcdf.gz'
grid_filename_gz = download_url(grid_file_url)
with gzip.open(grid_filename_gz) as gz:
    grid_dataset = netCDF4.Dataset('dummy', mode='r', memory=gz.read())
    lat_grid_raw = grid_dataset['latitude'][:]
    lon_grid_raw = grid_dataset['longitude'][:]
    grid_dataset.close()

# Function to process a single day
def process_day(syear, smonth, sday):
    product = 'rain_rate'
    prefix = f'{product}/{syear}/{smonth}/{sday}/'
    generator = ghe_container_client.list_blobs(name_starts_with=prefix)
    blob_urls = [f'{ghe_blob_root}{blob.name}' for blob in generator]

    if not blob_urls:
        print(f"No data found for {syear}-{smonth}-{sday}")
        return None

    rainfall = np.zeros(lat_grid_raw.shape)
    for blob_url in tqdm(blob_urls, desc=f'Processing {syear}-{smonth}-{sday}'):
        try:
            filename = download_url(blob_url)
            with gzip.open(filename) as gz:
                dataset = netCDF4.Dataset('dummy', mode='r', memory=gz.read())
                rainfall_sample = dataset['rain'][:]
                rainfall_sample[rainfall_sample < 0] = 0
                rainfall += rainfall_sample
                dataset.close()
        except Exception as e:
            print(f"Error processing {blob_url}: {e}")
            continue

    return rainfall

# Iterate over a date range
start_date = datetime.strptime('2022-05-31', '%Y-%m-%d')
end_date = datetime.strptime('2023-06-01', '%Y-%m-%d')  # exclusive

rainfall_dict = {}

current_date = start_date
while current_date < end_date:
    syear = current_date.strftime('%Y')
    smonth = current_date.strftime('%m')
    sday = current_date.strftime('%d')
    rainfall = process_day(syear, smonth, sday)
    if rainfall is not None:
        rainfall_dict[current_date.strftime('%Y-%m-%d')] = rainfall
    current_date += datetime.timedelta(days=1)

# Example: print min/max for a date
for date_str, rainfall in rainfall_dict.items():
    print(f"{date_str}: min={np.min(rainfall)}, max={np.max(rainfall)}")


# Initialize the total rainfall raster with zeros, same shape as the first raster
all_dates = list(rainfall_dict.keys())
if not all_dates:
    raise ValueError("rainfall_dict is empty — nothing to sum!")

raster_shape = rainfall_dict[all_dates[0]].shape
total_rainfall = np.zeros(raster_shape)

# Sum all daily rasters
for date in all_dates:
    total_rainfall += rainfall_dict[date]

# Print stats
print(f"Total rainfall raster:")
print(f"  Shape: {total_rainfall.shape}")
print(f"  Min: {np.min(total_rainfall)}")
print(f"  Max: {np.max(total_rainfall)}")

# Parallel execution
# rainfall = np.zeros(lat_grid_raw.shape)
# n_valid = np.zeros(lat_grid_raw.shape)

# rain_units = None
# variable_description = None

# with ThreadPoolExecutor() as executor:
#     results = list(tqdm(executor.map(process_blob, blob_urls), total=len(blob_urls)))

# # Accumulate results
# for rainfall_sample, units, description in results:
#     rainfall += rainfall_sample
#     rain_units = units
#     variable_description = description  # Overwritten, but should be same for all

# min_rf = np.min(rainfall)
# max_rf = np.max(rainfall)

# print('Rainfall ranges from {}{} to {}{}'.format(min_rf, rain_units, max_rf, rain_units))

# rainfall_raw = rainfall.copy()

import rasterio
from rasterio.transform import from_bounds
# Get the dimensions of your data
height, width = total_rainfall.shape

# Get the extent of your data from the grid files
# Assuming the grid_dataset variables define the corners of the grid
# For accurate georeferencing, you'd ideally use the actual corner coordinates
# or a more robust method to derive the transform.
# Here, we'll approximate using the min/max of your grid arrays.
min_lon = np.nanmin(lon_grid_raw)
max_lon = np.nanmax(lon_grid_raw)
min_lat = np.nanmin(lat_grid_raw)
max_lat = np.nanmax(lat_grid_raw)

# Create a geotransform: This defines how the image coordinates
# map to geographic coordinates.
# from_bounds(west, south, east, north, width, height)
# The order of arguments is crucial.
transform = from_bounds(west=min_lon, south=min_lat, east=max_lon, north=max_lat,
                        width=width, height=height)

# Define the coordinate reference system (CRS).
# WGS84 (latitude/longitude) is a common choice for global data.
# Its EPSG code is 4326.
crs = 'EPSG:4326'
output_raster_filename = "data/raster/rainfall/rainfall_global_year2022_2023.tif"
# Write the raster file
with rasterio.open(
    output_raster_filename,
    'w',
    driver='GTiff',
    height=height,
    width=width,
    count=1,  # Number of bands (e.g., 1 for single-band rainfall data)
    dtype=total_rainfall.dtype, # Use the data type of your NumPy array
    crs=crs,
    transform=transform,
    nodata=-9999.0 # Set the nodata value if you have one, or comment out
) as dst:
    dst.write(total_rainfall, 1) # Write the rainfall_raw array to band 1

print(f"\nRaster data saved to: {output_raster_filename}")

# Define the projection
proj = ccrs.PlateCarree()  # or use ccrs.PlateCarree() if more appropriate

# Set up figure and axis
fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': proj})

# Set extent: (lon_min, lon_max, lat_min, lat_max)
ax.set_extent([
    np.nanmin(lon_grid_raw), np.nanmax(lon_grid_raw),
    np.nanmin(lat_grid_raw), np.nanmax(lat_grid_raw)
], crs=proj)


# Plot the data
total_rainfall += 1
norm = colors.LogNorm(vmin=total_rainfall.min(), vmax=total_rainfall.max(), clip=True)
cmap = plt.cm.Blues

# Plot the color mesh (assuming lon_grid and lat_grid are 2D)
mesh = ax.pcolormesh(
    lon_grid_raw, lat_grid_raw, total_rainfall,
    transform=proj,  # specify that the data is in lat/lon
    cmap=cmap, norm=norm, shading='auto'
)

# Add coastlines, borders, gridlines
ax.coastlines(resolution='110m')
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
gl.top_labels = False
gl.right_labels = False

# Add colorbar
cb = plt.colorbar(mesh, ax=ax, orientation='vertical', pad=0.05)
cb.set_label(rain_units)

# Title
plt.title(f'Global rainfall ({rain_units})')
plt.tight_layout()
plt.show()
