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


import netCDF4
from azure.storage.blob import ContainerClient
import cartopy.crs as ccrs

# The grid spacing for all GHE files is defined in a separate NetCDF file.
grid_file_url = 'https://ghe.blob.core.windows.net/noaa-ghe/NPR.GEO.GHE.v1.Navigation.netcdf.gz'

temp_dir = os.path.join(tempfile.gettempdir(),'ghe')
os.makedirs(temp_dir,exist_ok=True)
# Azure storage constants
storage_account_name = 'ghe'
container_name = 'noaa-ghe'
storage_account_url = 'https://' + storage_account_name + '.blob.core.windows.net'
ghe_blob_root = storage_account_url + '/' + container_name + '/'

ghe_container_client = ContainerClient(account_url=storage_account_url, 
                                         container_name=container_name,
                                         credential=None)
# Functions
# Download a URL to a temporary file
def download_url(url):
    url_as_filename = url.replace('://', '_').replace('/', '_')    
    destination_filename = os.path.join(temp_dir,url_as_filename)    
    urllib.request.urlretrieve(url, destination_filename)  
    return destination_filename
# Download the grid spacing file
# This file is ~150MB, so best to cache this
grid_filename_gz = download_url(grid_file_url)
with gzip.open(grid_filename_gz) as gz:
    grid_dataset = netCDF4.Dataset('dummy', mode='r', memory=gz.read())
    print(grid_dataset.variables)
    lat_grid_raw = grid_dataset['latitude']
    lon_grid_raw = grid_dataset['longitude']

# Select data
# Data are stored as product/year/month/day/filename
product = 'rain_rate'

# Grab data from April 9, 2020
syear = '2023'; smonth = '04'; sday = '09'

# Filenames look like:
#
# NPR.GEO.GHE.v1.S202001170000.nc.gz
#
# ...where the last four digits represent time, n increments of 15 minutes from 0000

# We can either sum over a whole day, or take a single 15-minute window
single_time_point = False

if single_time_point:
    
    # Pick an arbitrary time of day to plot
    stime = '0200'
    
    filename = 'NPR.GEO.GHE.v1.S' + syear + smonth + sday + stime + '.nc.gz'
    blob_urls = [ghe_blob_root + product + '/' + syear + '/' + smonth + '/' + sday + '/' \
                 + filename]
    
else:
    
    prefix = product + '/' + syear  + '/' #+ smonth + '/' + sday + '/'
    print('Finding blobs matching prefix: {}'.format(prefix))
    generator = ghe_container_client.list_blobs(name_starts_with=prefix)
    blob_urls = []
    for blob in generator:
        blob_urls.append(ghe_blob_root + blob.name)
    print('Found {} matching scans'.format(len(blob_urls)))

####
def process_blob(blob_url):
    filename = download_url(blob_url)
    with gzip.open(filename) as gz:
        dataset = netCDF4.Dataset('dummy', mode='r', memory=gz.read())
        rainfall_sample = dataset['rain'][:]
        rainfall_sample[rainfall_sample < 0] = 0  # remove invalid
        rain_units = dataset['rain'].units
        variable_description = str(dataset.variables)
        dataset.close()
    return rainfall_sample, rain_units, variable_description

# Parallel execution
rainfall = np.zeros(lat_grid_raw.shape)
n_valid = np.zeros(lat_grid_raw.shape)

rain_units = None
variable_description = None

with ThreadPoolExecutor() as executor:
    results = list(tqdm(executor.map(process_blob, blob_urls), total=len(blob_urls)))

# Accumulate results
for rainfall_sample, units, description in results:
    rainfall += rainfall_sample
    rain_units = units
    variable_description = description  # Overwritten, but should be same for all

min_rf = np.min(rainfall)
max_rf = np.max(rainfall)

print('Rainfall ranges from {}{} to {}{}'.format(min_rf, rain_units, max_rf, rain_units))

rainfall_raw = rainfall.copy()
###

# Take a look at what's in each NetCDF file
print(variable_description)

# Prepare indices, downsample for faster plotting
image_size = np.shape(rainfall_raw)
nlat = image_size[0]; nlon = image_size[1]

assert(np.shape(rainfall_raw)==np.shape(lat_grid_raw))
assert(np.shape(rainfall_raw)==np.shape(lon_grid_raw))

# Downsample by decimation
ds_factor = 10

lon_grid = lon_grid_raw[::ds_factor,::ds_factor,]
lat_grid = lat_grid_raw[::ds_factor,::ds_factor,]
rainfall = rainfall_raw[::ds_factor,::ds_factor,]
print(lon_grid)

#Plot rainfall
plt.figure(figsize=(20,20))

# Prepare a matplotlib Basemap so we can render coastlines and borders
m = Basemap(projection='merc',
  llcrnrlon=np.nanmin(lon_grid),urcrnrlon=np.nanmax(lon_grid),
  llcrnrlat=np.nanmin(lat_grid),urcrnrlat=np.nanmax(lat_grid),
  resolution='c')

# Convert lat/lon to a 2D grid
# lon_grid,lat_grid = np.meshgrid(lon,lat)
x,y = m(lon_grid,lat_grid)

# Clip our plot values to an upper threshold, and leave anything
# below the lower threshold as white (i.e., unplotted)
n_files = len(blob_urls)
upper_plot_threshold = n_files*10
lower_plot_threshold = n_files*0.01

Z = rainfall.copy()
Z[Z > upper_plot_threshold] = upper_plot_threshold
Z[Z < lower_plot_threshold] = np.nan
Z = np.ma.masked_where(np.isnan(Z),Z)

# Choose normalization and color mapping
norm = mpl.colors.LogNorm(vmin=Z.min(), vmax=Z.max(), clip=True)
cmap = plt.cm.Blues

# Plot as a color mesh
cs = m.pcolormesh(x,y,Z,norm=norm,cmap=cmap,shading='auto')

# Draw extra stuff to make our plot look fancier... sweeping clouds on a plain background
# are great, but sweeping clouds on contentinal outlines are *very* satisfying.
m.drawcoastlines()
m.drawmapboundary()
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
m.colorbar(cs)

plt.title('Global rainfall ({})'.format(rain_units))
plt.show()