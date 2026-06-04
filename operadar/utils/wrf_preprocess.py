"""
Eulalia Busquets - 2/6/2026
Preprocess for WRF outputs (.nc format) before being read by Operadar

Necessary libraries: xarray, wrf, pyproj, netCDF4
Notes:
- The 'RHO' variable needs to be an output of WRF.
  The simpler option to get it is using the Runtime I/O option:
  http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.0/users_guide_chap5.html#runtimeio
- A single time output needs to be selected
"""

import wrf
from netCDF4 import Dataset
import pyproj

import numpy as np
import xarray as xr

#-----------------------------------------------
#INTRODUCE HERE THE PATH OF THE SIMULATION, THE SELECTED TIME COORDINATE AND THE OUTPUT NAME
path = "/cnrm/precip/users/busquets/data/wrfout_d03_20241025_08" #Path to the file 
t = 131 #Time that will be saved
output_name = "wrfout_d03_20241025_08_processed_t131.nc" #Name of the processed simulation
#-----------------------------------------------
print(f"Preprocessing: \
      file: {path}\
      time_output:{t}\
      output name: {output_name}")
print('Opening datasets...')
#Opening datasets
ds = Dataset(path)
xa =  xr.open_dataset(path)

#Extract the variables needed by operadar and add them at the simulation:
print('computing tc...')
tc = wrf.getvar(ds, 'tc', wrf.ALL_TIMES)
print('computing z...')
z = wrf.getvar(ds, 'z', wrf.ALL_TIMES)
print('Adding variables...')
xa['tc'] = tc
xa['z'] = z

#-----------------------------------------------
# Code extracted from Jthielen "wrf_python_extraction" 
# https://gist.github.com/jthielen/8881d32d08d625c75c906a7e8ad7583f
#-----------------------------------------------

print('opening WRF projection')
# Define the WRF projection (see https://fabienmaussion.info/2018/01/06/wrf-projection/)
wrf_proj = pyproj.Proj(proj='lcc',
                       lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2,
                       lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON,
                       a=6370000, b=6370000)
# Easting and Northing of the domain center point
wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
e, n = pyproj.transform(wgs_proj, wrf_proj, ds.CEN_LON, ds.CEN_LAT)

# Grid parameters
dx, dy = ds.DX, ds.DY
nx, ny = ds.dimensions['west_east'].size, ds.dimensions['south_north'].size

# Lower left corner of the domain
x0 = -(nx-1) / 2. * dx + e
y0 = -(ny-1) / 2. * dy + n

# Get grid values
x, y = np.arange(nx) * dx + x0, np.arange(ny) * dy + y0

# Add in dimension coordinates
eta_attrs = {attr: ds['ZNU'].getncattr(attr) for attr in ds['ZNU'].ncattrs()}
eta_attrs['axis'] = 'Z'
xa['bottom_top'] = xr.DataArray(ds['ZNU'][17], dims='bottom_top', attrs=eta_attrs)
xa['south_north'] = xr.DataArray(y, dims='south_north', attrs={'axis': 'Y', 'units': 'm'})
xa['west_east'] = xr.DataArray(x, dims='west_east', attrs={'axis': 'X', 'units': 'm'})

# Define the grid_mapping
xa['LambertConformal'] = xr.DataArray(np.array(0), attrs={
    'grid_mapping_name': 'lambert_conformal_conic',
    'earth_radius': 6370000,
    'standard_parallel': (ds.TRUELAT1, ds.TRUELAT2),
    'longitude_of_central_meridian': ds.STAND_LON,
    'latitude_of_projection_origin': ds.MOAD_CEN_LAT
})

for var in xa.data_vars:
    xa[var].attrs['grid_mapping'] = 'LambertConformal'
    if 'projection' in xa[var].attrs:
        del xa[var].attrs['projection']
    if 'coordinates' in xa[var].attrs:
        del xa[var].attrs['coordinates']

#-----------------------------------------------
print('renaming time...')
xa_cf = xa.rename({'Time':'time'}) #I need to rename Time, to follow the CF1.6-compliant temporal unit

xa_t=xa_cf.isel(time=t) #I need to pick just one time

xa_t.to_netcdf(output_name)