#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""
import time as tm
import numpy as np
from pathlib import Path
from netCDF4 import Dataset
from typing import Sequence

from operadar.utils.make_links import link_keys_with_available_hydrometeors, link_varname_with_wrf_name
import operadar.read.with_netCDF4wrf as readncw


import wrf
from netCDF4 import Dataset
import pyproj
import xarray as xr
import os
import uuid
import shutil

def preprocess_wrf(filePath:Path):

    """
    Eulalia Busquets - 2/6/2026
    Preprocess for WRF outputs (.nc format) before being read by Operadar

    Necessary libraries: xarray, wrf, pyproj, netCDF4
    Notes:
    - The 'RHO' variable needs to be an output of WRF.
    The simpler option to get it is using the Runtime I/O option:
    http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.0/users_guide_chap5.html#runtimeio
    """


    print("Running WRF-Preprocessing")

    #Opening datasets
    ds = Dataset(filePath)
    xa =  xr.open_dataset(filePath)

    #Extract the variables needed by operadar and add them at the simulation:
    tc = wrf.getvar(ds, 'tc', wrf.ALL_TIMES)
    z = wrf.getvar(ds, 'z', wrf.ALL_TIMES)

    xa=xa.isel(Time=0) #To avoid mismatch between Dataset and Datarray time coordinates

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
    xa['bottom_top'] = xr.DataArray(ds['ZNU'][0], dims='bottom_top', attrs=eta_attrs)
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

    #output_file = Path("/tmp/preprocessed_WRF.nc") #Save in a temporary directory to avoid writting problems
    #                                               #This may be a problem. However, I belive that the solution is worse. 
     
    #Check if the folder modelFiles/WRF/tmp exist. If not, create the dir:
    dir_tmp = Path("modelFiles/WRF/tmp/") 
    dir_tmp.mkdir(parents=True, exist_ok=True)
    
    output_file = Path(f"modelFiles/WRF/tmp/preprocessed_WRF_{uuid.uuid4().hex}.nc")
    
    xa_cf.to_netcdf(output_file, mode='w', format='NETCDF4')
    return output_file



def read_wrfthom(filePath:Path,
                micro:str,
                subDomain:Sequence[int]|Sequence[float]|None,
                hydrometeorMoments: dict[str,int],real_case: bool,verbose:bool,
               )-> tuple[np.ndarray,
                         np.ndarray,
                         np.ndarray,
                         np.ndarray,
                         np.ndarray,
                         dict[str,np.ndarray],
                         dict[str,np.ndarray],
                         np.ndarray,
                         np.ndarray,
                         np.ndarray,
                         ]:
    """Read and extract data from an wrf.nc file
    !!! CAUTION! In this function, a previous .py has been run to obtain with wrf.getvar
        the Z and Tc variables. They are included in the .nc as 'z' and tc'
    Args:  
      filePath (str): input file path.
    Returns:
        X (ndarray): 1D horizontal coordinates in m
        Y (ndarray): 1D horizontal coordinates in m
        Z (ndarray): 3D array of altitude values for each model level
        lon (ndarray): 2D array of longitude coordinates
        lat (ndarray): 2D array of latitude coordinates
        M (dict[ndarray]): dictionary of 3D contents for each hydrometeor 
        Nc (dict[ndarray]): dictionary of 3D number concentrations for each hydrometeor
        Tc (ndarray) : 3D temperature in Celsius
        p (ndarray) : 3D pressure field in Pa
        qv (ndarray) : 3D specific humidity (kg/kg)
    """
    #Open file
    print("\tWRF .nc file:", )
    if verbose : deb=tm.time()
    wrf_file = Dataset(filePath,'r')

    if subDomain != None :
        i_min, i_max, j_min, j_max = readncw.get_subdomain_indices(wrfFile=wrf_file,
                                                            subDomain=subDomain,
                                                            real_case=real_case
                                                            )
        
    if verbose : print('\t\tSubdomain indices extracted.'); deb=tm.time()

    else:
        i_min, i_max, j_min, j_max = 0, -1, 0, -1
        if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
                                                        
      
    X, Y, Z, LAT, LON = readncw.get_geometry(wrfFile=wrf_file,
                                        real_case=real_case,
                                        i_min=i_min, i_max=i_max,
                                        j_min=j_min, j_max=j_max,
                                        )    
    if verbose : print('\t\tGot geometry in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    p, Tc, rho3D = readncw.get_pressure_temperature_density(wrfFile=wrf_file,
                                            i_min=i_min, i_max=i_max,
                                            j_min=j_min, j_max=j_max,
                                        )
    if verbose : print('\t\tGot pressure, temperature and dry air density in',round(tm.time()-deb,6),'seconds');deb=tm.time()


    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments, datatype='model')
    name_hydro = link_varname_with_wrf_name()

    M, qv = readncw.get_contents(wrfFile=wrf_file, hydrometeors=hydromet_list,
                        name_var_hydro=name_hydro, temperature=Tc, rho3D=rho3D,
                        i_min=i_min, i_max=i_max,
                        j_min=j_min, j_max=j_max,
                        )
    if verbose : print('\t\tGot 3D contents in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    Nc = readncw.get_concentrations(wrfFile=wrf_file, microphysics_scheme='ICE3',
                            hydrometeors=hydromet_list, temperature=Tc, rho3D=rho3D,
                            i_min=i_min, i_max=i_max,
                            j_min=j_min, j_max=j_max,
                            )
    if verbose : print('\t\tGot 3D number concentrations in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    
    return X, Y, Z, LON, LAT, M, Nc, Tc, p, qv

