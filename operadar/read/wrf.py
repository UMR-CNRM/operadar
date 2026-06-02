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

from operadar.utils.make_links import link_keys_with_available_hydrometeors, link_varname_with_mesonh_name, link_varname_with_wrf_name
import operadar.read.with_netCDF4wrf as readncw


import xarray as xr 

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
    print("\tWRF .nc file:",filePath)
    if verbose : deb=tm.time()
    mnh_file = Dataset(filePath,'r')

    if subDomain != None :
        i_min, i_max, j_min, j_max = readncw.get_subdomain_indices(mnhFile=mnh_file,
                                                            subDomain=subDomain,
                                                            real_case=real_case
                                                            )
        
    if verbose : print('\t\tSubdomain indices extracted.'); deb=tm.time()

    else:
        i_min, i_max, j_min, j_max = 0, -1, 0, -1
        if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
                                                        
      
    X, Y, Z, LAT, LON = readncw.get_geometry(mnhFile=mnh_file,
                                        real_case=real_case,
                                        i_min=i_min, i_max=i_max,
                                        j_min=j_min, j_max=j_max,
                                        )    
    if verbose : print('\t\tGot geometry in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    p, Tc, rho3D = readncw.get_pressure_temperature_density(mnhFile=mnh_file,
                                            i_min=i_min, i_max=i_max,
                                            j_min=j_min, j_max=j_max,
                                        )
    if verbose : print('\t\tGot pressure, temperature and dry air density in',round(tm.time()-deb,6),'seconds');deb=tm.time()


    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments, datatype='model')
    name_hydro = link_varname_with_wrf_name()

    M, qv = readncw.get_contents(mnhFile=mnh_file, hydrometeors=hydromet_list,
                        name_var_hydro=name_hydro, temperature=Tc, rho3D=rho3D,
                        i_min=i_min, i_max=i_max,
                        j_min=j_min, j_max=j_max,
                        )
    if verbose : print('\t\tGot 3D contents in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    Nc = readncw.get_concentrations(mnhFile=mnh_file, microphysics_scheme='ICE3',
                            hydrometeors=hydromet_list, temperature=Tc, rho3D=rho3D,
                            i_min=i_min, i_max=i_max,
                            j_min=j_min, j_max=j_max,
                            )
    if verbose : print('\t\tGot 3D number concentrations in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    
    return X, Y, Z, LON, LAT, M, Nc, Tc, p, qv
