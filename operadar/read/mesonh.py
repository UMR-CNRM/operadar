#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros & davidcl
"""

import time as tm
import numpy as np
from netCDF4 import Dataset
from operadar.utils.make_links import link_keys_with_available_hydrometeors, link_varname_with_mesonh_name
from operadar.read.with_netCDF4 import *



def read_mesonh(filePath: str,micro: str,subDomain:list[float]|None,
                hydrometeorMoments: dict[int],real_case: bool,verbose:bool,
               )-> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,dict[np.ndarray],dict[np.ndarray],np.ndarray]:
    """Read and extract data from an MesoNH.nc file

    Args:
        filePath (str): input file path.
        micro (str): microphysics scheme.
        subDomain (list[float] | None): either a list of 4 float or None.
        hydrometeorMoments (dict[int]): dictionary of form {'hydrometeor key':number of moments}.
        real_case (bool): if False, consider it is an idealized case.
        verbose (bool): show more messages to the user.

    Returns:
        X (ndarray): 1D horizontal coordinates in m
        Y (ndarray): 1D horizontal coordinates in m
        Z (ndarray): 3D array of altitude values for each model level
        LON (ndarray): 2D array of longitude coordinates
        LAT (ndarray): 2D array of latitude coordinates
        M (dict[ndarray]): dictionary of 3D contents for each hydrometeor 
        Nc (dict[ndarray]): dictionary of 3D number concentrations for each hydrometeor
        Tc (ndarray) : 3D temperature in Celsius
    """

    print("\tMesoNH .nc file:",filePath)
    if verbose : deb=tm.time()
    mnh_file = Dataset(filePath,'r')
    if verbose : print('\t\tLoaded file in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    check_variable_is_in_dataset(mnh_file)
    
    if subDomain != None :
        i_min, i_max, j_min, j_max = get_subdomain_indices(mnhFile=mnh_file,
                                                           subDomain=subDomain,
                                                           real_case=real_case,
                                                           )
        if verbose : print('\t\tSubdomain indices extracted.'); deb=tm.time()
    else:
        i_min, i_max, j_min, j_max = 0, -1, 0, -1
        if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()

    X, Y, Z, LAT, LON = get_geometry(mnhFile=mnh_file,
                                     real_case=real_case,
                                     i_min=i_min, i_max=i_max,
                                     j_min=j_min, j_max=j_max,
                                     )    
    if verbose : print('\t\tGot geometry in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    Tc, rho3D = get_temperature_and_density(mnhFile=mnh_file,
                                            i_min=i_min, i_max=i_max,
                                            j_min=j_min, j_max=j_max,
                                            )
    if verbose : print('\t\tGot temperature and dry air density in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments, datatype='model')
    name_hydro = link_varname_with_mesonh_name()
    
    M = get_contents(mnhFile=mnh_file, hydrometeors=hydromet_list,
                     name_var_hydro=name_hydro, temperature=Tc, rho3D=rho3D,
                     i_min=i_min, i_max=i_max,
                     j_min=j_min, j_max=j_max,
                     )
    if verbose : print('\t\tGot 3D contents in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    Nc = get_concentrations(mnhFile=mnh_file, microphysics_scheme=micro,
                            hydrometeors=hydromet_list, temperature=Tc, rho3D=rho3D,
                            i_min=i_min, i_max=i_max,
                            j_min=j_min, j_max=j_max,
                            )
    if verbose : print('\t\tGot 3D number concentrations in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    return X, Y, Z, LON, LAT, M, Nc, Tc
