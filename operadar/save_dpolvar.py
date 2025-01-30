#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:00:34 2023

@author: augros
"""

import sys
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path



def save_dpolvar(M:dict[np.ndarray], Nc:dict[np.ndarray], Vm_k:dict[np.ndarray], Tc:np.ndarray,
                 Z:np.ndarray, X:np.ndarray, Y:np.ndarray,
                 lat:np.ndarray, lon:np.ndarray,
                 datetime:pd.Timestamp, outfile:Path):
    """
    Save synthetic dual-polarization variables and other model fields 
    (Arome or MesoNH) in a netcdf file.
    
    INPUTS
    - X,Y : 1D horizontal grid coordinates
    - Z : model levels indices (1D)
    - lat, lon : 2D fields with dimension = to (X,Y)
    - datetime : date and time information
    - M : dictionary of array for hydrometeor contents (conversion from kg to g/m3)
    - Nc : dictionary of array for hydrometeor number concentrations
    - Tc : temperature field
    - Vm_k : dictionary of array for all dual-pol variables
    - outfile
    """
    # M and Nc dict formatting for dataset backup
    hydromet_list = list(M.keys())
    contents = np.array([M[hydromet]*1000 for hydromet in hydromet_list]).astype('f4') # from kg to g/m3
    concentrations = np.array([Nc[hydromet] for hydromet in hydromet_list]).astype('f4')
    
    ds=xr.Dataset(
        data_vars=dict(Zh     = (["level","y","x"],Vm_k["Zhh"].astype('f4'), {"units": "dBZ"}),
                       Zdr    = (["level","y","x"],Vm_k["Zdr"].astype('f4'), {"units": "dB"}),
                       Kdp    = (["level","y","x"],Vm_k["Kdp"].astype('f4'), {"units": "°/km"}),
                       Rhohv  = (["level","y","x"],Vm_k["Rhohv"].astype('f4'), {"units": "1"}),
                       M      = (["hydrometeor","level","y","x"],contents, {"units": "g/m3"}),
                       Nc     = (["hydrometeor","level","y","x"],concentrations, {"units": "kg^-1"}),
                       T      = (["level","y","x"],Tc.astype('f4'), {"units": "°C"}),
                       Alt    = (["level","y","x"],Z.astype('i4'), {"units": "m"}),
                       ),
        coords=dict(y   = (["y"], Y.astype('f4')),
                    x   = (["x"], X.astype('f4')),
                    lon = (["y","x"], lon.astype('f4')),
                    lat = (["y","x"], lat.astype('f4')),
                    level=(["level"], np.arange(Z.shape[0]).astype('i4')),
                    hydrometeor = (["hydrometeor"],hydromet_list),
		            time = (datetime),
		            #Radloc = (["radpos"],Radpos),
                    ),
        )
     
    ds.to_netcdf(outfile.with_suffix('.nc'))
    ds.close() ; del ds
    print("Model and dpol variables saved at :",outfile)



def create_tree_structure_outFiles(output_dir:Path):
    if not output_dir.exists():
        try:
            output_dir.mkdir(exist_ok=True, parents=True)
            print ('Creating output directories :',output_dir)
        except:    
            print ('Error in creation of',output_dir) ; sys.exit()
    else:
        print ('Output directories exist :',output_dir)   