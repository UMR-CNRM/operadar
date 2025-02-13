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
import operadar.operad_conf as cf



def save_netcdf(X:np.ndarray,
                 Y:np.ndarray,
                 Z:np.ndarray,
                 lat:np.ndarray,
                 lon:np.ndarray,
                 datetime:pd.Timestamp,
                 dpolDict:dict[np.ndarray],
                 contentsDict:dict[np.ndarray],
                 concentrationsDict:dict[np.ndarray],
                 temperature:np.ndarray,
                 outfile:Path,
                 ) :
    """Save synthetic dual-polarization variables and other model fields (Arome or MesoNH) in a netcdf file.

    Args:
        X (np.ndarray): 1D grid horizontal coordinates
        Y (np.ndarray): 1D grid horizontal coordinates
        Z (np.ndarray): 1D grid vertical coordinates
        lat (np.ndarray): 2D latitude coordinates
        lon (np.ndarray): 2D longitude coordinates
        datetime (pd.Timestamp): date and time coordinate
        dpolDict (dict[np.ndarray]): dictionnary containing the synthetic 3D radar fields
        contentsDict (dict[np.ndarray]): dictionnary containing the 3D model's content fields
        concentrationsDict (dict[np.ndarray]): dictionnary containing the 3D model's concentration fields
        temperature (np.ndarray): 3D model's temperature field
        outfile (Path): output path
    """

    # Contents and concentrations dict formatting for dataset backup
    hydromet_list = list(contentsDict.keys())
    contents = np.array([contentsDict[hydromet]*1000 for hydromet in hydromet_list]).astype('f4') # from kg/m3 to g/m3
    concentrations = np.array([concentrationsDict[hydromet] for hydromet in hydromet_list]).astype('f4')
    
    ds=xr.Dataset(
        data_vars=dict(Zh = (["level","y","x"],dpolDict["Zhh"].astype('f4'), {"units": "dBZ"}),
                       Zdr = (["level","y","x"],dpolDict["Zdr"].astype('f4'), {"units": "dB"}),
                       Kdp = (["level","y","x"],dpolDict["Kdp"].astype('f4'), {"units": "°/km"}),
                       Rhohv = (["level","y","x"],dpolDict["Rhohv"].astype('f4'), {"units": "1"}),
                       Contents = (["hydrometeor","level","y","x"],contents, {"units": "g/m3"}),
                       Concentrations = (["hydrometeor","level","y","x"],concentrations, {"units": "kg^-1"}),
                       T = (["level","y","x"],temperature.astype('f4'), {"units": "°C"}),
                       Alt = (["level","y","x"],Z.astype('i4'), {"units": "m"}),
                       ),
        coords=dict(y = (["y"], Y.astype('f4')),
                    x = (["x"], X.astype('f4')),
                    lon = (["y","x"], lon.astype('f4')),
                    lat = (["y","x"], lat.astype('f4')),
                    level =(["level"], np.arange(Z.shape[0]).astype('i4')),
                    hydrometeor = (["hydrometeor"],hydromet_list),
		            time = (datetime),
		            #Radloc = (["radpos"],Radpos),
                    ),
        attrs=dict(horizontal_resolution = int(Y[1]-Y[0]),
                   model = cf.model,
                   microphysics = cf.micro_scheme,
                   radar_band = cf.radar_band,
                   mixed_phase_type = cf.MixedPhase,
                   )
        )
     
    ds.to_netcdf(outfile.with_suffix('.nc'))
    ds.close() ; del ds
    print("Model and dpol variables saved at :",outfile.with_suffix('.nc'))



def create_tree_structure_outFiles(output_dir:Path):
    """Generate (if not existing yet) the tree structure to host output files."""
    if not output_dir.exists():
        try:
            output_dir.mkdir(exist_ok=True, parents=True)
            print ('Creating output tree structure :',output_dir)
        except:    
            print ('Error in creation of',output_dir) ; sys.exit()
    else:
        print ('Tree structure exists :',output_dir)   