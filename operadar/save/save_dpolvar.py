#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path



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
                 dpol2add:list,
                 model:str,
                 micro_scheme:str,
                 radar_band:str,
                 MixedPhase:str,
                 outfile:Path,
                 radarloc:str|list,
                 scattering_method:str,
                 hydrometeors_moments:dict,
                 real_case:bool,
                 ) :
    """Save synthetic dual-polarization variables and other model fields (Arome or MesoNH) in a netcdf file.

    Args:
        X (np.ndarray): 1D grid horizontal coordinates
        Y (np.ndarray): 1D grid horizontal coordinates
        Z (np.ndarray): 3D altitude coordinates
        lat (np.ndarray): 2D latitude coordinates
        lon (np.ndarray): 2D longitude coordinates
        datetime (pd.Timestamp): date and time coordinate
        dpolDict (dict[np.ndarray]): dictionary containing the synthetic 3D radar fields
        contentsDict (dict[np.ndarray]): dictionary containing the 3D model's content fields
        concentrationsDict (dict[np.ndarray]): dictionary containing the 3D model's concentration fields
        temperature (np.ndarray): 3D model's temperature field
        outfile (Path): output path
    """

    # Contents and concentrations dict formatting for dataset backup
    hydromet_list = list(contentsDict.keys())
    contents = np.array([contentsDict[hydromet]*1000 for hydromet in hydromet_list]).astype('f4') # from kg/m3 to g/m3
    concentrations = np.array([concentrationsDict[hydromet] for hydromet in hydromet_list]).astype('f4')
    dataset_variables = dict(Contents = (["hydrometeor","level","y","x"],contents, {"units": "g/m3"}),
                             Concentrations = (["hydrometeor","level","y","x"],concentrations, {"units": "kg^-1"}),
                             T = (["level","y","x"],temperature.astype('f4'), {"units": "°C"}),
                             Alt = (["level","y","x"],Z.astype('i4'), {"units": "m"}),
                             )
    dataset_variables = add_dualPol_variables(ds_variables=dataset_variables,
                                              dpolDict=dpolDict,
                                              dpolvar2add=dpol2add,
                                              ) 
    dataset_coordinates = dict(y = (["y"], Y.astype('i4')),
                               x = (["x"], X.astype('i4')),
                               level =(["level"], np.arange(Z.shape[0]).astype('i4')),
                               hydrometeor = (["hydrometeor"],hydromet_list),
                               time = (datetime),
                               )
    dataset_coordinates = add_lat_lon_coordinates(ds_coordinates=dataset_coordinates,
                                                  lat=lat,
                                                  lon=lon,
                                                  real_case=real_case,
                                                  )
    dataset_attributes = dict(horizontal_resolution = int(Y[1]-Y[0]),
                              model = model,
                              microphysics = micro_scheme,
                              radar_band = radar_band,
                              mixed_phase_type = MixedPhase,
                              radar_location = str(radarloc),
                              scattering_method = scattering_method,
                              )
    dataset_attributes.update({f'{key}_moment':value for key,value in hydrometeors_moments.items()})
        
    ds=xr.Dataset(data_vars = dataset_variables,
                  coords = dataset_coordinates,
                  attrs=dataset_attributes,
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



def add_dualPol_variables(ds_variables:dict,dpolDict:dict,dpolvar2add:list):
    units = {'Zh'    : {"units": "dBZ"  },
             'Zdr'   : {"units": "dB"   },
             'Kdp'   : {"units": "°/km" },
             'Rhohv' : {"units": "1"    },
             'Ah'    : {"units": "dB/km"},
             'Av'    : {"units": "dB/km"},
             }
    for var in dpolvar2add :
        ds_variables[var] = (["level","y","x"],dpolDict[var].astype('f4'), units[var])
    
    return ds_variables



def add_lat_lon_coordinates(ds_coordinates:dict,lat:np.ndarray,lon:np.ndarray,real_case:bool):
    if real_case :
        ds_coordinates['lon'] = (["y","x"], lon.astype('f4'))
        ds_coordinates['lat'] = (["y","x"], lat.astype('f4'))
    return ds_coordinates
