#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import epygram
import pandas as pd
from pathlib import Path
from numpy import ndarray
from typing import Sequence
from netCDF4 import Dataset
from pandas import Timestamp



def format_temporal_variable(filePath:Path,
                             model_type:str,
                             real_case:bool,
                             )-> Timestamp | int:
    """Extract the temporal variable from an Arome or MesoNH file and format it if necessary."""

    if model_type=='Arome':
        epygram.init_env()
        epygram_file = epygram.formats.resource(filename=filePath, openmode = 'r', fmt = 'FA')
        date_time_file = epygram_file.validity.get()
        epygram_file.close()
        return Timestamp(date_time_file)
    elif model_type=='MesoNH' and real_case :
        epygram.init_env()
        epygram_file = epygram.formats.resource(filename=filePath, openmode = 'r',fmt='netCDFMNH') 
        field = epygram_file.readfield('ZWS')
        date_time_file = field.validity.get()
        epygram_file.close()
        return Timestamp(date_time_file)
    elif model_type=='MesoNH' and not real_case :
        mnh_file = Dataset(filePath,'r')
        date_time_file = mnh_file.variables['time'][:]
        return date_time_file[0]



def get_lat_lon_from_subdomain(domain:Sequence[float])-> tuple[float,float,float,float]:
    """Get latitude/longitude minimum and maximun from a list of
    float [lonmin, lonmax, latmin, latmax]"""
    lon_min = domain[0] ; lon_max = domain[1]
    lat_min = domain[2] ; lat_max = domain[3]
    
    return lon_min, lon_max, lat_min, lat_max



def Fw_or_Nc(momentsDict:dict[str,int],
             hydrometeor:str,
             concentration:ndarray,
             Fw:ndarray,
             micro:str,
             ):
    """ Select the right column in the lookup table depending of the number of moment for the hydrometeor.
    TODO : change the reading of the table so that the parameters are kept the same !

    * pristine ice (treated as 2-moment => the concentration can vary in MesoNH and is registered)
    * rainwater : 2-moment for all LIMA versions / 1-moment for ICE3
    * all other species (graupel / snow / wet graupel / wet snow ? )=> treated as 1-moment but can have a variable wet fraction
    => 3d parameter = Fw
    """
    # 2-moment species, and pristine ice whatever its number of moments: the tables of ice are computed as a
    # function of the concentration for all the schemes (ICE3, ICJW, LIMA, LIMC), see read/lookup_tables.py.
    # (Before Nov. 2026 only 'ICE3' and a misspelled 'IJW' used the concentration for 1-moment ice: ICJW and
    # LIMC ice were interpolated with Fw on the concentration axis.)
    if momentsDict[hydrometeor] == 2 or hydrometeor == 'ii':
        col_name = 'Nc'
        field_temp = concentration
    else:
        col_name = "Fw"
        field_temp = Fw
    return field_temp, col_name
        
        

def define_output_path(out_dir_path:str,
                       model:str,
                       scheme:str,
                       radar_band:str,
                       temporal_variable:Timestamp|int):
    """Define output path depending on the temporal variable type."""

    if type(temporal_variable) is pd.Timestamp :
        outName = f"dpolvar_{model}_{scheme}_{radar_band}band_{temporal_variable.strftime('%Y%m%d_%H%M')}"
    else :
        outName = f"dpolvar_{model}_{scheme}_{radar_band}band_{int(temporal_variable)}"
    return Path(os.path.join(out_dir_path,outName))
