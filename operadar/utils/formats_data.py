#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import epygram
import pandas as pd
from pathlib import Path
from numpy import ndarray
from netCDF4 import Dataset
from pandas import Timestamp



def format_temporal_variable(filePath:Path,model_type:str,real_case:bool)-> Timestamp:
    """Extract the temporal variable from an Arome or MesoNH file and format it if necessary."""

    if model_type=='Arome':
        epygram.init_env()
        epygram_file = epygram.formats.resource(filename=filePath, openmode = 'r', fmt = 'FA')
        date_time_file = epygram_file.validity.get()
        epygram_file.close()
        return Timestamp(date_time_file)
    elif model_type=='MesoNH' and real_case :
        epygram_file = epygram.formats.resource(filename=filePath, openmode = 'r',fmt='netCDFMNH') 
        field = epygram_file.readfield('ZWS')
        date_time_file = field.validity.get()
        epygram_file.close()
        return Timestamp(date_time_file)
    elif model_type=='MesoNH' and not real_case :
        mnh_file = Dataset(filePath,'r')
        date_time_file = mnh_file.variables['time'][:]
        return date_time_file[0]



def get_lat_lon_from_subdomain(domain:list[float])-> tuple[float,float,float,float]:
    """Get latitude/longitude minimum and maximun from a list of
    float [lonmin, lonmax, latmin, latmax]"""
    lon_min = domain[0] ; lon_max = domain[1]
    lat_min = domain[2] ; lat_max = domain[3]
    
    return lon_min, lon_max, lat_min, lat_max



def select_Tmatrix_column(momentsDict:dict[int],
                          hydrometeor:str,
                          concentration:ndarray,
                          Fw:ndarray,
                          tmatrix_param:dict):
    """ Compute 3d parameter in Tmatrix table 
    TODO : change the reading of Tmatrix table so that the parameters are kept the same !

    * pristine ice (treated as 2-moment => the concentration can vary in MesoNH and is registered)
    * rainwater : 2-moment for all LIMA versions / 1-moment for ICE3
    * all other species (graupel / snow / wet graupel / wet snow ? )=> treated as 1-moment but can have a variable wet fraction
    => 3d parameter = Fw
    """
    
    if momentsDict[hydrometeor] == 2 :
        col_name = 'Nc'
        field_temp = concentration
        col = 'expCC'
    else:
        col_name = "Fw"
        field_temp = Fw
        col = 'Fw'
    col_min = tmatrix_param[f'{col}min'][hydrometeor]
    col_step = tmatrix_param[f'{col}step'][hydrometeor]
    col_max = tmatrix_param[f'{col}max'][hydrometeor]
    return field_temp, col_min, col_max, col_step, col_name

        
        

def define_output_path(out_dir_path,model,scheme,radar_band,temporal_variable):
    """Define output path depending on the temporal variable type."""
    if type(temporal_variable) is pd.Timestamp :
        outPath = f"{out_dir_path}dpolvar_{model}_{scheme}_{radar_band}band_{temporal_variable.strftime('%Y%m%d_%H%M')}"
    elif type(temporal_variable)==int:
        outPath = f"{out_dir_path}dpolvar_{model}_{scheme}_{radar_band}band_{temporal_variable}"
    return outPath
