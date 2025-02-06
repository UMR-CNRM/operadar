#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pandas import Timestamp
from numpy import ndarray



# A SUPPRIMER PROCHAINEMENT
def get_vortex_experiments(csvRow: str, microphysics_scheme: str) :
    listExpe = [str(x) for x in csvRow.expeNames.strip().split(',')]
    if microphysics_scheme == 'ICE3' :
        return listExpe[0]
    elif microphysics_scheme == 'ICE4' :
        return listExpe[1]
    elif microphysics_scheme == 'LIMASG' :
        return listExpe[2]
    elif microphysics_scheme == 'LIMAAG' :
        return listExpe[3]
    elif microphysics_scheme == 'LIMA49t' :
        return listExpe[4]



def format_date_time_argument(date_time:str|Timestamp):
    """Transform a str datetime to a pandas Timestamp (if not a Timestamp yet)"""
    if type(date_time)==str :
        return Timestamp(date_time)
    else :
        return date_time



def get_lat_lon_from_subdomain(domain:list[float])-> float:
    """Get latitude/longitude minimum and maximun from a list of
    float [lonmin, lonmax, latmin, latmax]"""
    lon_min = domain[0] ; lon_max = domain[1]
    lat_min = domain[2] ; lat_max = domain[3]
    
    return lon_min, lon_max, lat_min, lat_max



def select_Tmatrix_column(momentsDict:dict[int],hydrometeor:str,mask:ndarray,concentration:ndarray,Fw:ndarray,tmatrix_param:dict):
    """ Compute 3d parameter in Tmatrix table 
    TODO : change the reading of Tmatrix table so that the parameters are kept the same !

    * pristine ice (treated as 2-moment => the concentration can vary in MesoNH and is registered)
    * rainwater : 2-moment for all LIMA versions / 1-moment for ICE3
    * all other species (graupel / snow / wet graupel / wet snow ? )=> treated as 1-moment but can have a variable wet fraction
    => 3d parameter = Fw
    """
    
    if momentsDict[hydrometeor] == 2 :
        col_name = 'Nc'
        field_temp = concentration[mask]
        col = 'expCC'
    else:
        col_name = "Fw"
        field_temp = Fw
        col = 'Fw'
    col_min = tmatrix_param[f'{col}min'][hydrometeor]
    col_step = tmatrix_param[f'{col}step'][hydrometeor]
    col_max = tmatrix_param[f'{col}max'][hydrometeor]
    return field_temp, col_min, col_max, col_step, col_name