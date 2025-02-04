#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd



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



def format_date_time_argument(date_time:str|pd.Timestamp):
    if type(date_time)==str :
        return pd.Timestamp(date_time)
    else :
        return date_time



def get_lat_lon_from_subdomain(domain:list[float])-> float:
    """Get latitude/longitude minimum and maximun from a list of
    float [lonmin, lonmax, latmin, latmax]"""
    lon_min = domain[0] ; lon_max = domain[1]
    lat_min = domain[2] ; lat_max = domain[3]
    
    return lon_min, lon_max, lat_min, lat_max