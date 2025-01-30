#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd



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



def hydrometeorModel_from_hydrometeorDict(hydrometeors:dict,quiet=False) -> list:
    """
    Make the correspondance between available hydrometeor keys in the model and
    the desired keys given in the configuration file.
    
    Available keys/hydrometeors in the model :
    - 'vv' : water vapour
    - 'cc' : cloud water
    - 'rr' : rain
    - 'ss' : snow
    - 'gg' : graupel
    - 'ii' : pristine ice
    - 'hh' : hail
    """
    model_hydrometeors=['vv','cc','rr','ss','gg','ii','hh']
    hydrometeors_to_extract = []
    for key in hydrometeors.keys() :
        if key in model_hydrometeors:
            hydrometeors_to_extract += [key]
    if not quiet :
        print('\tTo extract in model hydrometeor fields :',hydrometeors_to_extract)
    return(hydrometeors_to_extract)



def hydrometeorTmatrix_from_hydrometeorDict(hydrometeors:dict,quiet=False) -> list:
    """
    Make the correspondance between available hydrometeor keys in the Tmatrix 
    tables and the desired keys given in the configuration file.
    
    Available keys/hydrometeors in the model :
    - 'rr' : rain
    - 'ss' : snow
    - 'ws' : wet snow
    - 'ii' : pristine ice
    - 'gg' : graupel
    - 'wg' : wet graupel
    - 'hh' : hail
    - 'wh' : wet hail
    """
    Tmatrix_hydrometeors = ['rr','ss','ws','ii','gg','wg','hh','wh']
    hydrometeors_to_extract = []
    for key in hydrometeors.keys() :
        if key in Tmatrix_hydrometeors:
            hydrometeors_to_extract += [key]
    if not quiet :
        print('\tTo extract in Tmatrix tables :',hydrometeors_to_extract)
    return(hydrometeors_to_extract)



def format_date_time_argument(date_time:str|pd.Timestamp):
    if type(date_time)==str :
        return pd.Timestamp(date_time)
    else :
        return date_time