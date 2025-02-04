#!/usr/bin/env python3
# -*- coding: utf-8 -*-



def link_varname_with_arome_name ():
    """Make the correspondance between the hydrometeor keys and the variable names commonly used in AROME."""
    model_hydrometeors=['vv','cc','rr','ii','ss','gg','hh']
    model_name=['HUMI.SPECIFI','CLOUD_WATER','RAIN','ICE_CRYSTAL','SNOW','GRAUPEL','HAIL']
    name_hydro_linked={t:model_name[it] for it,t in enumerate(model_hydrometeors)}
    return name_hydro_linked



def link_keys_with_available_hydrometeors(hydrometeors:dict,datatype:str,quiet=False) -> list:
    """Make the correspondance between available hydrometeor keys in the model or Tmatrix tables
    and the desired keys given in the configuration file.
    
    Available keys/hydrometeors in the model : 'vv','cc','rr','ss','gg','ii','hh'
    Available keys/hydrometeors in the Tmatrix tables : 'rr','ss','ws','ii','gg','wg','hh','wh'
    
    Legend:
    - 'vv' : water vapour
    - 'cc' : cloud water
    - 'rr' : rain
    - 'ss' : snow
    - 'ws' : wet snow
    - 'gg' : graupel
    - 'wg' : wet graupel
    - 'ii' : pristine ice
    - 'hh' : hail
    - 'wh' : wet hail
    """
    model_hydrometeors=['vv','cc','rr','ss','gg','ii','hh']
    Tmatrix_hydrometeors = ['rr','ss','ws','ii','gg','wg','hh','wh']
    hydrometeors_to_extract = []
    for key in hydrometeors.keys() :
        if key in model_hydrometeors:
            hydrometeors_to_extract += [key]
    if not quiet :
        print('\tTo extract in model hydrometeor fields :',hydrometeors_to_extract)
    return(hydrometeors_to_extract)

