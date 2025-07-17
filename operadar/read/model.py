#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time as tm
from pathlib import Path
from numpy import ndarray



def read_model_file(filePath:Path,
                    modelname:str,
                    micro_scheme:str,
                    real_case:bool,
                    domain:list[float]|None,
                    hydrometeorMoments:dict[int],
                    verbose:bool,
                    )-> tuple[ndarray,ndarray,ndarray,ndarray,ndarray,dict[ndarray],dict[ndarray],ndarray]:
    """Read model file (either Arome or MesoNH)"""
    
    print("Reading model variables")
    deb_timer = tm.time()
    
    if (modelname=="MesoNH"):
        from operadar.read.mesonh import read_mesonh
        [X, Y, Alt, lon, lat, M, Nc, Tc] = read_mesonh(filePath=filePath,
                                                         micro=micro_scheme,
                                                         subDomain=domain,
                                                         hydrometeorMoments=hydrometeorMoments,
                                                         real_case = real_case,
                                                         verbose = verbose)
    elif (modelname=="Arome"):
        from operadar.read.arome import read_arome
        [X, Y, Alt, lon, lat, M, Nc, Tc] = read_arome(filePath=filePath,
                                                      hydrometeorMoments=hydrometeorMoments,
                                                      subDomain=domain,
                                                      verbose=verbose,
                                                      )   
    
    else :
        print('_____________')
        print('/!\ ERROR /!\ :',modelname,'is not a valid name. Must be either "Arome" or "MesoNH".')
        sys.exit()
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    
    return [X, Y, Alt, lon, lat, M, Nc, Tc]
