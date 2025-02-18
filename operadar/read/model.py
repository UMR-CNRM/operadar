#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time as tm
from pathlib import Path
from numpy import ndarray
from pandas import Timestamp
import operadar.operadar_conf as cf



def read_model_file(filePath:Path,
                    domain:list[float]|None,
                    verbose:bool,
                    )-> tuple[ndarray,ndarray,ndarray,ndarray,ndarray,dict[ndarray],dict[ndarray],ndarray]:
    """Read model file (either Arome or MesoNH)"""
    
    print("Reading model variables")
    deb_timer = tm.time()
    
    if (cf.model=="MesoNH"):
        from operadar.read.mesonh import read_mesonh
        [M, Tc, CC, CCI, lat,lon, X, Y, Z] = read_mesonh(filePath=filePath,
                                                         micro=cf.micro_scheme,
                                                         subDomain=domain,
                                                         hydrometeorMoments=cf.hydrometeors_moments,
                                                         real_case = cf.real_case)
    elif (cf.model=="Arome"):
        from operadar.read.arome import read_arome
        [X, Y, Z, lon, lat, M, Nc, Tc] = read_arome(filePath=filePath,
                                                    hydrometeorMoments=cf.hydrometeors_moments,
                                                    subDomain=domain,
                                                    verbose=verbose,
                                                    )   
    
    else :
        print('_____________')
        print('/!\ ERROR /!\ :',cf.model,'is not a valid name. Must be either "Arome" or "MesoNH".')
        sys.exit()
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    
    return [X, Y, Z, lon, lat, M, Nc, Tc]