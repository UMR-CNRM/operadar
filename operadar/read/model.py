#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time as tm
from pathlib import Path
import operad_conf as cf



def read_model_file(filePath:Path, domain:list[float]|None, extract_once:bool=True):
    
    print("Reading model variables")
    deb_timer = tm.time()
    
    if (cf.model=="MesoNH"):
        from read.mesonh import read_mesonh
        [M, Tc, CC, CCI, lat,lon, X, Y, Z] = read_mesonh(filePath=filePath, micro=cf.micro_scheme,
                                                         subDomain=domain,
                                                         hydrometeors = cf.moments,
                                                         real_case = cf.real_case)
    elif (cf.model=="Arome"):
        from read.arome import read_arome
        [X, Y, Z, lon, lat, M, Nc, Tc] = read_arome(filePath=filePath, micro=cf.micro_scheme,
                                                    extract_once=extract_once, hydrometeors=cf.moments,
                                                    subDomain=domain,
                                                    )   
    
    else :
        print('_____________')
        print('/!\ ERROR /!\ :',cf.model,'is not a valid name. Must be either "Arome" or "MesoNH".')
        sys.exit()
    
    print("--> Done in",round(tm.time()- deb_timer,2),"seconds")
    
    return [X, Y, Z, lon, lat, M, Nc, Tc]