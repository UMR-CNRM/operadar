#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018

@author: augrosc & lebastardt & montangonm & sinhorin & davidcl

Dual Polarization Radar Forward operator :
calculation of dpol radar variables for each model grid point
- INPUT : 
    * a model file (MesoNH netcdf or AROME fa): modelfile => contains Z (altitude), Tc (temperature), hydrometeor contents M and concentration CC
    * the scattering coefficients tables (output files of Tmatrix) recorded for
    a range of T, M, CC (if 2 moments) 
- OUTPUT : 
    * netcdf or .npz file with modele points i,j,k and dpol var Zh, Zdr, Kdp, Rhohv 
    
The modifications are recorded via git:
    git status (to see the modified and not recorded files)
    git log (for the historic of all modifications)
    git gui (visual tool for git)
"""

# External modules
import sys
import math
import time as tm
import numpy as np
import pandas as pd
import datetime as dt

from pathlib import Path
sys.path.insert(0, "./lib")


# 0perad modules
import operad_conf as cf
import operad_lib as ope_lib
import read_tmatrix as read_tmat
import save_dpolvar as save
from operadar.operadar_utils import (
    hydrometeorTmatrix_from_hydrometeorDict,
    format_date_time_argument,
)



def operad(filename:Path, date_time:str|pd.Timestamp,
           read_tmatrix:bool=True, extract_once:bool=True,
           outPath:str=cf.outPath):
    
    begin_program_timer = tm.time()
    
    # Imports
    if (cf.model=="MesoNH"): import read_mesonh as meso
    elif (cf.model=="Arome"): import read_arome as aro

    # Create or check tree structure of the output directory path
    save.create_tree_structure_outFiles(Path(outPath))
    
    # Format datetime argument for the output file name and check existence
    date_time = format_date_time_argument(date_time)
    outFilePath = Path(outPath + f"dpolvar_{cf.model}_{cf.micro_scheme}_{cf.radar_band}band_{date_time.strftime('%Y%m%d_%H%M')}")
    if not outFilePath.exists():
        # do something
    else :
        print("Exists :",outFilePath)


    elapsed_time = tm.time() - begin_program_timer
    print("Elapsed time :",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")


if __name__ == '__main__':
    filepath = Path(sys.argv[1])
    date_time = sys.argv[2]
    operad(filename=filepath,date_time=date_time)
