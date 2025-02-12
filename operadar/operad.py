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


# 0perad modules
import operadar.operad_conf as cf
from operadar.read.model import read_model_file
from operadar.utils.masking import mask_precipitations
from operadar.radar.geometry import compute_radar_geometry
from operadar.read.tmatrix_tables import read_Tmatrix_Clotilde
from operadar.utils.formats_data import format_date_time_argument
from operadar.microphysics.mixed_phase import compute_mixed_phase
from operadar.radar.dualpol_variables import compute_dualpol_variables
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.save.save_dpolvar import create_tree_structure_outFiles, save_dpolvar




def operad(filename:Path, date_time:str|pd.Timestamp,
           read_tmatrix:bool=True, extract_once:bool=True,
           outPath:str=cf.outPath, radar_band:str=cf.radar_band,
           Tmatrix_params:dict={}, exp_name:str=cf.experience_name,
           subDomain:list[float]|None=cf.subDomain,
           ):
    
    begin_program_timer = tm.time()

    # Create or check tree structure of the output directory path
    create_tree_structure_outFiles(Path(outPath))
    
    # Format datetime argument for the output file name and check existence
    date_time = format_date_time_argument(date_time)
    outFilePath = outPath + f"dpolvar_{cf.model}_{cf.micro_scheme}_{radar_band}band_{date_time.strftime('%Y%m%d_%H%M')}"
    
    if not Path(outFilePath).with_suffix('.nc').exists():
        
        # Read Tmatrix tables (files from Clotilde)
        if read_tmatrix :
            Tmatrix_hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.moments,datatype='tmatrix')
            Tmatrix_params = read_Tmatrix_Clotilde(band=radar_band,hydrometeors=Tmatrix_hydromet_list)
        
        # Read model variables
        input_file_path = Path(cf.input_file_dir) / Path(exp_name) / filename
        if extract_once:
            [X, Y, Z, lon, lat, M, Nc, Tc] = read_model_file(filePath=input_file_path,
                                                             date_time=date_time,
                                                             domain=subDomain,
                                                             )
        else :
            [X, Y, Z, _, _, M, Nc, Tc] = read_model_file(filePath=input_file_path,
                                                         date_time=date_time,
                                                         domain=subDomain,
                                                         extract_once=False,
                                                         )
        # Compute radar geometry
        mask_dist_max, elevations = compute_radar_geometry(X=X, Y=Y, Z=Z, Tc=Tc,
                                                           elev_max=Tmatrix_params['ELEVmax']["rr"])
        
        # Mask precipitations
        mask_precip = mask_precipitations(contents=M,
                                          expMmin=Tmatrix_params['expMmin']["rr"],
                                          hydrometeorMoments=cf.moments)
        
        # Combine masks
        partial_mask = (mask_precip & mask_dist_max)
        
        # Compute mixed phase parametrization
        [M, Nc,Fw] = compute_mixed_phase(contents=M,
                                         concentrations=Nc,
                                         expMmin=Tmatrix_params['expMmin']["rr"]) 
        
        # Compute dual-pol radar variables
        dpolDict = compute_dualpol_variables(temperature=Tc,
                                             mask_precip_dist=partial_mask,
                                             elev=elevations, Fw=Fw,
                                             contents=M,
                                             concentrations=Nc,
                                             tmatrix_param=Tmatrix_params,
                                             save_single=outFilePath,
                                             X=X, Y=Y, Z=Z, lat=lat,
                                             lon=lon, date_time=date_time)
        
        # Saving file
        save_dpolvar(M, Nc, dpolDict, Tc, Z, X, Y,lat,lon,date_time,Path(outFilePath))
    else :
        print("File exists at :",outFilePath+'.nc')

    elapsed_time = tm.time() - begin_program_timer
    print("Elapsed time :",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")




if __name__ == '__main__':
    filepath = Path(sys.argv[1])
    date_time = sys.argv[2]
    operad(filename=filepath,date_time=date_time)
