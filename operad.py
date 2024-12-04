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
    * netcdf file with modele points i,j,k and dpol var Zh, Zdr, Kdp, Rhohv 
    
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
import operad_lib as ope_lib
import read_tmatrix as read_tmat
import save_dpolvar as save
import csv_lib as csv_lib
from compute_var_pol import compute_individual_hydrometeor_with_Tmatrix


# ===== Configuration files ===== #
#import common_settings as settings
import operad_conf as cf


# ===== Begining of the program ===== #
begining_program_timer = tm.time()

model = 'Arome'#sys.argv[1]
dateconf = '20220620'#sys.argv[2]
micro = 'ICE3'#sys.argv[3]


# ===== Reading dates in case(s) study csv file
studyCases = csv_lib.read_csv_file(csv_file_path=cf.csvPath,
                                   which_date=dateconf,
                                   time_columns_name=cf.csv_time_columns,
                                   csv_datetime_format=cf.csv_datetime_format,
                                   csv_delimiter=cf.csv_delimiter,
                                   )


# ===== Loop over study cases ===== #
previous_radar_band = None

for _,csv_row in studyCases.iterrows():
    deb,fin,run,radar_band,output_path = csv_lib.extract_csv_info(csv_row=csv_row,
                                                                  time_columns_name=cf.csv_time_columns,
                                                                  model_run_column_name=cf.csv_model_run_column,
                                                                  microphysics_scheme=micro,
                                                                  )
    lat_min,lat_max,lon_min,lon_max = csv_lib.extract_csv_domain(csv_row=csv_row,
                                                                 domain_columns_name=cf.csv_domain_columns)
    
    # ----- Testing existence of the output directory (or creates it) ----- #
    ope_lib.create_tree_structure_outFiles(output_path)
    if cf.singletype :
        outpath_singleType = f"{cf.outPath}/{deb.strftime('%Y%m%d')}/{str(run).zfill(2)}Z_{micro}_single_type_hydrometeor"
        ope_lib.create_tree_structure_outFiles(outpath_singleType)
        
    # ----- Create the datetime list for the following loop ----- #
    datetimelist = ope_lib.create_datetime_list(deb,fin,cf.step)
    
    
    liste_var_pol = ["Zhh", "Zdr", "Kdp","Rhohv"] # A SUPPRIMER ?
    print([[hydro,method] for hydro,method in zip(cf.list_types_tot,cf.method)])

    # ----- Reading Tmatrix tables ----- #
    if radar_band != previous_radar_band :
        print(f"\nReading Tmatrix tables for {micro} in {radar_band} band")
        deb_timer = tm.time()
        dict_Tmatrix = read_tmat.Read_TmatrixClotilde(cf.pathTmat,
                                                      radar_band,
                                                      micro,
                                                      )
        LAM = dict_Tmatrix['LAMmin']["rr"]/1000.
        previous_radar_band = radar_band
        print("End reading Tmatrix tables in",round(tm.time()- deb_timer,2),"seconds")
    else :
        print(f"Tmatrix tables for {micro} in {radar_band} band already in memory")
        
    # ----- Loop over timesteps ----- #
    extract_once = True
    for datetime in datetimelist: 
        outFile = ope_lib.define_out_filepath(datetime, output_path,
                                              micro, model, radar_band,
                                              )
        model_file_path= ope_lib.define_model_path(model, datetime,
                                                   run, csv_row,
                                                   micro, deb,
                                                   )

        # ----- Testing existence of the output file ----- #
        if Path(outFile).exists() and (cf.singletype == False) :
            print("netcdf file for",datetime.strftime('%H:%M'),"already exists")
            continue   
        
        # ----- If not, creation of the file ----- #
        print("\n-------",datetime,"-------")
        # ----- Reading model variables ----- #
        print("\nReading model variables") ; deb_timer = tm.time()  
        if extract_once :
            [M, Tc, CC, CCI, lat,lon, X, Y, Z] = ope_lib.read_model_variables(model, datetime,
                                                                              run, micro,
                                                                              lon_min, lon_max,
                                                                              lat_min, lat_max,
                                                                              model_file_path,
                                                                              extract_once,
                                                                              )
            extract_once = False
        else :
            [M, Tc, CC, CCI, _,_, X, Y, Z] = ope_lib.read_model_variables(model, datetime,
                                                                          run, micro,
                                                                          lon_min, lon_max,
                                                                          lat_min, lat_max,
                                                                          model_file_path,
                                                                          extract_once,
                                                                          )
        print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
        
        # ----- Compute radar geometry variables ----- #
        print("\nCompute radar geometry and mixed phase") ; deb_timer = tm.time()
        np.seterr(invalid='ignore') # silence warning of invalid division 0 by 0 (result in a nan)
        el = np.zeros(Tc.shape)
        mask_distmax = (el >= 0.)
        if (model=="MesoNH"):
            [mask_distmax, el] = ope_lib.compute_radargeo(X, Y, Z,dict_Tmatrix['ELEVmax']["rr"])
        mask_precip_dist = ope_lib.mask_precip(mask_distmax, M, dict_Tmatrix['expMmin'], micro)
        
        # ------ Mixed phase parametrization -------- #
        [M, Fw] = ope_lib.compute_mixedphase(M,Tc, cf.MixedPhase, dict_Tmatrix['expMmin'], micro) 
        print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    
        # Initialization of dict(dpol_var_dict) --> contains all 3D dpol variables (all hydrometeor included)
        dpol_var_dict = {var:np.zeros(Tc.shape) for var in cf.dpol_var_to_calc}
            
        # ----- Loop over hydrometeor types ----- #
        print("\nLoop over hydrometeor types :") ; deb_timer = tm.time()
        for hydro,moment,method in zip(cf.list_types_tot,cf.n_moments_model,cf.method):
            print(f"\t{hydro} is {moment}-moment with {micro}.\n\tChosen method : {method}")
            
            dict_save_singleType = {}
            if (cf.singletype):
                outFileType = outpath_singleType + f"/dpolvar_{model}_{micro}_{radar_band}_{datetime.strftime('%Y%m%d_%H%M')}_{hydro}_{method}.nc"
                for arg,val in zip(['X','Y','Z','lat','lon','echeance','path_save_singleType'],[X,Y,Z,lat,lon,datetime,outFileType]):
                    dict_save_singleType[arg] = val
                if Path(outFile).exists() and Path(outFileType).exists() :
                    print("netcdf file for",datetime.strftime('%H:%M'),"and individual",hydro,"already exist")
                    continue

            if method == 'Tmatrix':
                dpol_var_dict = compute_individual_hydrometeor_with_Tmatrix(
                                hydrometeor=hydro, nmoment=moment, Tmatrix=dict_Tmatrix,
                                contents = M, elevation= el, temperature=Tc, waterFraction=Fw,
                                N_rain=CC, N_ice=CCI, mask_precip_dist=mask_precip_dist,
                                dpolVar_dict=dpol_var_dict, radar_wavelenght=LAM,
                                args_for_saving_singleType=dict_save_singleType,
                                                                            )
                
            #elif method == 'Rayleigh' :
            #    dpol_var_dict = compute_individual_hydrometeor_with_Rayleigh()
            else :
                print('Not a valid methodology :',method)
                
            
        for var in cf.dpol_var_to_calc:
            dpol_var_dict[var][~mask_precip_dist] = np.nan
        
        # ----- dpol var calculation
        dpol_var_dict["Zhh"] = np.copy(dpol_var_dict["Zhhlin"])
        dpol_var_dict["Zhh"][dpol_var_dict["Zhhlin"]>0] = ope_lib.Z2dBZ(dpol_var_dict["Zhhlin"][dpol_var_dict["Zhhlin"]>0])
        dpol_var_dict["Zdr"] = np.copy(dpol_var_dict["Zhhlin"])
        dpol_var_dict["Zdr"][(dpol_var_dict["Zhhlin"]>0) & (dpol_var_dict["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                (dpol_var_dict["Zhhlin"]/dpol_var_dict["Zvvlin"])[(dpol_var_dict["Zhhlin"]>0) & (dpol_var_dict["Zvvlin"]>0)])
        dpol_var_dict["Rhohv"] = np.sqrt(np.divide(dpol_var_dict["S11S22"], dpol_var_dict["S11S11"]*dpol_var_dict["S22S22"]))
        print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")
        
        # ----- Save dpol var for all hydro in netcdf and/or npz file
        if not Path(outFile).exists():
            save.save_dpolvar(M, CC, CCI, dpol_var_dict, Tc, Z, X, Y,lat,lon,datetime,outFile)
    
        del dpol_var_dict

end_program_timer = tm.time()
elapsed_time = end_program_timer - begining_program_timer
print("End of the program in",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")
