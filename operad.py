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

    # ----- Reading Tmatrix tables ----- #
    if radar_band != previous_radar_band :
        print(f"\nReading Tmatrix tables for {micro} in {radar_band} band")
        deb_timer = tm.time()
        dict_Tmatrix = read_tmat.Read_TmatrixClotilde(cf.pathTmat,
                                                      radar_band,
                                                      micro,
                                                      cf.list_types_tot,
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
    
        # Initialization of dict(Vm_k) --> contains all 3D dpol variables (all hydrometeor included)
        Vm_k = {var:np.zeros(Tc.shape) for var in cf.liste_var_calc}
            
        # ----- Loop over hydromet types ----- #
        print("\nLoop over hydrometeor types :") ; deb_timer = tm.time()
        for hydromet,moment,method in zip(cf.list_types_tot,cf.n_moments_model,cf.method):
            print(f"\t{hydromet} is {moment}-moment with {micro}.\n\tChosen method : {method}")
            
            if (cf.singletype):
                outFileType = outpath_singleType + f"/dpolvar_{model}_{micro}_{radar_band}_{datetime.strftime('%Y%m%d_%H%M')}_{hydromet}_{method}.nc"         
                Vm_t = {var:np.zeros(Tc.shape) for var in cf.liste_var_calc}
            
                if Path(outFile).exists() and Path(outFileType).exists() :
                    print("netcdf file for",datetime.strftime('%H:%M'),"and individual",hydromet,"already exist")
                    continue
            
            #if method == 'Rayleigh':
            #else :    
            
            # Compute single type mask
            [mask_tot, M_masked, elevation_masked, Tc_masked, P3_masked] = ope_lib.singletype_mask(M[hydromet],
                                                                                                el, Tc, Fw,
                                                                                                CC, CCI,
                                                                                                mask_precip_dist,
                                                                                                dict_Tmatrix['expMmin'],
                                                                                                hydromet, moment,
                                                                                                )
            
            
            # Extract scattering coefficients for singletype
            [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11] = read_tmat.get_scatcoef(dict_Tmatrix, hydromet, moment,
                                                                                        elevation_masked,
                                                                                        Tc_masked,
                                                                                        P3_masked,
                                                                                        M_masked,
                                                                                        cf.n_interpol,
                                                                                        shutdown_warnings=True,
                                                                                        )
                
            # Single type dpol var computation       
            Vm_temp = {var:Vm_k[var][mask_tot] for var in cf.liste_var_calc}
            Vm_temp["Zhhlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S22carre #lin = linear
            Vm_temp["Zvvlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S11carre
            Vm_temp["Kdp"] = 180.*1e3/math.pi*LAM*ReS22fmS11f
            Vm_temp["S11S22"] = ReS22S11**2+ImS22S11**2
            Vm_temp["S11S11"] = np.copy(S11carre)
            Vm_temp["S22S22"] = np.copy(S22carre)
            
            # Addition of scattering coef for all hydromet
            for var in cf.liste_var_calc:
                Vm_k[var][mask_tot]+=Vm_temp[var]
                if (cf.singletype):
                    Vm_t[var][mask_tot]=Vm_temp[var]
                    Vm_t[var][~mask_tot] = np.nan 
            
            del S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
            del Vm_temp, mask_tot

            #  Dpol variables for single hydrometeor types 
            if (cf.singletype):
                Vm_t["Zhh"] = np.copy(Vm_t["Zhhlin"])
                Vm_t["Zhh"][Vm_t["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_t["Zhhlin"][Vm_t["Zhhlin"]>0])
                Vm_t["Zdr"] = np.copy(Vm_t["Zhhlin"])
                Vm_t["Zdr"][(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                        (Vm_t["Zhhlin"]/Vm_t["Zvvlin"])[(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)])
                Vm_t["Rhohv"] = np.sqrt(np.divide(Vm_t["S11S22"], Vm_t["S11S11"]*Vm_t["S22S22"]))
                
                # Writing dpol var for a single hydrometeor type hydromet
                save.save_dpolvar({hydromet:M[hydromet]}, CC, CCI,  Vm_t, Tc, Z, X, Y, lat,lon,datetime,outFileType,singleType=True)
                del Vm_t
                
        # ------------ end
        for var in cf.liste_var_calc:
            Vm_k[var][~mask_precip_dist] = np.nan 
        
        # ----- dpol var calculation
        Vm_k["Zhh"] = np.copy(Vm_k["Zhhlin"])
        Vm_k["Zhh"][Vm_k["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_k["Zhhlin"][Vm_k["Zhhlin"]>0])
        Vm_k["Zdr"] = np.copy(Vm_k["Zhhlin"])
        Vm_k["Zdr"][(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                (Vm_k["Zhhlin"]/Vm_k["Zvvlin"])[(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)])
        Vm_k["Rhohv"] = np.sqrt(np.divide(Vm_k["S11S22"], Vm_k["S11S11"]*Vm_k["S22S22"]))
        print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")
        
        # ----- Save dpol var for all hydromet in netcdf and/or npz file
        if not Path(outFile).exists():
            save.save_dpolvar(M, CC, CCI, Vm_k, Tc, Z, X, Y,lat,lon,datetime,outFile)
    
        del Vm_k

end_program_timer = tm.time()
elapsed_time = end_program_timer - begining_program_timer
print("End of the program in",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")
