#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018

@author: davidcl & augrosc

Computation of dpol radar variables for each hydrometeor type and with different method (Rayleigh or T-matrix)
- INPUT : 
    * a model file (MesoNH netcdf or AROME fa): modelfile => contains Z (altitude), temperature (temperature), hydrometeor contents contents and concentration N_rain
    * the scattering coefficients tables (output files of Tmatrix) recorded for
    a range of T, contents, N_rain (if 2 moments) 
- OUTPUT : 
    * netcdf file with modele points i,j,k and dpol var Zh, Zdr, Kdp, Rhohv
"""
# External modules
import sys
import math
import numpy as np
import pandas as pd
from pathlib import Path

# 0perad modules
sys.path.insert(0, "./lib")
import operad_lib as ope_lib
import read_tmatrix as read_tmat
import save_dpolvar as save
import csv_lib as csv_lib
import operad_conf as cf

def compute_dict_for_individual_hydrometeor(method,
                                            hydrometeor:str,
                                            nmoment:int,
                                            Tmatrix:dict,
                                            contents:dict,
                                            elevation:np.array,
                                            temperature,waterFraction,
                                            N_rain,N_ice,
                                            mask_precip_dist,
                                            dpolVar_dict:dict,
                                            radar_wavelenght,
                                            Z, X, Y, lat,lon,echeance,path_save_singleType,
                                            ):
    if method == 'Tmatrix':
        compute_individual_hydrometeor_with_Tmatrix(hydrometeor, nmoment, Tmatrix,
                                                    contents, elevation, temperature, waterFraction, N_rain,N_ice,
                                                    mask_precip_dist, dpolVar_dict,
                                                    radar_wavelenght,
                                                    )


def compute_individual_hydrometeor_with_Tmatrix(
    hydrometeor:str, nmoment:int, Tmatrix:dict, contents:dict, elevation:np.ndarray,
    temperature:np.ndarray, waterFraction:np.ndarray, N_rain:np.ndarray,
    N_ice:np.ndarray, mask_precip_dist:np.ndarray, dpolVar_dict:dict,
    radar_wavelenght:float, **args_for_saving_singleType:dict,
    )-> dict :
    
    # Compute single type mask
    [mask_tot, M_masked, elevation_masked, Tc_masked, P3_masked] = ope_lib.singletype_mask(contents[hydrometeor],
                                                                                           elevation, temperature, waterFraction,
                                                                                           N_rain, N_ice,
                                                                                           mask_precip_dist,
                                                                                           Tmatrix['expMmin'],
                                                                                           hydrometeor, nmoment,
                                                                                           )
    
    # Extract scattering coefficients for singletype
    [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11] = read_tmat.get_scatcoef(Tmatrix, hydrometeor, nmoment,
                                                                                   elevation_masked,
                                                                                   Tc_masked,
                                                                                   P3_masked,
                                                                                   M_masked,
                                                                                   cf.n_interpol,
                                                                                   shutdown_warnings=True,
                                                                                   )
        
    # Single type dpol var computation      
    temp_dict = {var:dpolVar_dict[var][mask_tot] for var in cf.dpol_var_to_calc}
    temp_dict["Zhhlin"]= 1e18*radar_wavelenght**4./(math.pi**5.*0.93)*4.*math.pi*S22carre #lin = linear
    temp_dict["Zvvlin"]= 1e18*radar_wavelenght**4./(math.pi**5.*0.93)*4.*math.pi*S11carre
    temp_dict["Kdp"] = 180.*1e3/math.pi*radar_wavelenght*ReS22fmS11f
    temp_dict["S11S22"] = ReS22S11**2+ImS22S11**2
    temp_dict["S11S11"] = np.copy(S11carre)
    temp_dict["S22S22"] = np.copy(S22carre)
    
    # Addition of scattering coef for all hydrometeor
    if (cf.singletype):
        dpolVar_dict_1hydro = {var:np.zeros(temperature.shape) for var in cf.dpol_var_to_calc}
    for var in cf.dpol_var_to_calc:
        dpolVar_dict[var][mask_tot]+=temp_dict[var]
        if (cf.singletype):
            dpolVar_dict_1hydro[var][mask_tot]=temp_dict[var]
            dpolVar_dict_1hydro[var][~mask_tot] = np.nan 
    
    del S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
    del elevation_masked, Tc_masked, M_masked, temp_dict, mask_tot

    #  Dpol variables for single hydrometeor types
    if (cf.singletype) and not Path(args_for_saving_singleType['path_save_singleType']).exists():
        dpolVar_dict_1hydro = calculate_var_pol(dpolVar_dict_1hydro)
        # Writing dpol var for a single hydrometeor type hydrometeor
        save.save_dpolvar(M={hydrometeor:contents[hydrometeor]}, CC=N_rain, CCI=N_ice, Vm_k=dpolVar_dict_1hydro,Tc=temperature,
                          Z=args_for_saving_singleType['Z'], X=args_for_saving_singleType['X'],
                          Y=args_for_saving_singleType['Y'], lat=args_for_saving_singleType['lat'],
                          lon=args_for_saving_singleType['lon'], datetime=args_for_saving_singleType['echeance'],
                          outfile=args_for_saving_singleType['path_save_singleType'],singleType=True,
                          )
        del dpolVar_dict_1hydro
        
    return dpolVar_dict

    
#def compute_individual_hydrometeor_with_Rayleigh() : return ?


def calculate_var_pol(dpolVar_dict):
    for var in cf.liste_var_pol :
        if var=="Zhh" or var=="Zh" or var=="zh" :
            dpolVar_dict = calculate_Zh(dpolVar_dict)
        elif var=="Zdr" or var=="zdr" :
            dpolVar_dict = calculate_Zdr(dpolVar_dict)
        elif var=="Kdp" or var=="kdp" :
            pass
            #dpolVar_dict = calculate_Kdp(dpolVar_dict)
        elif var=="Rhohv" or var=="rhohv" or var=="rho" :
            dpolVar_dict = calculate_RhoHV(dpolVar_dict)
            
    return dpolVar_dict


def calculate_Zh(dpolVar_dict):
    mask_zhh_linear_pos = dpolVar_dict["Zhhlin"]>0
    dpolVar_dict["Zhh"] = np.copy(dpolVar_dict["Zhhlin"])
    dpolVar_dict["Zhh"][mask_zhh_linear_pos] = ope_lib.Z2dBZ(dpolVar_dict["Zhhlin"][mask_zhh_linear_pos])
    return dpolVar_dict
        
        
def calculate_Zdr(dpolVar_dict):
    mask_zhh_zvv_linear_pos = (dpolVar_dict["Zhhlin"]>0) & (dpolVar_dict["Zvvlin"]>0)
    dpolVar_dict["Zdr"] = np.copy(dpolVar_dict["Zhhlin"])
    dpolVar_dict["Zdr"][mask_zhh_zvv_linear_pos] = ope_lib.Z2dBZ((dpolVar_dict["Zhhlin"]/dpolVar_dict["Zvvlin"])[mask_zhh_zvv_linear_pos])
    return dpolVar_dict
    
"""
def calculate_Kdp(method):
"""
def calculate_RhoHV(dpolVar_dict):
    dpolVar_dict["Rhohv"] = np.sqrt(np.divide(dpolVar_dict["S11S22"], dpolVar_dict["S11S11"]*dpolVar_dict["S22S22"]))
    return dpolVar_dict