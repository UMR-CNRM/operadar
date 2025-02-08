#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import time as tm
import numpy as np
from pathlib import Path
from pandas import Timestamp

import operadar.operad_conf as cf
from operadar.read.tmatrix_tables import get_scatcoef
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.utils.masking import mask_hydrometeor
from operadar.utils.formats_data import select_Tmatrix_column
from operadar.save.save_dpolvar import save_dpolvar



def compute_dualpol_variables(temperature:np.ndarray, mask_precip_dist:np.ndarray,
                              elev:np.ndarray, Fw:np.ndarray, contents:dict[np.ndarray],
                              concentrations:dict[np.ndarray], tmatrix_param:dict,
                              save_single:str, X:np.ndarray, Y:np.ndarray, Z:np.ndarray,
                              lon:np.ndarray, lat:np.ndarray, date_time:Timestamp)-> dict[np.ndarray] :
    
    dpol_var = ["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.moments,datatype='tmatrix')
    
    print("Computation of",dpol_var,"for :") ; deb_timer = tm.time()
    dpol_var_dict = {var:np.zeros(temperature.shape) for var in dpol_var}
    
    for h in hydrometeors :
        print('\t- hydrometeor :',h)
        
        # Mask single hydrometeor type
        mask_content = mask_hydrometeor(contents[h],expMmin=tmatrix_param['expMmin'][h])
        mask_tot= mask_precip_dist & mask_content 

        dpol_h = compute_dualpol_single_hydrometeor(hydrometeor=h,
                                                    dpol_var=dpol_var,
                                                    mask_tot=mask_tot,
                                                    Tc=temperature,
                                                    el=elev, Fw=Fw,
                                                    content_h=contents[h],
                                                    concentration_h=concentrations[h],
                                                    tmatrix_param=tmatrix_param,
                                                    )
        
        # Addition of scattering coef for all hydromet
        for var in dpol_var:
            dpol_var_dict[var][mask_tot]+=dpol_h[var]
            if cf.singletype :
                dpol_h[var][mask_tot]=dpol_h[var]
                dpol_h[var][~mask_tot] = np.NaN
        
        # If saving single type, compute final polarimetric values
        if cf.singletype :
            dpol_h = compute_dpol_var(dpolDict=dpol_h)
            outFilePath = Path(f'{save_single}_{h}')
            save_dpolvar(contents[h], concentrations[h], dpol_h,
                         temperature, Z, X, Y, lat,lon,date_time,outFilePath)
        del dpol_h   
        
    for var in dpol_var:
        dpol_var_dict[var][~mask_precip_dist] = np.NaN 
        
    # Dpol var calculation over the sum of scatering coefficients and linear Z
    dpol_var_dict = compute_dpol_var(dpolDict=dpol_var_dict)
    
    print("  --> Done in",round(tm.time() - deb_timer,2),"seconds")    
    return dpol_var_dict  



def compute_dualpol_single_hydrometeor(hydrometeor:str, dpol_var:list, Tc:np.ndarray,
                                       el:np.ndarray, Fw:np.ndarray, content_h:np.ndarray,
                                       mask_tot:np.ndarray, concentration_h:np.ndarray,
                                       tmatrix_param:dict):
    
    elev_temp=el[mask_tot] 
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    content_temp=content_h[mask_tot]
    concentration_temp=concentration_h[mask_tot]
    
    # Define P3 : Nc (2 moments) or Fw (1 moment) hydrometeor_moment = cf.moments[hydrometeor]
    field_temp, colMin, colMax, colStep, colName = select_Tmatrix_column(momentsDict=cf.moments,
                                                                         hydrometeor=hydrometeor,
                                                                         concentration=concentration_temp,
                                                                         Fw=Fw_temp,
                                                                         tmatrix_param=tmatrix_param)
    
    # Extract scattering coefficients for singletype
    [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11] = get_scatcoef(tmatDict=tmatrix_param,
                                                                         hydrometeor=hydrometeor,
                                                                         colName_tmatrix=colName,
                                                                         colMin=colMin,
                                                                         colStep=colStep,
                                                                         colMax=colMax,
                                                                         el_temp=elev_temp,
                                                                         Tc_temp=Tc_temp,
                                                                         colTmat=field_temp,
                                                                         M_temp=content_temp,
                                                                         n_interpol=cf.n_interpol)
      
    # Compute dualpol variables from scattering coefficients
    dpolDict_h = {var:np.zeros(Tc.shape)[mask_tot] for var in dpol_var}
    dpolDict_h = dpol_var_from_scatcoefs(temp_dict=dpolDict_h,
                                         wavelength=tmatrix_param['LAMmin']['rr']/1000.,
                                         S11carre=S11carre, S22carre=S22carre,
                                         ReS22fmS11f=ReS22fmS11f, ReS22S11=ReS22S11,
                                         ImS22S11=ImS22S11)
    
    del S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
    del elev_temp, Tc_temp, content_temp, Fw_temp, concentration_temp
    
    return dpolDict_h



def dpol_var_from_scatcoefs(temp_dict,wavelength,S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11):
    temp_dict["Zhhlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*S22carre
    temp_dict["Zvvlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*S11carre
    temp_dict["Kdp"] = 180.*1e3/math.pi*wavelength*ReS22fmS11f
    temp_dict["S11S22"] = ReS22S11**2+ImS22S11**2
    temp_dict["S11S11"] = S11carre
    temp_dict["S22S22"] = S22carre
    
    return temp_dict



def compute_dpol_var(dpolDict:dict[np.ndarray]):
    # ZH
    dpolDict["Zhh"] = dpolDict["Zhhlin"]
    dpolDict["Zhh"][dpolDict["Zhhlin"]>0] = linear_to_dBZ(dpolDict["Zhhlin"][dpolDict["Zhhlin"]>0])
    # ZDR
    dpolDict["Zdr"] = dpolDict["Zhhlin"]
    mask_Zdr = (dpolDict["Zhhlin"]>0) & (dpolDict["Zvvlin"]>0)
    dpolDict["Zdr"][mask_Zdr] = linear_to_dBZ((dpolDict["Zhhlin"]/dpolDict["Zvvlin"])[mask_Zdr])
    # RHOHV
    dpolDict["Rhohv"] = np.sqrt(np.divide(dpolDict["S11S22"], dpolDict["S11S11"]*dpolDict["S22S22"]))
    
    return dpolDict



def linear_to_dBZ(Z:np.ndarray):
    """Linear to dBZ conversion.""" 
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Z[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp
 