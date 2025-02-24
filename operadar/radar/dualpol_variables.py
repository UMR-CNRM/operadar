#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import time as tm
import numpy as np
from pathlib import Path
from pandas import Timestamp

from operadar.operadar_conf import save_netcdf_single_hydrometeor, n_interpol,dpol2add
from operadar.save.save_dpolvar import save_netcdf
from operadar.utils.masking import mask_hydrometeor
from operadar.read.tmatrix_tables import get_scatcoef
from operadar.utils.formats_data import select_Tmatrix_column
from operadar.utils.make_links import link_keys_with_available_hydrometeors



def compute_dualpol_variables(temperature:np.ndarray,
                              mask_precip_dist:np.ndarray,
                              elev:np.ndarray, Fw:np.ndarray,
                              contents:dict[np.ndarray],
                              concentrations:dict[np.ndarray],
                              tmatrix_param:dict,
                              hydrometeorMoments:dict[int],
                              X:np.ndarray,
                              Y:np.ndarray,
                              Z:np.ndarray,
                              lon:np.ndarray,
                              lat:np.ndarray,
                              date_time:Timestamp,
                              output_file_path:Path,
                              append_in_fa:bool,
                              )-> dict[np.ndarray] :
    """Compute synthetic radar dual-polarimetrization variables for a given wavelength,
    microphysics, and mixed phase parametrization.

    Args:
        temperature (np.ndarray): 3D array.
        mask_precip_dist (np.ndarray): 3D array of bool.
        elev (np.ndarray): 3D array.
        Fw (np.ndarray): 3D array.
        contents (dict[np.ndarray]): dict of 3D array (one per hydrometeor).
        concentrations (dict[np.ndarray]): dict of 3D array (one per hydrometeor).
        tmatrix_param (dict): dict of Tmatrix parameters.
        hydrometeorMoments (dict of form {str : int})): dict containing the number of moments for each hydrometeor of the microphysics scheme.
        X (np.ndarray): 1D array of horizontal grid coordinates.
        Y (np.ndarray): 1D array of horizontal grid coordinates.
        Z (np.ndarray): vertical coordinates.
        lon (np.ndarray): 2D array.
        lat (np.ndarray): 2D array.
        date_time (Timestamp): temporal variable.
        output_file_path (Path): out file path.
        append_in_fa (bool): if True, delete contents and concentrations field when used in the calculation of the dpol variables to save space.

    Returns:
        dict[np.ndarray]: dictionary containing all the dual-polarimetrization fields.
    """

    dpol_var, scat_coefs = build_list_dpol_var_and_scatcoefs(dpol2add=dpol2add)
    
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                         datatype='tmatrix'
                                                         )
    print("Computation of",dpol_var,"and",scat_coefs,"for :") ; deb_timer = tm.time()
    dpol_var_dict = {var:np.zeros(temperature.shape) for var in dpol_var}
    
    for h in hydrometeors :
        print('\t- hydrometeor :',h)
        
        # Mask single hydrometeor type
        mask_content = mask_hydrometeor(content=contents[h],
                                        expMmin=tmatrix_param['expMmin'][h]
                                        )
        mask_tot = (mask_precip_dist & mask_content) 

        dpolDict = compute_scatcoeffs_single_hydrometeor(hydrometeor=h,
                                                         dpol_var=dpol_var,
                                                         scat_coefs=scat_coefs,
                                                         mask_tot=mask_tot,
                                                         Tc=temperature,
                                                         el=elev, Fw=Fw,
                                                         content_h=contents[h],
                                                         concentration_h=concentrations[h],
                                                         tmatrix_param=tmatrix_param,
                                                         hydrometeorMoments=hydrometeorMoments,
                                                    )

        # Addition of scattering coef for all hydromet
        for var in dpolDict.keys():
            dpol_var_dict[var][mask_tot]+=dpolDict[var]
            if save_netcdf_single_hydrometeor :
                dpol_h = {var:np.zeros(temperature.shape) for var in dpol_var}
                dpol_h[var][mask_tot]=dpolDict[var]
                dpol_h[var][~mask_tot]= np.nan
        del dpolDict
        
        # If saving single type, compute final polarimetric values
        if save_netcdf_single_hydrometeor :
            dpol_h = compute_dpol_var(dpolDict=dpol_h)
            outFilePath = Path(f'{output_file_path}_{h}')
            save_netcdf(X=X, Y=Y, Z=Z, lat=lat, lon=lon, 
                        datetime=date_time, dpolDict=dpol_h,
                        contentsDict={h:contents[h]},
                        concentrationsDict={h:concentrations[h]},
                        temperature=temperature,
                        outfile=outFilePath,
                        )
            del dpol_h   
        if append_in_fa : del concentrations[h],contents[h]
        
    for var in dpol_var:
        dpol_var_dict[var][~mask_precip_dist] = np.nan 
        
    # Dpol var calculation over the sum of scatering coefficients and linear Z
    dpol_var_dict = compute_dpol_var(dpolDict=dpol_var_dict)
    
    print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")    
    return dpol_var_dict  



def compute_scatcoeffs_single_hydrometeor(hydrometeor:str,
                                          dpol_var:list,
                                          scat_coefs:list,
                                          Tc:np.ndarray,
                                          el:np.ndarray,
                                          Fw:np.ndarray,
                                          content_h:np.ndarray,
                                          mask_tot:np.ndarray,
                                          concentration_h:np.ndarray,
                                          tmatrix_param:dict,
                                          hydrometeorMoments:dict[int],
                                        ) -> dict[np.ndarray]:
    """Compute radar scattering coefficients for a single hydrometeor class."""
    
    elev_temp=el[mask_tot]
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    content_temp=content_h[mask_tot]
    concentration_temp=concentration_h[mask_tot]
    
    # Define P3 : Nc (2 moments) or Fw (1 moment) 
    field_temp, colMin, colMax, colStep, colName = select_Tmatrix_column(momentsDict=hydrometeorMoments,
                                                                         hydrometeor=hydrometeor,
                                                                         concentration=concentration_temp,
                                                                         Fw=Fw_temp,
                                                                         tmatrix_param=tmatrix_param,
                                                                         )
    # Extract scattering coefficients for save_netcdf_single_hydrometeor
    scat_coefs_dict = get_scatcoef(tmatDict=tmatrix_param,
                                   scat_coefs=scat_coefs,
                                   hydrometeor=hydrometeor,
                                   colName_tmatrix=colName,
                                   colMin=colMin,
                                   colStep=colStep,
                                   colMax=colMax,
                                   el_temp=elev_temp,
                                   Tc_temp=Tc_temp,
                                   colTmat=field_temp,
                                   M_temp=content_temp,
                                   n_interpol=n_interpol,
                                   )
    # Compute dualpol variables from scattering coefficients
    dpolDict_h = dpol_var_from_scatcoefs(wavelength=tmatrix_param['LAMmin']['rr']/1000.,
                                         dpol_var=dpol_var,
                                         scatCoefDict=scat_coefs_dict,
                                         )
    del scat_coefs_dict
    del elev_temp, Tc_temp, content_temp, Fw_temp, concentration_temp
    
    return dpolDict_h



def dpol_var_from_scatcoefs(wavelength:float,
                            dpol_var:list,
                            scatCoefDict:dict,
                            ) -> dict[np.ndarray]:
    """Compute linear polarimetric variables."""
    temp_dict = {}
    if "Zhhlin" in dpol_var :
        temp_dict["Zhhlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*scatCoefDict['S22S22']
    if "Zvvlin" in dpol_var :
        temp_dict["Zvvlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*scatCoefDict['S11S11']
    if "Kdp" in dpol_var :
        temp_dict["Kdp"] = 180.*1e3/math.pi*wavelength*scatCoefDict['ReS22fmS11f']
    if "Rhohv" in dpol_var :
        temp_dict["S11S22"] = scatCoefDict['ReS22S11']**2+scatCoefDict['ImS22S11']**2
        temp_dict["S11S11xS22S22"] = scatCoefDict['S11S11'] * scatCoefDict['S22S22']
    return temp_dict



def compute_dpol_var(dpolDict:dict[np.ndarray]) -> dict[np.ndarray]:
    """Compute ZH, ZDR and RHOHV radar variables."""
    if 'Zh' in dpol2add :
        dpolDict["Zh"] = np.copy(dpolDict["Zhhlin"])
        dpolDict["Zh"][dpolDict["Zhhlin"]>0] = linear_to_dBZ(dpolDict["Zhhlin"][dpolDict["Zhhlin"]>0])
    if 'Zdr' in dpol2add :
        dpolDict["Zdr"] = np.copy(dpolDict["Zhhlin"])
        mask_Zdr = (dpolDict["Zhhlin"]>0) & (dpolDict["Zvvlin"]>0)
        dpolDict["Zdr"][mask_Zdr] = linear_to_dBZ((dpolDict["Zhhlin"]/dpolDict["Zvvlin"])[mask_Zdr])
    if 'Rhohv' in dpol2add :
        dpolDict["Rhohv"] = np.sqrt(np.divide(dpolDict["S11S22"], dpolDict["S11S11xS22S22"]))
    
    return dpolDict



def linear_to_dBZ(Z:np.ndarray) -> np.ndarray:
    """Linear to dBZ conversion.""" 
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Z[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp



def build_list_dpol_var_and_scatcoefs(dpol2add:list):
    """Depending on the variables user wants to add, creating lists of linear variables and scattering coefficients to compute."""
    dpol_var,scat_coefs = [],[]
    if 'Zh' in dpol2add :
        dpol_var += ['Zhhlin']
        scat_coefs += ['S22S22']
    if 'Zdr' in dpol2add :
        if 'Zhhlin' not in dpol_var : dpol_var += ['Zhhlin','Zvvlin']
        else : dpol_var += ['Zvvlin']
        scat_coefs += ['S11S11']
    if 'Kdp' in dpol2add :
        dpol_var += ['Kdp']
        scat_coefs += ['ReS22fmS11f']
    if 'Rhohv' in dpol2add :
        dpol_var += ['Rhohv','S11S22','S11S11xS22S22']
        if 'S11S11' not in scat_coefs : scat_coefs += ['S11S11']
        if 'S22S22' not in scat_coefs : scat_coefs += ['S22S22']
        scat_coefs += ['ReS22S11','ImS22S11']
    
    return dpol_var, scat_coefs