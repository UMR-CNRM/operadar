#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import time as tm
import numpy as np
from pathlib import Path
from pandas import Timestamp

import operadar.operadar_conf as cf
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
                              X:np.ndarray,
                              Y:np.ndarray,
                              Z:np.ndarray,
                              lon:np.ndarray,
                              lat:np.ndarray,
                              date_time:Timestamp,
                              output_file_path:Path,
                              )-> dict[np.ndarray] :
    """Compute synthetic radar dual-polarimetrization variables for a given wavelength,
    microphysics, and mixed phase parametrization.

    Args:
        temperature (np.ndarray): 3D array
        mask_precip_dist (np.ndarray): 3D array of bool
        elev (np.ndarray): 3D array
        Fw (np.ndarray): 3D array
        contents (dict[np.ndarray]): dict of 3D array (one per hydrometeor)
        concentrations (dict[np.ndarray]): dict of 3D array (one per hydrometeor)
        tmatrix_param (dict): dict of Tmatrix parameters
        X (np.ndarray): 1D array of horizontal grid coordinates
        Y (np.ndarray): 1D array of horizontal grid coordinates
        Z (np.ndarray): vertical coordinates
        lon (np.ndarray): 2D array
        lat (np.ndarray): 2D array
        date_time (Timestamp): temporal variable
        output_file_path (Path): out file path

    Returns:
        dict[np.ndarray]: dictionnary containing all the dual-polarimetrization fields.
    """
    
    dpol_var = ["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.hydrometeors_moments,
                                                         datatype='tmatrix'
                                                         )
    print("Computation of",dpol_var,"for :") ; deb_timer = tm.time()
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
                                                       mask_tot=mask_tot,
                                                       Tc=temperature,
                                                       el=elev, Fw=Fw,
                                                       content_h=contents[h],
                                                       concentration_h=concentrations[h],
                                                       tmatrix_param=tmatrix_param,
                                                       final_dpolDict=dpol_var_dict,
                                                    )

        # Addition of scattering coef for all hydromet
        for var in dpol_var:
            dpol_var_dict[var][mask_tot]+=dpolDict[var]
            if cf.save_netcdf_single_hydrometeor :
                dpol_h = {var:np.zeros(temperature.shape) for var in dpol_var}
                dpol_h[var][mask_tot]=dpolDict[var]
                dpol_h[var][~mask_tot]= np.nan
        del dpolDict
        
        # If saving single type, compute final polarimetric values
        if cf.save_netcdf_single_hydrometeor :
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
        
    for var in dpol_var:
        dpol_var_dict[var][~mask_precip_dist] = np.nan 
        
    # Dpol var calculation over the sum of scatering coefficients and linear Z
    dpol_var_dict = compute_dpol_var(dpolDict=dpol_var_dict)
    
    print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")    
    return dpol_var_dict  



def compute_scatcoeffs_single_hydrometeor(hydrometeor:str,
                                        dpol_var:list,
                                        Tc:np.ndarray,
                                        el:np.ndarray,
                                        Fw:np.ndarray,
                                        content_h:np.ndarray,
                                        mask_tot:np.ndarray,
                                        concentration_h:np.ndarray,
                                        tmatrix_param:dict,
                                        final_dpolDict
                                        ) -> dict[np.ndarray]:
    """Compute radar scattering coefficients for a single hydrometeor class."""
    dpolDict_h = {var:np.zeros(Tc.shape)[mask_tot] for var in dpol_var}

    elev_temp=el[mask_tot]
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    content_temp=content_h[mask_tot]
    concentration_temp=concentration_h[mask_tot]
    
    # Define P3 : Nc (2 moments) or Fw (1 moment) hydrometeor_moment = cf.hydrometeors_moments[hydrometeor]
    field_temp, colMin, colMax, colStep, colName = select_Tmatrix_column(momentsDict=cf.hydrometeors_moments,
                                                                         hydrometeor=hydrometeor,
                                                                         concentration=concentration_temp,
                                                                         Fw=Fw_temp,
                                                                         tmatrix_param=tmatrix_param,
                                                                         )
    # Extract scattering coefficients for save_netcdf_single_hydrometeor
    [dpolDict_h["S11S11"], dpolDict_h["S22S22"],
     ReS22fmS11f, ReS22S11, ImS22S11] = get_scatcoef(tmatDict=tmatrix_param,
                                                     hydrometeor=hydrometeor,
                                                     colName_tmatrix=colName,
                                                     colMin=colMin,
                                                     colStep=colStep,
                                                     colMax=colMax,
                                                     el_temp=elev_temp,
                                                     Tc_temp=Tc_temp,
                                                     colTmat=field_temp,
                                                     M_temp=content_temp,
                                                     n_interpol=cf.n_interpol,
                                                     )
    # Compute dualpol variables from scattering coefficients
    dpolDict_h = dpol_var_from_scatcoefs(temp_dict=dpolDict_h,
                                         wavelength=tmatrix_param['LAMmin']['rr']/1000.,
                                         ReS22fmS11f=ReS22fmS11f, ReS22S11=ReS22S11,
                                         ImS22S11=ImS22S11)
    
    del ReS22fmS11f, ReS22S11, ImS22S11
    del elev_temp, Tc_temp, content_temp, Fw_temp, concentration_temp
    
    return dpolDict_h



def dpol_var_from_scatcoefs(temp_dict:dict[np.ndarray],
                            wavelength:float,
                            ReS22fmS11f, ReS22S11, ImS22S11,
                            ) -> dict[np.ndarray]:
    """Compute linear polarimetric variables."""
    temp_dict["Zhhlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*temp_dict["S22S22"]
    temp_dict["Zvvlin"]= 1e18*wavelength**4./(math.pi**5.*0.93)*4.*math.pi*temp_dict["S11S11"]
    temp_dict["Kdp"] = 180.*1e3/math.pi*wavelength*ReS22fmS11f
    temp_dict["S11S22"] = ReS22S11**2+ImS22S11**2
    
    return temp_dict



def compute_dpol_var(dpolDict:dict[np.ndarray]) -> dict[np.ndarray]:
    """Compute ZH, ZDR and RHOHV radar variables."""
    # ZH
    dpolDict["Zhh"] = np.copy(dpolDict["Zhhlin"])
    dpolDict["Zhh"][dpolDict["Zhhlin"]>0] = linear_to_dBZ(dpolDict["Zhhlin"][dpolDict["Zhhlin"]>0])
    # ZDR
    dpolDict["Zdr"] = np.copy(dpolDict["Zhhlin"])
    mask_Zdr = (dpolDict["Zhhlin"]>0) & (dpolDict["Zvvlin"]>0)
    dpolDict["Zdr"][mask_Zdr] = linear_to_dBZ((dpolDict["Zhhlin"]/dpolDict["Zvvlin"])[mask_Zdr])
    # RHOHV
    dpolDict["Rhohv"] = np.sqrt(np.divide(dpolDict["S11S22"], dpolDict["S11S11"]*dpolDict["S22S22"]))
    
    return dpolDict



def linear_to_dBZ(Z:np.ndarray) -> np.ndarray:
    """Linear to dBZ conversion.""" 
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Z[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp
 