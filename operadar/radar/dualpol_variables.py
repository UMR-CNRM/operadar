#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import time as tm
import numpy as np
from pathlib import Path
from pandas import Timestamp

from operadar.operadar_conf import save_netcdf_single_hydrometeor, dpol2add
from operadar.save.save_dpolvar import save_netcdf
from operadar.utils.masking import mask_hydrometeor
from operadar.read.lookup_tables import perform_nD_interpolation
from operadar.utils.formats_data import Fw_or_Nc
from operadar.utils.make_links import link_keys_with_available_hydrometeors



def compute_dualpol_variables(temperature:np.ndarray,
                              mask_precip_dist:np.ndarray,
                              elev:np.ndarray, Fw:np.ndarray,
                              contents:dict[np.ndarray],
                              concentrations:dict[np.ndarray],
                              tables_dict:dict,
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
        tables_dict (dict): dict containing the tables parameters and required columns.
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
    
    var2interpol = variables_to_interpolate(which_dpol=dpol2add)
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                         datatype='tables'
                                                         )
    print("Computation of",dpol2add,"for :") ; deb_timer = tm.time()
    
    initialize_dict = 0
    for h in hydrometeors :
        print('\t- hydrometeor :',h)
        
        # Mask single hydrometeor type
        mask_content = mask_hydrometeor(content=contents[h],
                                        expMmin=tables_dict['expMmin'][h]
                                        )
        mask_tot = (mask_precip_dist & mask_content) 

        dpolDict = compute_scatcoeffs_single_hydrometeor(hydrometeor=h,
                                                         variables_for_interpolation=var2interpol,
                                                         mask_tot=mask_tot,
                                                         Tc=temperature,
                                                         el=elev, Fw=Fw,
                                                         content_h=contents[h],
                                                         concentration_h=concentrations[h],
                                                         tables_dict=tables_dict,
                                                         hydrometeorMoments=hydrometeorMoments,
                                                    ) 
        
        # Addition of scattering coef for all hydrometeor
        if initialize_dict == 0 :
            fields2sum = {var:np.zeros(temperature.shape) for var in dpolDict.keys()}
            initialize_dict = 1
        for var in dpolDict.keys():
            fields2sum[var][mask_tot]+=dpolDict[var]
        
        # If saving single type, compute final polarimetric values
        if save_netcdf_single_hydrometeor :
            dpol_h = {var:np.zeros(temperature.shape) for var in dpolDict.keys()}
            for var in dpolDict.keys():
                dpol_h[var][mask_tot]=dpolDict[var]
                dpol_h[var][~mask_tot]= np.nan
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
        del dpolDict
        if append_in_fa : del concentrations[h],contents[h]
        
    for var in fields2sum.keys():
        fields2sum[var][~mask_precip_dist] = np.nan 
        
    # Dpol var calculation over the sum of scatering coefficients and linear Z
    fields2sum = compute_dpol_var(dpolDict=fields2sum)
    
    print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")    
    return fields2sum  



def compute_scatcoeffs_single_hydrometeor(hydrometeor:str,
                                          variables_for_interpolation:list,
                                          Tc:np.ndarray,
                                          el:np.ndarray,
                                          Fw:np.ndarray,
                                          content_h:np.ndarray,
                                          mask_tot:np.ndarray,
                                          concentration_h:np.ndarray,
                                          tables_dict:dict,
                                          hydrometeorMoments:dict[int],
                                        ) -> dict[np.ndarray]:
    """Compute radar scattering coefficients for a single hydrometeor class."""
    
    elev_temp=el[mask_tot]
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    content_temp=content_h[mask_tot]
    concentration_temp=concentration_h[mask_tot]
    
    # Define P3 : Nc (2 moments) or Fw (1 moment) 
    field_temp, colMin, colMax, colStep, colName = Fw_or_Nc(momentsDict=hydrometeorMoments,
                                                            hydrometeor=hydrometeor,
                                                            concentration=concentration_temp,
                                                            Fw=Fw_temp,
                                                            tables_dict=tables_dict,
                                                            )
    # Estimate for each grid point the scattering coefficients values based on the lookup tables
    # (return a dict containing 3D fields of scattering coefficients)
    fields3D_from_table = perform_nD_interpolation(tableDict=tables_dict,
                                   which_columns=variables_for_interpolation,
                                   hydrometeor=hydrometeor,
                                   colName=colName,
                                   colMin=colMin,
                                   colStep=colStep,
                                   colMax=colMax,
                                   el_temp=elev_temp,
                                   Tc_temp=Tc_temp,
                                   colTable=field_temp,
                                   M_temp=content_temp,
                                   )
    # Compute dualpol variables from scattering coefficients
    dpolDict_h = dpol_var_from_scatcoefs(wavelength=tables_dict['LAM'][hydrometeor],
                                         interpolated_from_table=fields3D_from_table,
                                         )
    del fields3D_from_table
    del elev_temp, Tc_temp, content_temp, Fw_temp, concentration_temp
    
    return dpolDict_h



def dpol_var_from_scatcoefs(wavelength:float,
                            interpolated_from_table:dict,
                            ) -> dict[np.ndarray]:
    """Compute linear polarimetric variables."""
    wavelength=wavelength/1000
    temp_dict = {}
    if "Zh" in dpol2add or "Zdr" in dpol2add :
        temp_dict["Zhhlin"]= ((1e3*wavelength)**4./(math.pi**5.*0.93))*interpolated_from_table['sighh']
    if "Zdr" in dpol2add :
        temp_dict["Zvvlin"]=((1e3*wavelength)**4./(math.pi**5.*0.93))*interpolated_from_table['sigvv']
    if "Kdp" in dpol2add :
        temp_dict["Kdp"] = interpolated_from_table['kdp']
    if "Rhohv" in dpol2add :
        temp_dict["numerator"] = interpolated_from_table['REdeltaco']**2+interpolated_from_table['IMdeltaco']**2
        temp_dict["denominator"] = interpolated_from_table['sigvv'] * interpolated_from_table['sighh']
    if "Ah" in dpol2add :
        temp_dict["Ah"] = interpolated_from_table['Ah']
    if "Av" in dpol2add :
        temp_dict["Av"] = interpolated_from_table['Av']
    return temp_dict



def compute_dpol_var(dpolDict:dict[np.ndarray]) -> dict[np.ndarray]:
    """Compute polarimetric variables."""
    finalDict = {}
    if 'Zh' in dpol2add :
        finalDict["Zh"] = np.copy(dpolDict["Zhhlin"])
        finalDict["Zh"][dpolDict["Zhhlin"]>0] = linear_to_dBZ(finalDict["Zh"][dpolDict["Zhhlin"]>0])
    if 'Zdr' in dpol2add :
        finalDict["Zdr"] = dpolDict["Zhhlin"]/dpolDict["Zvvlin"]
        mask_Zdr = (dpolDict["Zhhlin"]>0) & (dpolDict["Zvvlin"]>0)
        finalDict["Zdr"][mask_Zdr] = linear_to_dBZ(finalDict["Zdr"][mask_Zdr])
    if "Kdp" in dpol2add :
        finalDict["Kdp"] = np.copy(dpolDict["Kdp"])
    if 'Rhohv' in dpol2add :
        finalDict["Rhohv"] = 4*math.pi * np.sqrt(np.divide(dpolDict["numerator"], dpolDict["denominator"]))
    if "Ah" in dpol2add :
        finalDict["Ah"] = np.copy(dpolDict["Ah"])
    if "Av" in dpol2add :
        finalDict["Av"] = np.copy(dpolDict["Av"])
    
    return finalDict



def linear_to_dBZ(Z:np.ndarray) -> np.ndarray:
    """Linear to dBZ conversion.""" 
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Ztemp[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp



def variables_to_interpolate(which_dpol:list,scattering_method:str='Tmatrix'):
    """Depending on the variables user wants to add and the chosen scattering method, creating lists of linear variables and scattering coefficients to compute."""
    to_be_summed = []
    
    if scattering_method == 'Tmatrix' or scattering_method == 'both' :
        if 'Zh' in which_dpol :
            to_be_summed += ['sighh']
        if 'Zdr' in which_dpol :
            to_be_summed += ['sighh','sigvv']
        if 'Kdp' in which_dpol :
            to_be_summed += ['kdp']
        if 'Rhohv' in which_dpol :
            to_be_summed += ['REdeltaco','IMdeltaco','sighh','sigvv']
        if 'Ah' in which_dpol :
            to_be_summed += ['Ah']
        if 'Av' in which_dpol :
            to_be_summed += ['Av']
            
    if scattering_method == 'Rayleigh' or scattering_method == 'both' :
        if 'Zh' in which_dpol :
            to_be_summed += ['sighhR']
        if 'Zdr' in which_dpol :
            to_be_summed += ['sighhR','sigvvR']
        if 'Kdp' in which_dpol :
            to_be_summed += ['kdpR']
        if 'Rhohv' in which_dpol :
            to_be_summed += ['REdeltacoR','IMdeltacoR','sighhR','sigvvR']
        if 'Ah' in which_dpol :
            to_be_summed += ['AhR']
        if 'Av' in which_dpol :
            to_be_summed += ['AvR']
    
    return list(set(to_be_summed)) # return unique list of strings