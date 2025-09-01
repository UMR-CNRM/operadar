#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from hmac import new
import math
import time as tm
import numpy as np
from pathlib import Path
from pandas import Timestamp

from operadar.utils.formats_data import Fw_or_Nc
from operadar.save.save_dpolvar import save_netcdf
from operadar.utils.masking import mask_hydrometeor
from operadar.read.lookup_tables import retrieve_needed_columns
from operadar.utils.make_links import link_keys_with_available_hydrometeors



def compute_dualpol_variables(temperature:np.ndarray,
                              mask_precip_dist:np.ndarray,
                              elev:np.ndarray,
                              Fw:np.ndarray,
                              contents:dict[str,np.ndarray],
                              concentrations:dict[str,np.ndarray],
                              tables_dict:dict,
                              X:np.ndarray,
                              Y:np.ndarray,
                              Z:np.ndarray,
                              lon:np.ndarray,
                              lat:np.ndarray,
                              date_time:Timestamp,
                              output_file_path:Path,
                              append_in_fa:bool,
                              config,
                              )-> dict[str,np.ndarray] :
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
        
    var2interpol = variables_to_interpolate(which_dpol=config.dpol2add)
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=config.hydrometeors_moments,
                                                         datatype='tables'
                                                         )
    print("Computation of",config.dpol2add,"for :") ; deb_timer = tm.time()
    
    initialize_dict = 0
    for h in hydrometeors :
        print('\t- hydrometeor :',h)
        
        # Mask single hydrometeor type
        mask_content = mask_hydrometeor(content = contents[h],
                                        expMmin = tables_dict['expMmin'][h]
                                        )
        mask_tot = (mask_precip_dist & mask_content) 

        dpolDict = compute_scatcoeffs_single_hydrometeor(hydrometeor=h,
                                                         dpol2add=config.dpol2add,
                                                         variables_for_interpolation=var2interpol,
                                                         mask_tot=mask_tot,
                                                         Tc=temperature,
                                                         el=elev, Fw=Fw,
                                                         content_h=contents[h],
                                                         concentration_h=concentrations[h],
                                                         tables_dict=tables_dict,
                                                         hydrometeorMoments=config.hydrometeors_moments,
                                                    ) 
        
        # Addition of scattering coef for all hydrometeor
        if initialize_dict == 0 :
            fields2sum = {var:np.zeros(temperature.shape) for var in dpolDict.keys()} 
            initialize_dict = 1
        for var in dpolDict.keys():
            fields2sum[var][mask_tot]+=dpolDict[var]
        
        # If saving single type, compute final polarimetric values
        if config.save_netcdf_single_hydrometeor :
            dpol_h = {var:np.zeros(temperature.shape) for var in dpolDict.keys()}
            for var in dpolDict.keys():
                dpol_h[var][mask_tot]=dpolDict[var]
                dpol_h[var][~mask_tot]= np.nan
            dpol_h = compute_dpol_var(dpolDict=dpol_h,dpol2add=config.dpol2add)
            
            parent_directory = output_file_path.parent
            new_filename=Path(f"{output_file_path.stem}_{h}")
            path_single_netcdf = parent_directory.joinpath(new_filename)
            save_netcdf(X=X, Y=Y, Z=Z,
                        lat=lat, lon=lon,
                        datetime=date_time,
                        dpolDict=dpol_h,
                        contentsDict={h:contents[h]},
                        concentrationsDict={h:concentrations[h]},
                        temperature=temperature,
                        outfile=path_single_netcdf,
                        config=config,
                        )
                 
            del dpol_h   
        del dpolDict
        if append_in_fa : del concentrations[h],contents[h]
        
    for var in fields2sum.keys():
        fields2sum[var][~mask_precip_dist] = np.nan 
        
    # Dpol var calculation over the sum of scatering coefficients and linear Z
    fields2sum = compute_dpol_var(dpolDict=fields2sum,dpol2add=config.dpol2add)
    
    print("\t--> Done in",round(tm.time() - deb_timer,2),"seconds")    
    return fields2sum  



def compute_scatcoeffs_single_hydrometeor(hydrometeor:str,
                                          dpol2add:list,
                                          variables_for_interpolation:list,
                                          Tc:np.ndarray,
                                          el:np.ndarray,
                                          Fw:np.ndarray,
                                          content_h:np.ndarray,
                                          mask_tot:np.ndarray,
                                          concentration_h:np.ndarray,
                                          tables_dict:dict,
                                          hydrometeorMoments:dict[str,int],
                                        ) -> dict[str,np.ndarray]:
    """Compute radar scattering coefficients for a single hydrometeor class."""

    elev_temp=el[mask_tot]
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    content_temp=content_h[mask_tot]
    concentration_temp=concentration_h[mask_tot]
    
    # Define P3 : Nc (2 moments) or Fw (1 moment)
    field_temp, colName = Fw_or_Nc(momentsDict=hydrometeorMoments,
                                   hydrometeor=hydrometeor,
                                   concentration=concentration_temp,
                                   Fw=Fw_temp,
                                   )
        
    # Estimate for each grid point the scattering coefficients values based on the lookup tables
    # (return a dict containing 3D fields of scattering coefficients)
    fields3D_from_table = perform_nD_interpolation(tableDict=tables_dict,
                                                    hydrometeor=hydrometeor,
                                                    colName=colName,
                                                    elev=elev_temp,
                                                    T=Tc_temp,
                                                    P3=field_temp,
                                                    content=content_temp,
                                                    dpol2add=dpol2add,
                                                    )
    # Compute dualpol variables from scattering coefficients
    dpolDict_h = dpol_var_from_scatcoefs(wavelength=tables_dict['LAM'][hydrometeor],
                                         interpolated_from_table=fields3D_from_table,
                                         dpol2add=dpol2add,
                                         )
    del fields3D_from_table
    del elev_temp, Tc_temp, content_temp, Fw_temp, concentration_temp
    
    return dpolDict_h



def dpol_var_from_scatcoefs(wavelength:float,
                            interpolated_from_table:dict,
                            dpol2add:list,
                            ) -> dict[str,np.ndarray]:
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



def compute_dpol_var(dpolDict:dict[str,np.ndarray],dpol2add:list) -> dict[str,np.ndarray]:
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



def perform_nD_interpolation(tableDict:dict,
                             hydrometeor:str,
                             colName:str,
                             elev:np.ndarray,
                             T:np.ndarray,
                             P3:np.ndarray,
                             content:np.ndarray,
                             dpol2add:list[str],
                             ) -> dict :
    """Construct 3D fields of scattering coefficients, based on the tables. The interpolation is
    currently performed with 4 fields (elevation, temperature, P3 (=Fw or Nc), and content).

    Args:
        tableDict (dict): dictionnary containing the necessary columns of the table to perform interpolation
        hydrometeor (str): hydrometeor type
        colName (str): P3 name (either Fw or Nc)
        elev (np.ndarray): angle elevation field
        T (np.ndarray): temperature (°C) field
        P3 (np.ndarray): liquid water fraction or number concentration field 
        content (np.ndarray): hydrometeor content field of values (kg/m3)
        dpol2add (list[str]): used to select the appropriate columns in the lookup table

    Returns:
        dict: dictionnary containing fields of interpolated scattering coefficients over the grid
    """
    
    M = np.copy(content)*0.0-100
    M[content>0] = np.log10(content[content>0])
    
    P = np.copy(P3)
    if colName == 'Fw' :
        P_min, P_max, step_P = tableDict['Fwmin'][hydrometeor], tableDict['Fwmax'][hydrometeor],tableDict['Fwstep'][hydrometeor]
    elif colName == 'Nc' :
        P_min, P_max, step_P = tableDict['expCCmin'][hydrometeor], tableDict['expCCmax'][hydrometeor],tableDict['expCCstep'][hydrometeor]
        P[P3>0]=np.log10(P3[P3>0])
        P[P3<=0]=P_min

    T_min, T_max, step_T = tableDict['Tcmin'][hydrometeor], tableDict['Tcmax'][hydrometeor],tableDict['Tcstep'][hydrometeor]
    elev_min, elev_max, step_elev = tableDict['ELEVmin'][hydrometeor], tableDict['ELEVmax'][hydrometeor],tableDict['ELEVstep'][hydrometeor]
    M_min, M_max, step_M = tableDict['expMmin'][hydrometeor], tableDict['expMmax'][hydrometeor],tableDict['expMstep'][hydrometeor]
    
    # Nombre de pas (important pour le calcul de l'index)
    n_T = int((T_max - T_min) / step_T) + 1
    n_elev = int((elev_max - elev_min) / step_elev) + 1
    n_P = int((P_max - P_min) / step_P) + 1
    n_M = int((M_max - M_min) / step_M) + 1

    # Bornes inférieures
    T_inf, T_sup = get_bounds(T, T_min, T_max, step_T)
    elev_inf, elev_sup = get_bounds(elev, elev_min, elev_max, step_elev)
    P_inf, P_sup = get_bounds(P, P_min, P_max, step_P)
    M_inf, M_sup = get_bounds(M, M_min, M_max, step_M)
    
    # Index des points inf et sup
    idx_inf = get_linear_index(T_inf, elev_inf, P_inf, M_inf,
                               T_min, elev_min, P_min, M_min,
                               step_T, step_elev, step_P, step_M,
                               n_T, n_elev, n_P, n_M)

    idx_sup = get_linear_index(T_sup, elev_sup, P_sup, M_sup,
                               T_min, elev_min, P_min, M_min,
                               step_T, step_elev, step_P, step_M,
                               n_T, n_elev, n_P, n_M)
    
    # Coefficient alpha sur la diagonale 4D
    dx = np.stack([T - T_inf, elev - elev_inf, P - P_inf, M - M_inf], axis=0)
    d  = np.stack([step_T, step_elev, step_P, step_M], axis=0)

    # alpha = projection sur la diagonale
    norm2 = np.sum(d**2)
    broadcast_shape = M.shape
    alpha = np.sum(dx * d.reshape(-1,*([1]*len(broadcast_shape))), axis=0) / norm2
    
    scatCoefsDict = {}
    scatCoef_columns = retrieve_needed_columns(dpol2add=dpol2add)
    for column in scatCoef_columns:
        val_inf = tableDict[column][hydrometeor][idx_inf]
        val_sup = tableDict[column][hydrometeor][idx_sup]
        scatCoefsDict[column] = (1 - alpha) * val_inf + alpha * val_sup

    return scatCoefsDict 



def get_linear_index(T:np.ndarray, E:np.ndarray, P:np.ndarray, M:np.ndarray,
                     T_min:float, E_min:float, P_min:float, M_min:float,
                     step_T:float, step_E:float, step_P:float, step_M:float,
                     n_T:float, n_E:float, n_P:float, n_M:float,
                     ) -> np.ndarray :
    """Compute linear index into a flattened 4D grid (T, E, P, M).

    Args:
        T (np.ndarray): Temperature values.
        E (np.ndarray): Angle elevation field values.
        P (np.ndarray): Concentration or liquid water fraction field values.
        M (np.ndarray): Content values.
        T,E,P,M _min (float): Minimum T, E, P or M value in the grid.
        step_ T,E,P,M (float): Step for T, E, P or M
        n_ T,E,P,M (int): Number of T, E, P or M points in the grid.

    Returns:
        np.ndarray: Array of the indices to use in the lookup table.
    """
    iT = np.floor((T - T_min) / step_T).astype(int)
    iE = np.floor((E - E_min) / step_E).astype(int)
    iP = np.floor((P - P_min) / step_P).astype(int)
    iM = np.floor((M - M_min) / step_M).astype(int)

    return (
        iT * n_E * n_P * n_M +
        iE * n_P * n_M +
        iP * n_M +
        iM
    )



def get_bounds(val:np.ndarray,
               val_min:float,
               val_max:float,
               step:float,
               ) -> tuple[np.ndarray, np.ndarray]:
    """Compute lower and upper bounds from a regular grid.

    Args:
        val (np.ndarray): Input array of values.
        val_min (float): Minimum value in the reference grid.
        val_max (float): Maximum value in the reference grid.
        step (float): Step size between grid points.

    Returns:
        tuple[np.ndarray, np.ndarray]: Lower and upper bounds arrays (same shape as input).
    """

    val_inf = np.floor((val - val_min) / step) * step + val_min
    val_sup = val_inf + step
    val_sup = np.where(val == val_inf, val_inf, val_sup)
    
    # Ensure the possible values stay in the range of allowed values for the given field
    val_inf = np.clip(val_inf, val_min, val_max)
    val_sup = np.clip(val_sup, val_min, val_max)
    return val_inf, val_sup
