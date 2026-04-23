#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augros & davidcl        
"""

import sys
import time as tm
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


def initialize_table_dictionary(dpol2add:list,
                                method:str='Tmatrix',
                                test_mode:bool=False,
                                ) -> tuple[dict, list[str], list[str]] :
    """Initializes the dictionary to store all the necessary parameters and columns of the tables"""
    empty_table_dict = {}
    list_of_parameters = ['LAM',
                          'ELEVmin', 'ELEVstep', 'ELEVmax',
                          'Tcmin', 'Tcstep', 'Tcmax',
                          'Fwmin', 'Fwstep', 'Fwmax', 
                          'expMmin', 'expMstep', 'expMmax',
                          'expCCmin', 'expCCstep', 'expCCmax',
                          ]
    columns_to_loop_over = ['Tc', 'ELEV', 'M', 'Fw', 'N']
    columns_to_retrieve = retrieve_needed_columns(dpol2add=dpol2add,
                                                  scattering_method=method,
                                                  test_mode=test_mode,
                                                  )
    for key in list_of_parameters + columns_to_loop_over + columns_to_retrieve :
        empty_table_dict[key] = {}
        
    return empty_table_dict, list_of_parameters, columns_to_retrieve

def retrieve_needed_columns(dpol2add: List[str],
                            scattering_method: str = 'Tmatrix',
                            test_mode: bool = False,
                            ) -> List[str]:
    """Depending on the variables the user wants to compute and the chosen scattering method,
    creating a list of the table's column names to extract."""
    
    table_columnNames = []
    
    # Check for Tmatrix or Both
    if scattering_method in ['Tmatrix', 'both']:
        if 'Zh' in dpol2add:
            table_columnNames.append('sighh')
            if test_mode: table_columnNames.append('zhh')
        if 'Zdr' in dpol2add:
            table_columnNames.extend(['sighh', 'sigvv'])
        if 'Kdp' in dpol2add:
            table_columnNames.append('kdp')
        if 'Rhohv' in dpol2add:
            table_columnNames.extend(['REdeltaco', 'IMdeltaco', 'sighh', 'sigvv'])
        if 'Ah' in dpol2add:
            table_columnNames.append('Ah')
        if 'Av' in dpol2add:
            table_columnNames.append('Av')
            
    # Check for Rayleigh or Both
    if scattering_method in ['Rayleigh', 'both']:
        if 'Zh' in dpol2add:
            table_columnNames.append('sighhR')
        if 'Zdr' in dpol2add:
            # FIX: Corrected syntax error (double quote)
            table_columnNames.extend(['sighhR', 'sigvvR']) 
        if 'Kdp' in dpol2add:
            table_columnNames.append('kdpR')
        if 'Rhohv' in dpol2add:
            table_columnNames.extend(['REdeltacoR', 'IMdeltacoR', 'sighhR', 'sigvvR'])
        if 'Ah' in dpol2add:
            table_columnNames.append('AhR')
        if 'Av' in dpol2add:
            table_columnNames.append('AvR')    
    
    return list(set(table_columnNames)) # return unique list of strings


def cloud_water_species(hydrometeor:str,
                        cloud_water_over:str,
                        ) -> str :
    """Handle the different name of the cloud water lookup tables (there is
    one for cloud water over sea and one for cloud water over land)"""
    if cloud_water_over == 'land' :
        hydrometeor = 'cl'
    elif cloud_water_over == 'sea' :
        hydrometeor = 'cs'
    else :
        print('_____________')
        print('/!\\ ERROR /!\\ :',cloud_water_over,'is not a valid option for cloud water')
        print('                can only be cloud water over "sea" or "land"')
        sys.exit()
    return hydrometeor



def read_and_extract_tables_content(band: str,
                                    hydrometeors: List[str],
                                    moments: Dict,
                                    scheme: str,
                                    dpol2add: List[str],
                                    path_table: str,
                                    verbose: bool = True,
                                    cloud_water_over: str = 'land',
                                    test_interpolation: bool = False,
                                    ) -> Tuple[dict, List[str], List[str]]:
    """Extract min/step/max parameters and necessary columns in the table for later
    computation of the dual-pol variables.

    Args:
        band (str): radar band.
        hydrometeors (list): list of hydrometeors for which tables must be read.
        moments (dict): based on the dictionnary in the configuration file
        scheme (str): microphysics scheme.
        dpol2add (list): list of dual-polarization variables used to retrieve the relevant columns in the table
        path_table (str): table's directory path.
        verbose (bool): if True, print details about ongoing progress.
        cloud_water_over (str): wether to read cloud water tables over sea or over land
        test_interpolation: bool

    Returns table_dict, parameters_to_retrieve, columns_to_retrieve
    """
    print("Reading tables for", band, "band")
    deb_timer = tm.time()
    micro_for_table = scheme[0:4]
    
    table_dict, parameters_to_retrieve, columns_to_retrieve = initialize_table_dictionary(
        dpol2add=dpol2add,
        test_mode=test_interpolation,
    )
    
    # If testing interpolation, ensure we have the grid coordinates in the dict
    if test_interpolation: 
        columns_to_retrieve.extend(['Tc', 'ELEV', 'M', 'Fw', 'N'])
    
    for h in hydrometeors:
        if h == 'cc':
            hfile = cloud_water_species(hydrometeor=h, cloud_water_over=cloud_water_over)
        else:
            hfile = h
        
        nomfileCoefInt = f'{path_table}TmatCoefInt_{micro_for_table}_{moments[h]}M_{band}{hfile}'
        
        if verbose: print("\tReading min/step/max for", h, "in", nomfileCoefInt)
        df_params = pd.read_csv(nomfileCoefInt, sep=";", nrows=1)
        
        for value in parameters_to_retrieve:
            # Ensure we capture the value safely
            try:
                table_dict[value][h] = np.copy(df_params[value])[0]
            except KeyError:
                print(f"Warning: Parameter '{value}' not found in header for {h}")
                table_dict[value][h] = np.nan
        
        del df_params
        
        if verbose: print("\tRetrieving necessary columns in the table for", h)
        # Read data
        df_columns = pd.read_csv(
            nomfileCoefInt, 
            sep=";", 
            skiprows=[0, 1], 
            na_values=["Infinity", "-Infinity"],
            low_memory=False
        )
        df_columns = df_columns.astype(float)

        for columnName in columns_to_retrieve:
            # ---- Logic for P3 (Fw or N) ----
            if columnName == 'Fw' or columnName == 'N': 
                # Note: Ensure moments[h] logic matches your table structure
                if moments[h] == 1:
                    if h == 'ii': 
                        table_dict['N'][h] = df_columns['P3'].to_numpy()
                    else: 
                        table_dict['Fw'][h] = df_columns['P3'].to_numpy()
                elif moments[h] == 2: 
                    table_dict['N'][h] = df_columns['P3'].to_numpy()
            # ------------------------------- 
            else:
                # Standard columns
                if columnName in df_columns.columns:
                    table_dict[columnName][h] = df_columns[columnName].to_numpy()
                else:
                    # Fill with zeros or NaN if column missing
                    table_dict[columnName][h] = np.zeros(len(df_columns))
        
        del df_columns
    
    print("\t--> Done in", round(tm.time() - deb_timer, 2), "seconds")
    return table_dict, parameters_to_retrieve, columns_to_retrieve




# def variables_to_interpolate(which_dpol:list,scattering_method:str='Tmatrix'):
#     """Depending on the variables user wants to add and the chosen scattering method, creating lists of linear variables and scattering coefficients to compute."""
#     to_be_summed = []
    
#     if scattering_method == 'Tmatrix' or scattering_method == 'both' :
#         if 'Zh' in which_dpol :
#             to_be_summed += ['sighh']
#         if 'Zdr' in which_dpol :
#             to_be_summed += ['sighh','sigvv']
#         if 'Kdp' in which_dpol :
#             to_be_summed += ['kdp']
#         if 'Rhohv' in which_dpol :
#             to_be_summed += ['REdeltaco','IMdeltaco','sighh','sigvv']
#         if 'Ah' in which_dpol :
#             to_be_summed += ['Ah']
#         if 'Av' in which_dpol :
#             to_be_summed += ['Av']
            
#     if scattering_method == 'Rayleigh' or scattering_method == 'both' :
#         if 'Zh' in which_dpol :
#             to_be_summed += ['sighhR']
#         if 'Zdr' in which_dpol :
#             to_be_summed += ['sighhR','sigvvR']
#         if 'Kdp' in which_dpol :
#             to_be_summed += ['kdpR']
#         if 'Rhohv' in which_dpol :
#             to_be_summed += ['REdeltacoR','IMdeltacoR','sighhR','sigvvR']
#         if 'Ah' in which_dpol :
#             to_be_summed += ['AhR']
#         if 'Av' in which_dpol :
#             to_be_summed += ['AvR']
    
#     return list(set(to_be_summed)) # return unique list of strings
