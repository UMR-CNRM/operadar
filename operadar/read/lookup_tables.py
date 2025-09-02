#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augros & davidcl        
"""

import sys
import time as tm
import numpy as np
import pandas as pd



def initialize_table_dictionary(dpol2add:list,
                                method:str='Tmatrix',
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
                                                  )
    for key in list_of_parameters + columns_to_loop_over + columns_to_retrieve :
        empty_table_dict[key] = {}
        
    return empty_table_dict, list_of_parameters, columns_to_retrieve



def retrieve_needed_columns(dpol2add:list,
                            scattering_method:str='Tmatrix',
                            ) -> list :
    """Depending on the variables the user wants to compute and the chosen scattering method,
    creating a list of the table's column names to extract."""
    
    table_columnNames = []
    
    if scattering_method == 'Tmatrix' or scattering_method == 'both' :
        if 'Zh' in dpol2add :
            table_columnNames += ['sighh']
        if 'Zdr' in dpol2add :
            table_columnNames += ['sighh','sigvv']
        if 'Kdp' in dpol2add :
            table_columnNames += ['kdp']
        if 'Rhohv' in dpol2add :
            table_columnNames += ['REdeltaco','IMdeltaco','sighh','sigvv']
        if 'Ah' in dpol2add :
            table_columnNames += ['Ah']
        if 'Av' in dpol2add :
            table_columnNames += ['Av']
            
    if scattering_method == 'Rayleigh' or scattering_method == 'both' :
        if 'Zh' in dpol2add :
            table_columnNames += ['sighhR']
        if 'Zdr' in dpol2add :
            table_columnNames += ['sighhR'',sigvvR']
        if 'Kdp' in dpol2add :
            table_columnNames += ['kdpR']
        if 'Rhohv' in dpol2add :
            table_columnNames += ['REdeltacoR','IMdeltacoR','sighhR','sigvvR']
        if 'Ah' in dpol2add :
            table_columnNames += ['AhR']
        if 'Av' in dpol2add :
            table_columnNames += ['AvR']    
    
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



def read_and_extract_tables_content(band:str,
                                    hydrometeors:list,
                                    moments:dict,
                                    scheme:str,
                                    dpol2add:list,
                                    path_table:str,
                                    verbose:bool=True,
                                    cloud_water_over:str='land',
                                    ) -> dict :
    """Extract min/step/max parameters and necessary columns in the table for later
    computation of the dual-pol variables.
    
    Args:
        band (str): radar band.
        hydrometeors (list): list of hydrometeors for which tables must be read.
        moments (dict): based on the dictionnary in the configuration file
        cloud_water_over (str): wether to read cloud water tables over sea or over land
        scheme (str): microphysics scheme.
        dpol2add (list): list of dual-polarization variables used to retrieve the relevant columns in the table
        path_table (str): table's directory path.
        verbose (bool): if True, print details about ongoing progress.

    Returns:
        dict : dictionary containing min/step/max values for multiple parameters
    """
    print("Reading tables for",band,"band")
    deb_timer = tm.time()
    micro_for_table = scheme[0:4]
    table_dict, parameters_to_retrieve, columns_to_retrieve = initialize_table_dictionary(dpol2add=dpol2add)
    
    for h in hydrometeors:
        if h == 'cc':
            hfile = cloud_water_species(hydrometeor=h, cloud_water_over=cloud_water_over)
        else :
            hfile = h
        nomfileCoefInt = f'{path_table}TmatCoefInt_{micro_for_table}_{moments[h]}M_{band}{hfile}'
        
        if verbose : print("\tReading min/step/max for",h)
        df_params = pd.read_csv(nomfileCoefInt, sep=";",nrows = 1)
        for value in parameters_to_retrieve :
            table_dict[value][h] = np.copy(df_params[value])[0]
        del df_params
        
        if verbose : print("\tRetrieving necessary columns in the table for",h)
        df_columns = pd.read_csv(nomfileCoefInt, sep=";",skiprows = [0, 1])
        df_columns = df_columns.astype(float)

        for columnName in columns_to_retrieve :
            if columnName == 'Fw' or columnName == 'N' :
                if moments[h]==1 :
                    table_dict['Fw'][h] = df_columns['P3'].to_numpy()
                elif moments[h]==2 :
                    table_dict['N'][h] = df_columns['P3'].to_numpy()
            else :
                table_dict[columnName][h] = df_columns[columnName].to_numpy()
        del df_columns
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    return table_dict
