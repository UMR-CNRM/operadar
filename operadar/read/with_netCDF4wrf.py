#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: davidcl
"""
"""
File modified by Eulalia Busquets to read netCDF from WRF at 2/6/2026
Main changes: 
- Names of variables: e.g. XLAT, XLONG for latitude and longitude
- Time coordinates: Variables are called as wrfFile.variables['P'][:,:,:] instead of wrfFile.variables['P'][0,:,:,:] (a single value in the Time dimension)

Pending: Incroporate the Thompson scheme in the get number concentrations
!! CAUTION !! - Regarding the WRF outputs...
- WRF need to save the 'rho' variable in the output. 
- A preprocess is needed to select a single time value, as well as adding the 'Tc' and 'z' variables (computes through wrf.getvar function of wrf python library)
"""
import sys
import numpy as np
from netCDF4 import Dataset
from typing import Sequence
from operadar.utils.formats_data import get_lat_lon_from_subdomain



def check_variable_is_in_dataset(wrfFile:Dataset):
    necessary_variables = ['XLAT','XLONG','tc','z','RHO','P','PB']
    variables_list = wrfFile.variables.keys()
    missing_var = []
    for var in necessary_variables :
        if var not in variables_list :
            missing_var += [var]

    if len(missing_var)>0 :
        print('_____________')
        print('/!\\ ERROR /!\\ Missing variables :',missing_var)
        sys.exit()



def get_geometry(wrfFile:Dataset,
                 real_case:bool,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 ):
    X = wrfFile.variables['west_east'][:][j_min:j_max]
    Y = wrfFile.variables['south_north'][:][i_min:i_max]
    Z = wrfFile.variables['z'][:,i_min:i_max,j_min:j_max] # height above sea level (m)
    if real_case :
        ZS = wrfFile.variables['HGT'][:][i_min:i_max,j_min:j_max] # altitude (m) of model surface (ground), but not used
        LAT = wrfFile.variables['XLAT'][:][i_min:i_max,j_min:j_max]
        LON = wrfFile.variables['XLONG'][:][i_min:i_max,j_min:j_max]
    else :
        ZS = np.zeros((X.shape[0],Y.shape[0])) # is null for idealized cases
        LAT =float('nan')
        LON = float('nan')
    Z_3D = Z
    return X, Y, Z_3D, LAT, LON



def get_subdomain_indices(wrfFile:Dataset,
                          subDomain:Sequence[float]|Sequence[int],
                          real_case:bool,
                          ) -> tuple[int,int,int,int]:
    if real_case:
        LAT = wrfFile.variables['XLAT'][:]
        LON = wrfFile.variables['XLONG'][:]
        lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
        mask_zoom = ((LON>lon_min) &(LON<lon_max) & (LAT>lat_min) & (LAT<lat_max) )
        [ilon,jlat]=np.where(mask_zoom)
        i_min, i_max = np.nanmin(ilon), np.nanmax(ilon)
        j_min, j_max = np.nanmin(jlat), np.nanmax(jlat)  
    else :
        i_min, i_max = subDomain[0], subDomain[1]
        j_min, j_max = subDomain[2], subDomain[3]
    return i_min, i_max, j_min, j_max



def get_pressure_temperature_density(wrfFile:Dataset,
                                i_min:int, i_max:int,
                                j_min:int, j_max:int,
                                ):
    # Pressure
    p= (wrfFile.variables['P'][:,:,:][:,i_min:i_max,j_min:j_max] + #ChangeEB - As the simulation only have 1 time, change [var][0,:,:,:] for [var][:,:,:]
        wrfFile.variables['PB'][:,:,:][:,i_min:i_max,j_min:j_max])
    # Temperature
    temperature_celsius=wrfFile.variables['tc'][:,:,:][:,i_min:i_max,j_min:j_max]
    # Density
    density_3D = wrfFile.variables['RHO'][:,:,:][:,i_min:i_max,j_min:j_max]
    return p, temperature_celsius, density_3D



def get_contents(wrfFile:Dataset,
                 hydrometeors:list,
                 name_var_hydro:dict,
                 temperature:np.ndarray,
                 rho3D:np.ndarray,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 )->tuple[dict[str,np.ndarray],np.ndarray]: #ChangeEB - I changed tuple() for tuple[]
    contents = {}
    for key in hydrometeors:
        contents[key] = np.empty(temperature.shape)
        contents[key] = wrfFile.variables[name_var_hydro[key]][:,:,:][:,i_min:i_max,j_min:j_max]*rho3D[:,:,:] # kg/kg of dry air
        contents[key][contents[key]==999.] = float('nan')
    qv=np.empty(temperature.shape)
    wrfFile.variables[name_var_hydro['vv']][:,:,:][:,i_min:i_max,j_min:j_max] 
    return contents,qv



def get_concentrations(wrfFile:Dataset, #ChangeEB_Pending: I need to implement here the Thompson Microphysics scheme...
                                        #ChangeEB_Pending: I also need to do this a function of the Microphysical Moments instead of schemes...
                       microphysics_scheme:str,
                       hydrometeors:list,
                       temperature:np.ndarray,
                       rho3D:np.ndarray,
                       i_min:int, i_max:int,
                       j_min:int, j_max:int,):
    concentrations = {}
    for key in hydrometeors:
        concentrations[key] = np.zeros(temperature.shape)
    
    if microphysics_scheme[0:2]=="IC": #ChangeEB_Pending: Here I will have to add the Thompson scheme...
        concentrations['ii'] = wrfFile.variables['QNICE'][:,:,:][:,i_min:i_max,j_min:j_max]
        concentrations['ii'][concentrations['ii']==999.] = float('nan')
    if microphysics_scheme[0:3] =="LIM" :
        concentrations['rr'] = wrfFile.variables['QNRAIN'][0,:,:,:][:,i_min:i_max,j_min:j_max] #former name: QNRAIN
        concentrations['rr'][concentrations['rr']==999.]=float('nan')
        concentrations['ii'] = wrfFile.variables['QNICE'][0,:,:,:][:,i_min:i_max,j_min:j_max] #former name: QNICE
        concentrations['ii'][concentrations['ii']==999.]=float('nan')
    if microphysics_scheme[0:4] == "THOM" : #If I have the Thompson MP scheme, I get the rain and ice number concentration
        concentrations['rr'] = wrfFile.variables['QNRAIN'][:,:,:][:,i_min:i_max,j_min:j_max] #former name: QNRAIN
        concentrations['rr'][concentrations['rr']==999.]=float('nan')
        concentrations['ii'] = wrfFile.variables['QNICE'][:,:,:][:,i_min:i_max,j_min:j_max] #former name: QNICE
        concentrations['ii'][concentrations['ii']==999.]=float('nan')
        #For cloud in Thompson Scheme Nc is set be 100cm-3 = 10**8m-3. I need to put this here to multiply by the density
        concentrations['cc'] =  np.full(temperature.shape, 10**8)
        concentrations['cc'][concentrations['ii']==999.]=float('nan')
                        
    concentrations['rr']*=rho3D
    concentrations['ii']*=rho3D
    concentrations['cc']*=rho3D
    #Prevent hydrometeor concentration from being negative 
    for key in hydrometeors:
        concentrations[key] = np.where(concentrations[key]<0, 0, concentrations[key])
    
    return concentrations
