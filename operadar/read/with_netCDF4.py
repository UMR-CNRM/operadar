#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: davidcl
"""

import sys
import numpy as np
from netCDF4 import Dataset
from operadar.utils.formats_data import get_lat_lon_from_subdomain



def check_variable_is_in_dataset(mnhFile:Dataset):
    necessary_variables = ['XHAT','YHAT','ZHAT','PABST','THT','RHOREFZ']
    variables_list = mnhFile.variables.keys()
    missing_var = []
    for var in necessary_variables :
        if var not in variables_list :
            missing_var += [var]

    if len(missing_var)>0 :
        print('_____________')
        print('/!\ ERROR /!\ Missing variables :',missing_var)
        sys.exit()



def get_geometry(mnhFile:Dataset,
                 real_case:bool,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 ):
    X = mnhFile.variables['XHAT'][:][j_min:j_max]
    Y = mnhFile.variables['YHAT'][:][i_min:i_max]
    ZHAT = mnhFile.variables['ZHAT'][:] # height above ground (m)
    if real_case :
        ZS = mnhFile.variables['ZS'][:][i_min:i_max,j_min:j_max] # altitude (m) of model surface (ground)
        LAT = mnhFile.variables['latitude'][:][i_min:i_max,j_min:j_max]
        LON = mnhFile.variables['longitude'][:][i_min:i_max,j_min:j_max]
    else :
        ZS = np.zeros((X.shape[0],Y.shape[0])) # is null for idealized cases
        LAT =float('nan')
        LON = float('nan')
    Z_3D = np.empty((ZHAT.shape[0],ZS.shape[0],ZS.shape[1]))
    for level in range(ZHAT.shape[0]):
        Z_3D[level,:,:]=ZS[:,:]+ZHAT[level]
    return X, Y, Z_3D, LAT, LON



def get_subdomain_indices(mnhFile:Dataset,
                          subDomain:list[float],
                          real_case:bool,
                          ):
    if real_case:
        LAT = mnhFile.variables['latitude'][:]
        LON = mnhFile.variables['longitude'][:]
        lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
        mask_zoom = ((LON>lon_min) &(LON<lon_max) & (LAT>lat_min) & (LAT<lat_max) )
        [ilon,jlat]=np.where(mask_zoom)
        i_min, i_max = np.nanmin(ilon), np.nanmax(ilon)
        j_min, j_max = np.nanmin(jlat), np.nanmax(jlat)  
    else :
        i_min, i_max = subDomain[0], subDomain[1]
        j_min, j_max = subDomain[2], subDomain[3]
    return i_min, i_max, j_min, j_max



def get_temperature_and_density(mnhFile:Dataset,
                                i_min:int, i_max:int,
                                j_min:int, j_max:int,
                                ):
    # Temperature
    p=mnhFile.variables['PABST'][0,:,:,:][:,i_min:i_max,j_min:j_max]
    p[np.where(p==999.)] = float('nan')
    Th=mnhFile.variables['THT'][0,:,:,:][:,i_min:i_max,j_min:j_max]
    Th[np.where(Th==999.)] = float('nan')
    temperature_celsius = Th*((p/100000)**(0.4/1.4))-273.15
    # Density
    rhodref = mnhFile.variables['RHOREFZ'][:]
    density_3D = np.ones(temperature_celsius.shape)
    IKE = temperature_celsius.shape[0]
    for k in range(IKE):    
        density_3D[k,:,:]=rhodref[k]
    return temperature_celsius, density_3D



def get_contents(mnhFile:Dataset,
                 hydrometeors:list,
                 name_var_hydro:dict,
                 temperature:np.ndarray,
                 rho3D:np.ndarray,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 )->dict[np.ndarray]:
    contents = {}
    for key in hydrometeors:
        contents[key] = np.empty(temperature.shape)
        contents[key] = mnhFile.variables[name_var_hydro[key]][0,:,:,:][:,i_min:i_max,j_min:j_max]*rho3D[:,:,:] # kg/kg of dry air
        contents[key][contents[key]==999.] = float('nan')
    return contents



def get_concentrations(mnhFile:Dataset,
                       microphysics_scheme:str,
                       hydrometeors:list,
                       temperature:np.ndarray,
                       rho3D:np.ndarray,
                       i_min:int, i_max:int,
                       j_min:int, j_max:int,):
    concentrations = {}
    for key in hydrometeors:
        concentrations[key] = np.zeros(temperature.shape)
    
    if microphysics_scheme[0:3]=="ICE":
        concentrations['ii'] = mnhFile.variables['CIT'][0,:,:,:][:,i_min:i_max,j_min:j_max]
        concentrations['ii'][concentrations['ii']==999.] = float('nan')
    if microphysics_scheme[0:3] =="LIM" :
        concentrations['rr'] = mnhFile.variables['CRAIN'][0,:,:,:][:,i_min:i_max,j_min:j_max] #former name: CRAINT
        concentrations['rr'][concentrations['rr']==999.]=float('nan')
        concentrations['ii'] = mnhFile.variables['CICE'][0,:,:,:][:,i_min:i_max,j_min:j_max] #former name: CICET
        concentrations['ii'][concentrations['ii']==999.]=float('nan')
    concentrations['rr']*=rho3D
    concentrations['ii']*=rho3D
    return concentrations
