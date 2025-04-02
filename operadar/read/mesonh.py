#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros & davidcl
"""

import sys
import time as tm
import numpy as np
from netCDF4 import Dataset
from operadar.utils.make_links import link_keys_with_available_hydrometeors, link_varname_with_mesonh_name
from operadar.utils.formats_data import get_lat_lon_from_subdomain


#============== Read MesoNH variables ===============
"""
Read MesoNH 3D variables in ncfile: pressure, temperature, hydrometeor contents 
"""
def read_mesonh(filePath: str,micro: str,subDomain:list[float]|None,
                hydrometeorMoments: dict[int],real_case: bool,verbose:bool,
               )-> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,dict[np.ndarray],dict[np.ndarray],np.ndarray]:
    """Read and extract data from an MesoNH.nc file

    Args:
        filePath (str): input file path.
        micro (str): microphysics scheme.
        subDomain (list[float] | None): either a list of 4 float or None.
        hydrometeorMoments (dict[int]): dictionary of form {'hydrometeor key':number of moments}.
        real_case (bool): if False, consider it is an idealized case.
        verbose (bool): show more messages to the user.

    Returns:
        X (ndarray): 1D horizontal coordinates in m
        Y (ndarray): 1D horizontal coordinates in m
        Z (ndarray): 1D array of vertical coordinates in model pressure levels
        LON (ndarray): 2D array of longitude coordinates
        LAT (ndarray): 2D array of latitude coordinates
        M (dict[ndarray]): dictionary of 3D contents for each hydrometeor 
        Nc (dict[ndarray]): dictionary of 3D number concentrations for each hydrometeor
        Tc (ndarray) : 3D temperature in Celsius
    """

    print("\tMesoNH .nc file:",filePath)
    if verbose : deb=tm.time()
    mnh_file = Dataset(filePath,'r')
    if verbose : print('\t\tLoaded file in',round(tm.time()-deb,3),'seconds'); deb=tm.time()

    X, Y, Z = get_geometry(mnhFile=mnh_file)
    if verbose : print('\t\tGot geometry (X,Y,Z) in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    if real_case:
        LAT = mnh_file.variables['latitude'][:]
        LON = mnh_file.variables['longitude'][:]
        if subDomain != None :
            lonmin, lonmax, latmin, latmax = extract_subdomain(mnhFile=mnh_file,
                                                               longitude_field=LON,
                                                               latitude_field=LAT,
                                                               subDomain=subDomain,
                                                               )
            if verbose : print('\t\tSubdomain extracted in',round(tm.time()-deb,6),'seconds'); deb=tm.time()
        else :
            if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
    else:
        LAT=float('nan')
        LON = float('nan')
        if verbose : print('\t\tIdealized case, will use all the points.'); deb=tm.time()
    
    Tc, rho3D = get_temperature_and_density(mnhFile=mnh_file)
    if verbose : print('\t\tGot temperature and dry air density in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments, datatype='model')
    name_hydro = link_varname_with_mesonh_name()
    
    M = get_contents(mnhFile=mnh_file,
                     hydrometeors=hydromet_list,
                     name_var_hydro=name_hydro,
                     temperature=Tc,
                     rho=rho3D,
                     )
    if verbose : print('\t\tGot 3D contents in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    Nc = get_concentrations(mnhFile=mnh_file,
                            microphysics_scheme=micro,
                            hydrometeors=hydromet_list,
                            temperature=Tc,
                            rho=rho3D,
                            )
    if verbose : print('\t\tGot 3D number concentrations in',round(tm.time()-deb,6),'seconds');deb=tm.time()

    # Cropping data if necessary (in a function later ?)
    if real_case and subDomain != None :
        M_zoom,Nc_zoom={},{}
        for key in hydromet_list:
            M_zoom[key]=M[key][:,lonmin:lonmax,latmin:latmax]
            Nc_zoom[key]=Nc[key][:,lonmin:lonmax,latmin:latmax]
        
        Tc_zoom=Tc[:,lonmin:lonmax,latmin:latmax]
        LAT_zoom=LAT[lonmin:lonmax,latmin:latmax]
        LON_zoom=LON[lonmin:lonmax,latmin:latmax]
        X_zoom=X[latmin:latmax]
        Y_zoom=Y[lonmin:lonmax]
        Z_zoom=Z[:,lonmin:lonmax,latmin:latmax]
        
        return X_zoom, Y_zoom, Z_zoom, LON_zoom, LAT_zoom, M_zoom, Nc_zoom, Tc_zoom 
    else:
        return X, Y, Z, LON, LAT, M, Nc, Tc



def get_geometry(mnhFile:Dataset):
    X = mnhFile.variables['XHAT'][:]
    Y = mnhFile.variables['YHAT'][:]
    ZHAT = mnhFile.variables['ZHAT'][:] # height above ground (m)
    ZS = mnhFile.variables['ZS'][:] # altitude (m) of model surface (ground)
    Z = np.empty((ZHAT.shape[0],ZS.shape[0],ZS.shape[1]))
    for level in range(ZHAT.shape[0]):
        Z[level,:,:]=ZS[:,:]+ZHAT[level]
    return X, Y, Z



def extract_subdomain(mnhFile:Dataset,
                      longitude_field:np.ndarray,
                      latitude_field:np.ndarray,
                      subDomain:list[float],
                      ):
    lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
    mask__zoom = ((longitude_field>lon_min)
                 &(longitude_field<lon_max)
                 & (latitude_field>lat_min)
                 & (latitude_field<lat_max)
                 )
    [ilon,jlat]=np.where(mask__zoom)
    i_lonmin,i_lonmax=np.nanmin(ilon),np.nanmax(ilon)
    j_latmin,j_latmax=np.nanmin(jlat),np.nanmax(jlat)
    return i_lonmin, i_lonmax, j_latmin, j_latmax



def get_temperature_and_density(mnhFile:Dataset):
    # Temperature
    p=mnhFile.variables['PABST'][0,:,:,:]
    p[np.where(p==999.)]=float('nan')
    Th=mnhFile.variables['THT'][0,:,:,:]
    Th[np.where(Th==999.)]=float('nan')
    Tc=Th*((p/100000)**(0.4/1.4))-273.15
    # Density
    rhodref = mnhFile.variables['RHOREFZ'][:]
    rho3D = np.ones(Tc.shape)
    IKE = Tc.shape[0]
    for k in range(IKE):    
        rho3D[k,:,:]=rhodref[k]
    return Tc, rho3D



def get_contents(mnhFile:Dataset,
                 hydrometeors:list,
                 name_var_hydro:dict,
                 temperature:np.ndarray,
                 rho:np.ndarray,
                 )->dict[np.ndarray]:
    contents = {}
    for key in hydrometeors:
        contents[key] = np.empty(temperature.shape)
        contents[key] = mnhFile.variables[name_var_hydro[key]][0,:,:,:]*rho[:,:,:] # kg/kg of dry air
        contents[key][contents[key]==999.] = float('nan')
    return contents



def get_concentrations(mnhFile:Dataset,
                       microphysics_scheme:str,
                       hydrometeors:list,
                       temperature:np.ndarray,
                       rho:np.ndarray,):
    concentrations = {}
    for key in hydrometeors:
        concentrations[key] = np.zeros(temperature.shape)
    
    if microphysics_scheme[0:3]=="ICE":
        concentrations['ii'] = mnhFile.variables['CIT'][0,:,:,:]
        concentrations['ii'][concentrations['ii']==999.] = float('nan')
    if microphysics_scheme[0:3] =="LIM" :
        concentrations['rr'] = mnhFile.variables['CRAIN'][0,:,:,:] #former name: CRAINT
        concentrations['rr'][concentrations['rr']==999.]=float('nan')
        concentrations['ii'] = mnhFile.variables['CICE'][0,:,:,:] #former name: CICET
        concentrations['ii'][concentrations['ii']==999.]=float('nan')
    concentrations['rr']*=rho
    concentrations['ii']*=rho
    return concentrations
