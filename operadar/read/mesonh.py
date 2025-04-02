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

    X, Y, Z = get_geometry(mnh_file=mnh_file)
    if verbose : print('\t\tGot geometry (X,Y,Z) in',round(tm.time()-deb,6),'seconds');deb=tm.time()
    
    if real_case:
        LAT = mnh_file.variables['latitude'][:]
        LON = mnh_file.variables['longitude'][:]
        if subDomain != None :
            lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
            mask_zoom=((LON>lon_min) & (LON<lon_max) & (LAT>lat_min) & (LAT<lat_max))
            [ilon,jlat]=np.where(mask_zoom) ; print(ilon)
            i_lonmin,i_lonmax=np.nanmin(ilon),np.nanmax(ilon)
            j_latmin,j_latmax=np.nanmin(jlat),np.nanmax(jlat)
            if verbose : print('\t\tSubdomain extracted in',round(tm.time()-deb,6),'seconds'); deb=tm.time()
        else :
            if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
    else:
        LAT=float('nan')
        LON = float('nan')
        if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
    
    #    #This is for taking care of cropped netcdf which still contain XHAT and YHAT with non cropped indices
    #    if X.shape!=LON.shape : 
    #        ilatmin,ilatmax,ilonmin,ilonmax = ope_lib.crop_latlon(filePath,LAT[0],LAT[-1],LON[0],LON[-1])
    #        X=mnh_file.variables['XHAT'][ilonmin:ilonmax+1]
    #        Y=mnh_file.variables['YHAT'][ilatmin:ilatmax+1]
    
    # =======================
    
    # === Pressure and temperature and dry air density
    p=mnh_file.variables['PABST'][0,:,:,:]
    p[np.where(p==999.)]=float('nan')
    Th=mnh_file.variables['THT'][0,:,:,:]
    Th[np.where(Th==999.)]=float('nan')
    Tc=Th*((p/100000)**(0.4/1.4))-273.15 
    del Th, p
    
    rhodref = mnh_file.variables['RHOREFZ'][:]
    rho3D=np.ones(Tc.shape)
    
    IKE=Tc.shape[0]
    for k in range(IKE):    
        rho3D[k,:,:]=rhodref[k]
    
    # === Hydrometeors contents and concentrations
    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments, datatype='model')
    name_hydro = link_varname_with_mesonh_name()
    #list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    #list_hydro=['RVT','RCT','RRT','RIT','RST','RGT','RHT']
    #name_hydro={}
    
    M = {}
    for key in hydromet_list:
        M[key] = np.empty(Tc.shape)
        M[key] = mnh_file.variables[name_hydro[key]][0,:,:,:]*rho3D[:,:,:] # kg/kg of dry air
        M[key][M[key]==999.] = float('nan')
    
    Nc = {}
    for key in hydromet_list:
        Nc[key] = np.zeros(Tc.shape)
    
    if micro[0:3]=="ICE":
        Nc['ii'] = mnh_file.variables['CIT'][0,:,:,:]
        Nc['ii'][Nc['ii']==999.] = float('nan')
    if micro[0:3] =="LIM" :
        Nc['rr'] = mnh_file.variables['CRAIN'][0,:,:,:] #former name: CRAINT
        Nc['rr'][Nc['rr']==999.]=float('nan')
        Nc['ii'] = mnh_file.variables['CICE'][0,:,:,:] #former name: CICET
        Nc['ii'][Nc['ii']==999.]=float('nan')
    Nc['rr']*=rho3D
    Nc['ii']*=rho3D

    # ===== Calcul of the grid considering orography
    #Orography =mnh_file.variables['ZS'][:]
    #if np.any(Orography > 0):
    #    Z = ope_lib.compute_grid_alt(var3D,ztop,orography,level)
    #else:
    #    Z=mnh_file.variables['level'][:] # but in 3D shape ?
    
    print("End reading model variables")
    
    if real_case and subDomain != None :
        Mzoom,Nczoom={},{}
        for key in hydromet_list:
            Mzoom[key]=M[key][:,i_lonmin:i_lonmax,j_latmin:j_latmax]
            Nczoom[key]=Nczoom[key][:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        
        Tczoom=Tc[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        LATzoom=LAT[i_lonmin:i_lonmax,j_latmin:j_latmax]
        LONzoom=LON[i_lonmin:i_lonmax,j_latmin:j_latmax]
        Xzoom=X[j_latmin:j_latmax]
        Yzoom=Y[i_lonmin:i_lonmax]
        Zzoom=Z[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        
        return Xzoom, Yzoom, Zzoom, LONzoom, LATzoom, Mzoom, Nczoom, Tczoom 
    else:
        return X, Y, Z, LON, LAT, M, Nc, Tc



def get_geometry(mnh_file:Dataset):
    X = mnh_file.variables['XHAT'][:]
    Y = mnh_file.variables['YHAT'][:]
    ZHAT = mnh_file.variables['ZHAT'][:] # height above ground (m)
    ZS = mnh_file.variables['ZS'][:] # altitude (m) of model surface (ground)
    Z = np.empty((ZHAT.shape[0],ZS.shape[0],ZS.shape[1]))
    for level in range(ZHAT.shape[0]):
        Z[level,:,:]=ZS[:,:]+ZHAT[level]
    return X, Y, Z