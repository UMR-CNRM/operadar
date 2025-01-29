#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""

import numpy as np
"import operad_conf as cf" # NOT USED ANYMORE --> make call directly into functions
from netCDF4 import Dataset
from operad_utils import hydrometeorModel_from_hydrometeorDict


#============== Read MesoNH variables ===============
"""
Read MesoNH 3D variables in ncfile: pressure, temperature, hydrometeor contents 
"""
def read_mesonh(modelfile: str,
                micro: str,
                lon_min: float,
                lon_max: float,
                lat_min: float,
                lat_max: float,
                hydrometeors: dict,
                real_case: bool(),
               ):
    
    # === Model file
    print("Reading "+modelfile)
    
    # === Extract Dataset 
    print("Reading ncfile: ",modelfile)
    ncfile1 = Dataset(modelfile,'r')
    #print(ncfile1.variables.keys())
    
    
    # === Geometry: X, Y, Z coordinates and radar cover mask
    X=ncfile1.variables['XHAT'][:]
    Y=ncfile1.variables['YHAT'][:]
    ZHAT=ncfile1.variables['ZHAT'][:] # height above ground (m)
    ZS=ncfile1.variables['ZS'][:] # altitude (m) of model surface (ground)
    Z=np.empty((ZHAT.shape[0],ZS.shape[0],ZS.shape[1]))
    for level in range(ZHAT.shape[0]):
        Z[level,:,:]=ZS[:,:]+ZHAT[level]
    time=ncfile1.variables['time'][:]
    print("MesoNH ncfile time = ",time)
    
    # === Get lat lon if real case
    if (real_case):
        LAT = ncfile1.variables['latitude'][:]
        LON = ncfile1.variables['longitude'][:]
        mask_zoom=((LON>lon_min) & (LON<lon_max) & (LAT>lat_min) & (LAT<lat_max))
        [ilon,jlat]=np.where(mask_zoom)
        i_lonmin,i_lonmax=np.nanmin(ilon),np.nanmax(ilon)
        j_latmin,j_latmax=np.nanmin(jlat),np.nanmax(jlat)
    else:
        LAT=float('nan')
        LON = float('nan')
    
    
    #    #This is for taking care of cropped netcdf which still contain XHAT and YHAT with non cropped indices
    #    if X.shape!=LON.shape : 
    #        ilatmin,ilatmax,ilonmin,ilonmax = ope_lib.crop_latlon(modelfile,LAT[0],LAT[-1],LON[0],LON[-1])
    #        X=ncfile1.variables['XHAT'][ilonmin:ilonmax+1]
    #        Y=ncfile1.variables['YHAT'][ilatmin:ilatmax+1]
    
    
    # =======================
    
    # === Pressure and temperature and dry air density
    p=ncfile1.variables['PABST'][0,:,:,:]
    p[np.where(p==999.)]=float('nan')
    Th=ncfile1.variables['THT'][0,:,:,:]
    Th[np.where(Th==999.)]=float('nan')
    Tc=Th*((p/100000)**(0.4/1.4))-273.15 
    del Th, p
    
    rhodref = ncfile1.variables['RHOREFZ'][:]
    rho3D=np.ones(Tc.shape)
    
    IKE=Tc.shape[0]
    for k in range(IKE):    
        rho3D[k,:,:]=rhodref[k]
    # =====================
    
    # === Hydrometeors contents and concentrations
    
    hydromet_list = hydrometeorModel_from_hydrometeorDict(hydrometeors)
    
    list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    list_hydro=['RVT','RCT','RRT','RIT','RST','RGT','RHT']
    name_hydro={}
    M={}
    for t in hydromet_list:
        M[t] = np.empty(Tc.shape)    
    
    # Arrays initialisation
    for it,t in enumerate(list_t_full):
        name_hydro[t]=list_hydro[it]
    
    
    for t in hydromet_list:
        M[t]=ncfile1.variables[name_hydro[t]][0,:,:,:]*rho3D[:,:,:]
        M[t][M[t]==999.]=float('nan')
    
    if(micro =="ICE3" or micro =="ICE4"):
        CCI=ncfile1.variables['CIT'][0,:,:,:]
        CCI[CCI==999.]=float('nan')
        CC=np.empty(Tc.shape)
    if(micro =="LIMA" or micro =="LIMT" or micro =="LIMASG" or micro =="LIMAAG"):
        CC=ncfile1.variables['CRAIN'][0,:,:,:] #former name: CRAINT
        CC[CC==999.]=float('nan')
        CCI=ncfile1.variables['CICE'][0,:,:,:] #former name: CICET
        CCI[CCI==999.]=float('nan')
    CC*=rho3D
    CCI*=rho3D


    # ===== Calcul of the grid considering orography
    #Orography =ncfile1.variables['ZS'][:]
    #if np.any(Orography > 0):
    #    Z = ope_lib.compute_grid_alt(var3D,ztop,orography,level)
    #else:
    #    Z=ncfile1.variables['level'][:] # but in 3D shape ?
    
    # =====================================================
    print("Z.shape",Z.shape)
    print("End reading model variables")
    
    if (real_case):
        Mzoom={}
        for t in hydromet_list:
            Mzoom[t]=M[t][:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        
        Tczoom=Tc[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        CCzoom=Tc[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        CCIzoom=CCI[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        LATzoom=LAT[i_lonmin:i_lonmax,j_latmin:j_latmax]
        LONzoom=LON[i_lonmin:i_lonmax,j_latmin:j_latmax]
        Xzoom=X[j_latmin:j_latmax]
        Yzoom=Y[i_lonmin:i_lonmax]
        Zzoom=Z[:,i_lonmin:i_lonmax,j_latmin:j_latmax]
        
        return Mzoom, Tczoom, CCzoom, CCIzoom, LATzoom, LONzoom, Xzoom, Yzoom, Zzoom
    
    else:
        return M, Tc, CC, CCI, LAT, LON, X, Y, Z
#=====================================================================
