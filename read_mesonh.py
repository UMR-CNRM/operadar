#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""

import numpy as np
import operad_conf as cf
from netCDF4 import Dataset

#============== Read MesoNH variables ===============
"""
Read MesoNH 3D variables in ncfile: pressure, temperature, hydrometeor contents 
"""
def read_mesonh(micro,modelfile):
    
    # === Model file
    #time='00'+str(time)          
    #modelfile=cf.pathmodel+cf.filestart+time[-3:]+'.nc'
    print("Reading "+modelfile)
    
    # === Extract Dataset 
    print("Reading ncfile: ",modelfile)
    ncfile1 = Dataset(modelfile,'r')
    #print(ncfile1.variables.keys())

    
    # === Geometry: X, Y, Z coordinates and radar cover mask
    X=ncfile1.variables['XHAT'][:]
    Y=ncfile1.variables['YHAT'][:]
    Z=ncfile1.variables['ZHAT'][:]

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
    list_t_full=['vv','cc','rr','ii','ss','gg']
    list_hydro=['RVT','RCT','RRT','RIT','RST','RGT','RHT']
    name_hydro={}
    M={}
    for t in cf.list_types:
        M[t] = np.empty(Tc.shape)    

    # Arrays initialisation
    for it,t in enumerate(list_t_full):
        name_hydro[t]=list_hydro[it]
    

    for t in cf.list_types:
        M[t]=ncfile1.variables[name_hydro[t]][0,:,:,:]*rho3D[:,:,:]
        M[t][M[t]==999.]=float('nan')

    if(cf.micro =="ICE3" or cf.micro =="ICE4"):
        CCI=ncfile1.variables['CIT'][0,:,:,:]
        CCI[CCI==999.]=float('nan')
        CC=np.empty(Tc.shape)
    if(cf.micro =="LIMA" or cf.micro =="LIMT" or cf.micro =="LIMA_SG" or cf.micro =="LIMA_AG"):
        CC=ncfile1.variables['CRAINT'][0,:,:,:]
        CC[CC==999.]=float('nan')
        CCI=ncfile1.variables['CICET'][0,:,:,:]
        CCI[CCI==999.]=float('nan')
    CC*=rho3D
    CCI*=rho3D
    
    # =====================================================
    
    print("End reading model variables")
    
    return M, Tc, CC, CCI, X, Y, Z
#=====================================================================