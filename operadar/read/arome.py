#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""
import numpy as np
import epygram
import pickle as pkl
import time as tm

from read.util_with_epygram import *
from operadar_utils import hydrometeorModel_from_hydrometeorDict, get_lat_lon_from_subdomain


#================= Read Arome variables =======================================
"""
Read Arome 3D variables in ncfile: pressure, temperature, hydrometeor contents
and concentrations 
INPUT : 
    - micro
    - extract_once
    - hydrometeors : list of hydrometeors to read in AROME file
    - moments : dict with number of moments for each hydrometeor type
    - filePath
    - lon_min,lon_max,lat_min,lat_max
OUTPUT: 
    - M (kg/m3): dict of hydrometeor contents 
    - Nc dict (/): dict of concentrations
    - Tc (°C), CC, CCI, lon, lat, Z
"""

def read_arome(filePath: str,
               micro: str,
               extract_once: bool,
               hydrometeors: dict,
               #moments: dict,
               #CCIconst: float,
               subDomain:list[float]|None,
               ):   
    
    epygram.init_env() # mandatory
    
    hydromet_list = hydrometeorModel_from_hydrometeorDict(hydrometeors)
    
    print("\tAROME .fa file: ",filePath)
    loaded_epygram_file = epygram.formats.resource(filePath, openmode = 'r', fmt = 'FA')
    ps = loaded_epygram_file.readfield('SURFPRESSION')
    X_res=loaded_epygram_file.geometry.grid['X_resolution']
    Y_res=loaded_epygram_file.geometry.grid['Y_resolution']
    
    
    if extract_once :
        # Vertical levels values
        A = [level[1]['Ai'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        B = [level[1]['Bi'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        # Save in temporary pickle files
        tmpA = open('tmpA.obj', 'wb') ; pkl.dump(A,tmpA) ; tmpA.close()
        tmpB = open('tmpB.obj', 'wb') ; pkl.dump(B,tmpB) ; tmpB.close()
    else :
        # Read previously saved pickle files
        tmpA = open('tmpA.obj', 'rb') ; A = pkl.load(tmpA) ; tmpA.close()
        tmpB = open('tmpB.obj', 'rb') ; B = pkl.load(tmpB) ; tmpB.close()
    

    if subDomain != None :
        lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
        imin,jmin=(np.round(ps.geometry.ll2ij(lon_min,lat_min)).astype(int))
        imax,jmax=(np.round(ps.geometry.ll2ij(lon_max,lat_max)).astype(int))
        arome_file = epygram.resources.SubdomainResource(
                     resource=loaded_epygram_file, openmode='r', name='Subdomain',
                     subarray=dict(imin=imin, imax=imax, jmin=jmin, jmax=jmax),
                     )
    else:
        arome_file=loaded_epygram_file
    
    
    if extract_once : 
        [lon, lat] = get_2D_lat_lon_epygram(arome_file)
    [p, psurf, pdep, phis] = get_geometry(arome_file, A, B)
    
    # ======== Hydrometeor contents, concentrations and temperature
    [M, T, R]  = get_contents_and_T(arome_file, p, hydromet_list)
    Tc=T-273.15
    Nc = get_concentrations(arome_file, p, hydromet_list,moments,CCIconst)
    
    # ========= Horizontal grid
    Y = Y_res*np.arange(Tc.shape[1]).astype('i4')
    X = X_res*np.arange(Tc.shape[2]).astype('i4')

        
    # ========= Altitude z of each level
    Z = get_altitude(A, B, T, p, pdep, psurf, phis, R)

    # ========= Close file
    loaded_epygram_file.close()
    arome_file.close()
    
    if extract_once : return  X, Y, Z, lon, lat, M, Nc, Tc
    else : return X, Y, Z, None, None, M, Nc, Tc

# =============================================================================




 #================ Appel module directement ====================================    
if __name__ == '__main__':

    # ---- to be included in conf file
    micro="ICE3"
    #rep="/home/augros/DONNEES/AROME/20220816/PEAROME/R09/"
    rep="/cnrm/precip/users/augros/DONNEES/AROME/"
    file="historic.arome.franmg-01km30+0008:00.fa"
    filePath=rep+file
    
    # === Infos === #
    #loaded_epygram_file.listfields()
    #loaded_epygram_file.what()
     # -------------------------------------

 # =============================================================================
