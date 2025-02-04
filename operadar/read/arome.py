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
from operadar_utils import (
    hydrometeorModel_from_hydrometeorDict,
    get_lat_lon_from_subdomain
)



def read_arome(filePath: str, micro: str, extract_once: bool, hydrometeors: dict, subDomain:list[float]|None):   
    """Read and extract data from an AROME.fa file
    
    Args:
        filePath (str): _description_
        micro (str): _description_
        extract_once (bool): _description_
        hydrometeors (dict): _description_
        subDomain (list[float] | None): _description_

    Returns:
        X (np.ndarray): 1D horizontal coordinates in m
        Y (np.ndarray): 1D horizontal coordinates in m
        Z (np.ndarray): 1D array of vertical coordinates in model pressure levels
        lon (np.ndarray | None): 2D array of longitude coordinates or None if extract_once = False
        lat (np.ndarray | None): 2D array of latitude coordinates or None if extract_once = False
        M (dict[np.ndarray]): dictionnary of 3D contents for each hydrometeor 
        Nc (dict[np.ndarray]): dictionnary of 3D number concentrations for each hydrometeor
        Tc (np.ndarray) : 3D temperature in Celsius
        
    """
    
    epygram.init_env() # mandatory
    
    hydromet_list = hydrometeorModel_from_hydrometeorDict(hydrometeors)
    
    print("\tAROME .fa file: ",filePath)
    loaded_epygram_file = epygram.formats.resource(filename=filePath,
                                                   openmode = 'r',
                                                   fmt = 'FA',
                                                   )
    ps = loaded_epygram_file.readfield('SURFPRESSION')
    X_res=loaded_epygram_file.geometry.grid['X_resolution']
    Y_res=loaded_epygram_file.geometry.grid['Y_resolution']
    
    if extract_once :
        # Hybrid pressure coefficients
        A = [level[1]['Ai'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        B = [level[1]['Bi'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        # Save in temporary pickle files
        tmpA = open('tmp/tmpA.obj', 'wb') ; pkl.dump(A,tmpA) ; tmpA.close()
        tmpB = open('tmp/tmpB.obj', 'wb') ; pkl.dump(B,tmpB) ; tmpB.close()
    else :
        # Read previously saved pickle files
        tmpA = open('tmp/tmpA.obj', 'rb') ; A = pkl.load(tmpA) ; tmpA.close()
        tmpB = open('tmp/tmpB.obj', 'rb') ; B = pkl.load(tmpB) ; tmpB.close()
    
    if subDomain != None :
        lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
        imin,jmin=(np.round(ps.geometry.ll2ij(lon_min,lat_min)).astype(int))
        imax,jmax=(np.round(ps.geometry.ll2ij(lon_max,lat_max)).astype(int))
        arome_file = epygram.resources.SubdomainResource(resource=loaded_epygram_file,
                                                         openmode='r',
                                                         name='Subdomain',
                                                         subarray=dict(imin=imin,
                                                                       imax=imax,
                                                                       jmin=jmin,
                                                                       jmax=jmax) )
    else:
        arome_file=loaded_epygram_file
    loaded_epygram_file.close()
    
    # Extract 2D lon and lat fields only once if multiple iterations over the same (sub)domain
    if extract_once : 
        [lon, lat] = get_2D_lat_lon_epygram(epygram_file=arome_file)
    [p, psurf, pdep, geosurf] = get_geometry(epygram_file=arome_file,
                                             A=A, B=B )
    
    [M, T, R]  = get_contents_T_and_R(epygram_file=arome_file,
                                      pressure=p,
                                      hydrometeors=hydromet_list )
    Nc = get_concentrations(epygram_file=arome_file,
                            hydrometeorsConfig=hydrometeors,
                            content=M,
                            temperature=T )
    Tc=T-273.15
    
    X = X_res*np.arange(Tc.shape[2]).astype('i4')
    Y = Y_res*np.arange(Tc.shape[1]).astype('i4')
    Z = get_altitude(A, B, T, p, pdep, psurf, geosurf, R)

    arome_file.close()
    if extract_once : return  X, Y, Z, lon, lat, M, Nc, Tc
    else : return X, Y, Z, None, None, M, Nc, Tc



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
