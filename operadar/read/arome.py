#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""
import epygram
import time as tm
import numpy as np
import pickle as pkl
from pathlib import Path
from pandas import Timestamp

from operadar.read.with_epygram import *
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.utils.formats_data import get_lat_lon_from_subdomain, check_correspondance_datetime_and_file


def read_arome(filePath:str,
               date_time:Timestamp,
               extract_once: bool,
               hydrometeorMoments:dict[int],
               subDomain:list[float]|None,
               testing=False,
               )-> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,dict[np.ndarray],dict[np.ndarray],np.ndarray]:   
    """Read and extract data from an AROME.fa file
    
    Args:
        filePath (str): _description_
        extract_once (bool): _description_
        hydrometeorMoments (dict): _description_
        subDomain (list[float] | None): _description_

    Returns:
        X (ndarray): 1D horizontal coordinates in m
        Y (ndarray): 1D horizontal coordinates in m
        Z (ndarray): 1D array of vertical coordinates in model pressure levels
        lon (ndarray | None): 2D array of longitude coordinates or None if extract_once = False
        lat (ndarray | None): 2D array of latitude coordinates or None if extract_once = False
        M (dict[ndarray]): dictionnary of 3D contents for each hydrometeor 
        Nc (dict[ndarray]): dictionnary of 3D number concentrations for each hydrometeor
        Tc (ndarray) : 3D temperature in Celsius
    """
    
    epygram.init_env() # mandatory
    
    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                          datatype='model',
                                                          )
    print("\tAROME .fa file: ",filePath)
    if testing : deb=tm.time()
    loaded_epygram_file = epygram.formats.resource(filename=filePath,
                                                   openmode = 'r',
                                                   fmt = 'FA',
                                                   )
    check_correspondance_datetime_and_file(loaded_file=loaded_epygram_file,
                                           date_time_user=date_time,
                                           )
    if testing : print('load',tm.time()-deb); deb=tm.time()
    
    ps = loaded_epygram_file.readfield('SURFPRESSION')
    X_res=loaded_epygram_file.geometry.grid['X_resolution']
    Y_res=loaded_epygram_file.geometry.grid['Y_resolution']
    
    if testing : print('read ps, x, y',tm.time()-deb); deb=tm.time()
    
    if extract_once :
        # Hybrid pressure coefficients
        A = [level[1]['Ai'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        B = [level[1]['Bi'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
        # Save in temporary pickle files
        Path("./tmp/").mkdir(exist_ok=True, parents=True)
        tmpA = open('./tmp/hybrid_pressure_coefA.obj', 'wb') ; pkl.dump(A,tmpA) ; tmpA.close()
        tmpB = open('./tmp/hybrid_pressure_coefB.obj', 'wb') ; pkl.dump(B,tmpB) ; tmpB.close()
    else :
        # Read previously saved pickle files
        tmpA = open('./tmp/hybrid_pressure_coefA.obj', 'rb') ; A = pkl.load(tmpA) ; tmpA.close()
        tmpB = open('./tmp/hybrid_pressure_coefB.obj', 'rb') ; B = pkl.load(tmpB) ; tmpB.close()
        
    if testing : print('extract A B and save it',tm.time()-deb); deb=tm.time()
    
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
        if testing : print('extract subdomain',tm.time()-deb); deb=tm.time()
    else:
        arome_file=loaded_epygram_file
        if testing : print('no subdomain',tm.time()-deb); deb=tm.time()
    loaded_epygram_file.close()
    
    # Extract 2D lon and lat fields only once if multiple iterations over the same (sub)domain
    if extract_once : 
        [lon, lat] = get_2D_lat_lon_epygram(epygram_file=arome_file)
    [p, psurf, pdep, geosurf] = get_geometry(epygram_file=arome_file,
                                             hybrid_pressure_coefA=A,
                                             hybrid_pressure_coefB=B,
                                             )
    if testing : print('get lat lon and geometry pressure',tm.time()-deb); deb=tm.time()
    
    [M, T, R]  = get_contents_T_and_R(epygram_file=arome_file,
                                      pressure=p,
                                      hydrometeors=hydromet_list,
                                      )
    if testing : print('get contents T and R',tm.time()-deb); deb=tm.time()
    
    Nc = get_concentrations(epygram_file=arome_file,
                            hydrometeorsConfig=hydrometeorMoments,
                            content=M,
                            temperature=T,
                            )
    if testing : print('get concentrations',tm.time()-deb); deb=tm.time()
    
    Tc=T-273.15
    X = X_res*np.arange(Tc.shape[2]).astype('i4')
    Y = Y_res*np.arange(Tc.shape[1]).astype('i4')
    Z = get_altitude(hybrid_pressure_coefA=A,
                     hybrid_pressure_coefB=B,
                     temperature=T,
                     pressure_departure=pdep,
                     surface_pressure=psurf,
                     surface_geopotential=geosurf,
                     specific_gas_constant=R,
                     )
    if testing : print('get altitude',tm.time()-deb); deb=tm.time()
    
    arome_file.close()
    if extract_once : return  X, Y, Z, lon, lat, M, Nc, Tc
    else : return X, Y, Z, None, None, M, Nc, Tc
