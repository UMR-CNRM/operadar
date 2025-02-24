#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""
import epygram
import time as tm
import numpy as np
from pathlib import Path

from operadar.read.with_epygram import *
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.utils.formats_data import get_lat_lon_from_subdomain


def read_arome(filePath:Path,
               hydrometeorMoments:dict[int],
               subDomain:list[float]|None,
               verbose:bool,
               )-> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,dict[np.ndarray],dict[np.ndarray],np.ndarray]:   
    """Read and extract data from an AROME.fa file
    
    Args:
        filePath (str): input file path.
        hydrometeorMoments (dict): dictionary of form {'hydrometeor key':number of moments}.
        subDomain (list[float] | None): either a list of 4 float or None.
        verbose (bool) : show more messages to the user.

    Returns:
        X (ndarray): 1D horizontal coordinates in m
        Y (ndarray): 1D horizontal coordinates in m
        Z (ndarray): 1D array of vertical coordinates in model pressure levels
        lon (ndarray): 2D array of longitude coordinates
        lat (ndarray): 2D array of latitude coordinates
        M (dict[ndarray]): dictionary of 3D contents for each hydrometeor 
        Nc (dict[ndarray]): dictionary of 3D number concentrations for each hydrometeor
        Tc (ndarray) : 3D temperature in Celsius
    """
    
    epygram.init_env() # mandatory
    
    hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                          datatype='model',
                                                          )
    print("\tAROME .fa file:",filePath)
    if verbose : deb=tm.time()
    loaded_epygram_file = epygram.formats.resource(filename=filePath,
                                                   openmode = 'r',
                                                   fmt = 'FA',
                                                   )
    if verbose : print('\t\tLoaded file in ',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    X_res=loaded_epygram_file.geometry.grid['X_resolution']
    Y_res=loaded_epygram_file.geometry.grid['Y_resolution']
    
    if verbose : print('\t\tRed surface pressure, x, y in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    
    # Hybrid pressure coefficients
    A = [level[1]['Ai'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]
    B = [level[1]['Bi'] for level in loaded_epygram_file.geometry.vcoordinate.grid['gridlevels']][1:]     
    if verbose : print('\t\tExtracted A and B hybrid pressure coefficients in',round(tm.time()-deb,6),'seconds'); deb=tm.time()
    
    if subDomain != None :
        ps = loaded_epygram_file.readfield('SURFPRESSION')
        (lon_min, lon_max, lat_min, lat_max) = get_lat_lon_from_subdomain(subDomain) 
        imin,jmin=(np.round(ps.geometry.ll2ij(lon_min,lat_min)).astype(int))
        imax,jmax=(np.round(ps.geometry.ll2ij(lon_max,lat_max)).astype(int))
        arome_file = epygram.resources.SubdomainResource(resource=loaded_epygram_file, openmode='r',
                                                         name='Subdomain',
                                                         subarray=dict(imin=imin, imax=imax, jmin=jmin, jmax=jmax),
                                                         )
        if verbose : print('\t\tSubdomain extracted in',round(tm.time()-deb,6),'seconds'); deb=tm.time()
    else:
        arome_file=loaded_epygram_file
        if verbose : print('\t\tNo subdomain provided, will use all the points.'); deb=tm.time()
    
    # Extract 2D lon and lat fields only once if multiple iterations over the same (sub)domain
    [lon, lat] = get_2D_lat_lon_epygram(epygram_file=arome_file)
    if verbose : print('\t\tGot 2D lat/lon fields in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    [p, psurf, pdep, geosurf] = get_geometry(epygram_file=arome_file,
                                             hybrid_pressure_coefA=A,
                                             hybrid_pressure_coefB=B,
                                             )
    if verbose : print('\t\tGot 3D pressure departure, surface pressure, and geopotential in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    [M, T, R]  = get_contents_T_and_R(epygram_file=arome_file,
                                      pressure=p,
                                      hydrometeors=hydromet_list,
                                      )
    del p
    if verbose : print('\t\tGot 3D contents (kg/m3) and temperature in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    Nc = get_concentrations(epygram_file=arome_file,
                            hydrometeorsConfig=hydrometeorMoments,
                            content=M,
                            temperature=T,
                            )
    if verbose : print('\t\tGot 3D concentrations in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    T=T-273.15
    X = X_res*np.arange(T.shape[2]).astype('i4')
    Y = Y_res*np.arange(T.shape[1]).astype('i4')
    Alt = get_altitude(hybrid_pressure_coefA=A,
                       hybrid_pressure_coefB=B,
                       temperature=T,
                       pressure_departure=pdep,
                       surface_pressure=psurf,
                       surface_geopotential=geosurf,
                       specific_gas_constant=R,
                       )
    if verbose : print('\t\tComputed altitude 3D field in',round(tm.time()-deb,3),'seconds'); deb=tm.time()
    
    loaded_epygram_file.close()
    arome_file.close()
    return  X, Y, Alt, lon, lat, M, Nc, T
