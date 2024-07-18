#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""
import numpy as np
import epygram
import bronx
import pickle as pkl
import read_arome_lib as arolib
import time as tm


#================= Read Arome variables =======================================
"""
Read Arome 3D variables in ncfile: pressure, temperature, hydrometeor contents
and concentrations 
INPUT : micro,modelfile,lon_min,lon_max,lat_min,lat_max
OUTPUT: M, Tc, CC, CCI, lon, lat, Z
"""

def read_arome(modelfile: str,
               micro: str,
               extract_once: bool,
               hydrometeors_list: list,
               lon_min: float = -5.2,
               lon_max: float = 8.3,
               lat_min: float = 41.3,
               lat_max: float = 51.15,
               ):   
    
    epygram.init_env()
    
    print("  AROME fa file: ",modelfile)
    ficA = epygram.formats.resource(modelfile, openmode = 'r', fmt = 'FA')
    ps = ficA.readfield('SURFPRESSION')
    X_res=ficA.geometry.grid['X_resolution']
    Y_res=ficA.geometry.grid['Y_resolution']
    
    if extract_once :
        # Vertical levels values
        A = [level[1]['Ai'] for level in ficA.geometry.vcoordinate.grid['gridlevels']][1:]
        B = [level[1]['Bi'] for level in ficA.geometry.vcoordinate.grid['gridlevels']][1:]
        # Save in temporary pickle files
        tmpA = open('tmpA.obj', 'wb') ; pkl.dump(A,tmpA) ; tmpA.close()
        tmpB = open('tmpB.obj', 'wb') ; pkl.dump(B,tmpB) ; tmpB.close()
    else :
        # Read previously saved pickle files
        tmpA = open('tmpA.obj', 'rb') ; A = pkl.load(tmpA) ; tmpA.close()
        tmpB = open('tmpB.obj', 'rb') ; B = pkl.load(tmpB) ; tmpB.close()
    
    # === Infos === #
    #ficA.listfields()
    #ficA.what()
    
    # === Extrat subdomain === #
    imin,jmin=(np.round(ps.geometry.ll2ij(lon_min,lat_min)).astype(int))
    imax,jmax=(np.round(ps.geometry.ll2ij(lon_max,lat_max)).astype(int))
    ficsubdo = epygram.resources.SubdomainResource(resource=ficA, openmode='r', name='Subdomain',
                                                  subarray=dict(imin=imin, imax=imax, jmin=jmin, jmax=jmax))
    #ficsubdo.readfield('S089RAIN').cartoplot()[0].savefig('subdo.png')
    
    # ======== Horizontal, vertical coordinates, pressure
    if extract_once : 
        [lon, lat] = arolib.get_lat_lon_epygram(ficsubdo)
    [p, psurf, pdep, phis] = arolib.get_geometry(ficsubdo, A, B)
    
    # ======== Hydrometeor contents and temperature
    [M, T, R]  = arolib.get_contents_and_T(ficsubdo, p, hydrometeors_list)
    Tc=T-273.15
    [CC , CCI] = arolib.get_concentrations(ficsubdo, p,micro,hydrometeors_list)
    
    # ========= Horizontal grid
    Y = Y_res*np.arange(Tc.shape[1]).astype('i4')
    X = X_res*np.arange(Tc.shape[2]).astype('i4')

        
    # ========= Altitude z of each level
    Z = arolib.get_altitude(A, B, T, p, pdep, psurf, phis, R)

    # ========= Close file
    ficA.close()
    ficsubdo.close()
    
    if extract_once : return M, Tc, CC, CCI, Z, X, Y, lon, lat
    else : return M, Tc, CC, CCI, Z, X, Y

# =============================================================================




 #================ Appel module directement ====================================    
if __name__ == '__main__':

    # ---- to be included in conf file
    micro="ICE3"
    #rep="/home/augros/DONNEES/AROME/20220816/PEAROME/R09/"
    rep="/cnrm/precip/users/augros/DONNEES/AROME/"
    file="historic.arome.franmg-01km30+0008:00.fa"
    modelfile=rep+file
     # -------------------------------------

 # =============================================================================
