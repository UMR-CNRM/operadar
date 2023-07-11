#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""

import numpy as np
import epygram
import bronx
import read_arome_lib as arolib



#================= Read Arome variables =======================================
"""
Read Arome 3D variables in ncfile: pressure, temperature, hydrometeor contents
and concentrations 
INPUT : microphysics,modelfile,lon_min,lon_max,lat_min,lat_max
OUTPUT: M, Tc, CC, CCI, lon, lat, Z
"""

def read_arome(modelfile: str,
               microphysics: str,
               hydrometeors_list: list,
               lonmin: float = -5.2,
               lonmax: float = 8.3,
               latmin: float = 41.3,
               latmax: float = 51.15,
               ):   
    
    epygram.init_env()
    
    print("Reading AROME fa file: ",modelfile)
    ficA = epygram.formats.resource(modelfile, openmode = 'r', fmt = 'FA')
    ps = ficA.readfield('SURFPRESSION')
    
    # === Infos === #
    #ficA.listfields()
    #ficA.what()
    
    # === Zoom === #
    imin,jmin=(np.round(ps.geometry.ll2ij(lonmin,latmin)).astype(int))
    imax,jmax=(np.round(ps.geometry.ll2ij(lonmax,latmax)).astype(int))
    
    #Fichier FA avec SubdomainResource
    ficsubdo = epygram.resources.SubdomainResource(resource=ficA, openmode='r', name='Subdomain',
                                                  subarray=dict(imin=imin, imax=imax, jmin=jmin, jmax=jmax))
    #ficsubdo.readfield('S089RAIN').cartoplot()[0].savefig('subdo.png')
    
    
    # ======== Horizontal, vertical coordinates, pressure
    [p, psurf, pdep, phis, A, B, lon, lat] = arolib.get_geometry(ficA,ficsubdo)
      
    # ======== Hydrometeor contents and temperature
    [M, T, R]  = arolib.get_contents_and_T(ficsubdo, p, hydrometeors_list)
    Tc=T-273.15
    [CC , CCI] = arolib.get_concentrations(ficsubdo, p,microphysics,hydrometeors_list)
    
        
    # ========= Altitude z of each level
    Z=arolib.get_altitude(A, B, T, p, pdep, psurf, phis, R)
    
     
    # ========= Close file
    ficsubdo.close()
        
    return M, Tc, CC, CCI, lon, lat, Z

# =============================================================================




 #================ Appel module directement ====================================    
if __name__ == '__main__':

    # ---- to be included in conf file
    microphysics="ICE3"
    #rep="/home/augros/DONNEES/AROME/20220816/PEAROME/R09/"
    rep="/cnrm/precip/users/augros/DONNEES/AROME/"
    file="historic.arome.franmg-01km30+0008:00.fa"
    modelfile=rep+file
     # -------------------------------------

 # =============================================================================
