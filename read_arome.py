#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros
"""

import numpy as np

import epygram
import bronx


import operad_conf as cf
import read_arome_lib as arolib

epygram.init_env()


## ---- to be included in conf file
#micro="ICE3"
##rep="/home/augros/DONNEES/AROME/20220816/PEAROME/R09/"
#rep="/cnrm/precip/users/augros/DONNEES/AROME/"
#file="historic.arome.franmg-01km30+0008:00.fa"
##rep="/home/augros/DONNEES/AROME/20220816/3DVARFR/T00/" 
##file="historic.arome.franmg-01km30+0017:00.fa"
#lat_min,lat_max=47.0,50.
#lon_min,lon_max=1.,4.
#modelfile=rep+file
# -------------------------------------


#================= Read Arome variables =======================================
"""
Read Arome 3D variables in ncfile: pressure, temperature, hydrometeor contents 
input : micro,modelfile,lon_min,lon_max,lat_min,lat_max
output: M, Tc, CC, CCI, lon, lat, Z


"""
def read_arome(micro,timestr):   
    
    epygram.init_env()
    # ======== Open file 
    modelfile=cf.pathmodel+cf.filestart+timestr+".fa"
    
    
    print("Reading AROME fa file: ",modelfile)
    ficA = epygram.formats.resource(modelfile, openmode = 'r', fmt = 'FA')
    ps = ficA.readfield('SURFPRESSION')
    
    
    # === Infos
    #ficA.listfields()
    #ficA.what()
    
    # ======= Zoom
    imin,jmin=(np.round(ps.geometry.ll2ij(cf.lon_min,cf.lat_min)).astype(int))
    imax,jmax=(np.round(ps.geometry.ll2ij(cf.lon_max,cf.lat_max)).astype(int))
    
    #Fichier FA avec SubdomainResource
    ficsubdo = epygram.resources.SubdomainResource(resource=ficA, openmode='r', name='Subdomain',
                                                  subarray=dict(imin=imin, imax=imax, jmin=jmin, jmax=jmax))
    #ficsubdo.readfield('S089RAIN').cartoplot()[0].savefig('subdo.png')
    
    
    # ======== Horizontal, vertical coordinates, pressure
    [p, psurf, pdep, phis, A, B, lon, lat]=arolib.get_geometry(ficA,ficsubdo)
      
    # ======== Hydrometeor contents and temperature
    [M, T, R]=arolib.get_contents_t(ficsubdo, p)
    
    Tc=T-273.15
    CC=np.empty(Tc.shape)
    CCI=cf.CCIconst*np.ones(Tc.shape)
        
    # ========= Altitude z of each level
    Z=arolib.get_altitude(A, B, T, p, pdep, psurf, phis, R)
    
     
    # ========= Close file
    ficsubdo.close()
        
    return M, Tc, CC, CCI, lon, lat, Z

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
