#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 4 2023

@author: augrosc

Configuration file for operad.py
"""
import pandas as pd
import datetime as dt
import common_settings as settings


# ==========  Model simulation options ===============
model='Arome'
micro = "ICE3"
run  = settings.run
deb = settings.deb
fin = settings.fin
step = settings.step

save_npz    = settings.save_npz
save_netcdf = settings.save_netcdf

LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['vv','cc','rr','ii','ss','gg'] # available model variables
list_types_tot = ['rr','ii','ss','gg','wg']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type
#savetxt=False # True: output file in ascii format
#              # False: output file in native python compressed format (.npz) 

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)


# ==========  Radar options =============== 
band="S"
distmax_rad = 255.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 12000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)
#latrad=
#lonrad=

# ========= Zoom ==========================
#lat_min,lat_max=42,45
#lon_min,lon_max=1,5
lat_min = settings.lat_min ; lat_max = settings.lat_max
lon_min = settings.lon_min ; lon_max = settings.lon_max

# ========== Directories / files name options =========
# Time list
datetimelist=[]
ech = deb
while ech <= fin :
    datetimelist += [ech]
    ech += step

# Vortex experiment name
if settings.run == '00' :
    expeOLIVE = 'GN49'
elif settings.run == '12' :
    expeOLIVE = 'GOIM'

# Model files
commonPath_fa  = f"/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/"
commonFilename = "historic.arome.franmg-01km30+00" #08:00.fa"
pathmodel = commonPath_fa + f"{expeOLIVE}/{deb.strftime('%Y%m%dT')}{run}00P/forecast/"

# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/cnrm/precip/users/augros/Programmes/TMATRIX/DPOLSIMUL"
# Output files
pathfick = f"/cnrm/precip/SAVE/davidcl/THESE/operadar_files/{deb.strftime('%Y%m%d')}/{run}Z_{micro}_k{MixedPhase}/{settings.radar_ids}"
pathTmat = repTmat+"/OUTPUT/"


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius