#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 4 2023

@author: augrosc

Configuration file for operad.py
"""

# ==========  Model simulation options ===============
model='Arome'
micro = "ICE3" # CLOE ICE3 / LIMA_SG / LIMA_AG / ICE4
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
list_types=['vv','cc','rr','ii','ss','gg']
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
lat_min,lat_max=47.0,50.
lon_min,lon_max=1.,4.


# ========== Directories / files name options =========
timelist=[8] #range(1,36) #[5,6,7,8] #ech=[20] #[36]

pathmodel="/home/augros/DONNEES/AROME/20220816/PEAROME/R09/"
#pathmodel="/cnrm/precip/users/augros/DONNEES/AROME/"
filestart="historic.arome.franmg-01km30+00" #08:00.fa"

# Tmatrix directory
table_ind="" # number of the selected Tmatrix table 
#repTmat="/cnrm/precip/users/augros/DONNEES/TMATRIX" #"." 
repTmat="/home/augros/Programmes/OPERADAR/"

# Output files
pathfick=pathmodel+'k'+MixedPhase+'/'
pathTmat=repTmat+"/OUTPUT/"


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius

