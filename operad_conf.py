#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 4 2023

@author: augrosc

Configuration file for operad.py
"""
import numpy as np

# ==========  Model simulation options ===============
model='MesoNH'
micro = "ICE3" # CLOE ICE3 / LIMA_SG / LIMA_AG / ICE4
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
list_types=['vv','cc','rr','ii','ss','gg']
list_types_tot=['rr','ii','ss','gg','wg']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type
#savetxt=False # True: output file in ascii format
#              # False: output file in native python compressed format (.npz) 

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)

# ==========  Radar options =============== 
band="C"
distmax_rad = 1000.*1000 #255*1000 # Max distance for dual-pol variables computation
alt_max = 12000. # Max altitude for dual-pol variables computation
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)
latrad=float("nan")
lonrad=float("nan")
Radpos = np.array([latrad,lonrad])

# ========== Directories / files name options =========
timelist=["18","19","20","21","22","23","24","25","26","27","28","29","30"] #range(1,36) #[5,6,7,8] #ech=[20] #[36]

pathmodel="/cnrm/precip/users/augros/DONNEES/MESONH/CORSE/CT1KM/"
filestart="CT1KM.1.SEG01.0"

# Tmatrix directory
table_ind="" # number of the selected Tmatrix table 
pathTmat="/cnrm/precip/users/augros/Programmes/TMATRIX/DPOLSIMUL/OUTPUT/" 

# Output files
pathfick=pathmodel+'k'+MixedPhase+'/'


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius



         

