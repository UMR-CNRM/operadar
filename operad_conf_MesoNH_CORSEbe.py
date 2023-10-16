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
#radarloc="latlon" # radar location: center or latlon (if latlon ==> to be defined below buy user)
latrad=float("nan")
lonrad=float("nan")
Radpos = np.array([latrad,lonrad])

# ========== Directories / files name options =========
timelist=["18","20","22","24","26","28"] #range(1,36) #[5,6,7,8] #ech=[20] #[36]

pathmodel="/scratch/work/dricard/CORSE/1km/CT1KM/"
filestart="CT1KM.1.SEG01.0"

# Tmatrix directory
table_ind="" # number of the selected Tmatrix table 
pathTmat="/home/cnrm_other/ge/mrmp/augros/TmatCoefInt_SCXW/"

# Output files
pathfick="/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/CT1KM/"+"dpolvar/"
#pathfick=pathmodel+'k'+MixedPhase+'/'


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius



         

