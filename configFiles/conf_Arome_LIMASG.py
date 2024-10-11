#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime as dt

# ==========  Model simulation options ===============
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['vv','cc','rr','ii','ss','gg'] # available model variables
list_types_tot = ['rr','ii','ss','gg','wg']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)

step = dt.timedelta(minutes=5)
save_npz    = False # save file(s) in .npz format
save_netcdf = True # save file(s) in .nc format (netcdf4)

# Radar options
distmax_rad = 255.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 15000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)


# ========== Directories / file's name ========= #
# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/cnrm/precip/SAVE/davidcl/THESE/DPOLSIMUL"

# Output directory
pathTmat = repTmat+"/OUTPUT/"

# Model files paths
commonPath  = f"/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/"
commonFilename = "historic.arome.franmg-01km30+00"
outPath = f"/cnrm/precip/SAVE/davidcl/THESE/operadar_files_test"
csvPath = "/cnrm/precip/SAVE/davidcl/THESE/testCase.csv"

# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius