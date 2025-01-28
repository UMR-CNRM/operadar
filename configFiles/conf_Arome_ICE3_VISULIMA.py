#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 15 2023
@author: augrosc et davidcl

Configuration file for operad.py
"""
import datetime as dt

# ====  Microphysics scheme options 
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['rr','ss','gg','ii'] # available hydrometeors (model fields)
htypes_table = ['rr','ss','gg','ii','wg']
moments={'rr':1,'ss':1,'gg':1,'ii':1,'wg':1}

# ==== Subdomain
subDomain=False 

# ==== Real or Ideal case
real_case=True

# ==== Tmatrix options
MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)

# ==== Directories / files name options 
commonPath  = "/scratch/work/straussc/OPE_RADAR/ICE3/"
commonFilename = "historic.arome.franmg-01km30+00" #08:00.fa"
outPath = "/home/cnrm_other/ge/mrmp/augros/WKD/VISULIMA/dpolvar/ICE3/"
csvPath = "./study_cases/VISULIMA.csv"
pathTmat="/home/cnrm_other/ge/mrmp/augros/TmatCoefInt_SCXW/"

# ==== Forward operator options
singletype=False #False #True # if True: computes dpol var for each type
step = dt.timedelta(hours=1)
step_seconds = 3600.
save_npz    = False
save_netcdf = True

# ==== Radar options
distmax_rad = 1000.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 15000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)
#latrad=
#lonrad=


# ==== Constants 
RT = 6371.229*10**3 # Earth radius
