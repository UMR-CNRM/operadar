#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 4 2023

@author: augrosc

Configuration file for operad.py
"""

# ==========  Model simulation options ===============
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['vv','cc','rr','ii','ss','gg','hh'] # available model variables
list_types_tot = ['rr','ii','ss','gg','wg','hh','wh']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)


# ==========  Radar options =============== 
band="S"
distmax_rad = 255.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 12000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)
#latrad=
#lonrad=

# ========== Directories / file's name ========= #
# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/cnrm/precip/users/augros/Programmes/TMATRIX/DPOLSIMUL"

# Output directory
pathTmat = repTmat+"/OUTPUT/"


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius