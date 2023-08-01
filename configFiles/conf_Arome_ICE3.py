#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr 4 2023

@author: augrosc

Configuration file for operad.py
"""

# ==========  Model simulation options =============== #
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['vv','cc','rr','ii','ss','gg'] # available model variables
list_types_tot = ['rr','ii','ss','gg','wg']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)


# ========== Directories / file's name ========= #
# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/cnrm/precip/users/augros/Programmes/TMATRIX/DPOLSIMUL"

# Output directory
pathTmat = repTmat+"/OUTPUT/"


# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius