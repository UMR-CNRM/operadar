#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime as dt

# ==========  Model simulation options ===============
LIMToption="" #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
CCIconst=800.
htypes_model=['vv','cc','rr','ii','ss','gg'] # available model variables
n_moments_model = [1,1,1,1,1,1]                 # corresponding moments
list_types_tot = ['rr','ii','ss','gg','wg']
method = ['Tmatrix','Tmatrix','Tmatrix','Tmatrix','Tmatrix']

MixedPhase="Fwposg" # 'Tpos' or 'Fwpos' or 'Fwposg' #
 
singletype=False #False #True # if True: computes dpol var for each type

n_interpol = 32      # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)

step = dt.timedelta(minutes=5)

# Radar options
distmax_rad = 255.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 15000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)


# ========== Directories / file's name ========= #
# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/home/davidcl/Programmation/data/DPOLSIMUL"

# Output directory
pathTmat = repTmat+"/OUTPUT/"

# Model files paths
commonPath  = "/home/davidcl/Programmation/data/expeOLIVE/arome/3dvarfr/"
arome_experience_name = None
commonFilename = "historic.arome.franmg-01km30+00"
outPath = "/home/davidcl/Programmation/output_test/operadar_test"
csvPath = "/home/davidcl/Programmation/data/testCase_convective.csv"

# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius

csv_time_columns = ["model_start_time", "model_end_time"]
csv_model_run_column = "model_run" #set to None in case of a MesoNH file
csv_datetime_format = "%Y-%m-%d %H:%M"
csv_delimiter = ";"
csv_domain_columns = "model_domain" # or ["latmin","latmax","lonmin","lonmax"]

liste_var_pol = ["Zhh", "Zdr", "Kdp","Rhohv"]
dpol_var_to_calc=["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]