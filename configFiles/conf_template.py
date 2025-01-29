#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 *   This is a template of the configuration file expected by operad.py 
 *   ------------------------------------------------------------------
 *
 *   You may want to adjust some parameters. 
 *   Please make a copy of this file before any changes
 *
"""
import datetime as dt

# ==========  Characteristics of the data (AROME or MESONH files) ===============
micro_scheme = 'ICE3' # can be ICE3, ICE4 or LIMA
moments={'cc':2,'rr':2,'ss':1,'gg':1,'ii':1,'wg':1}
step = dt.timedelta(minutes=5)
 

# ========== Forward operator options ===============
# Model options
LIMToption=""                               #"" or "cstmu" the model variables are taken from LIMT simulation # but a constant mu is applied in the PSD  for the dpol variables calculation 
n_interpol = 32                             # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)
MixedPhase="Fwposg"                         # 'Tpos' or 'Fwpos' or 'Fwposg'
singletype=False                            # if True: computes dpol var for each type
list_types_tot = ['rr','ii','ss','gg','wg'] # only for T-matrix computation

# Radar simulation options
distmax_rad = 255.*1000     # maximum radius of the radar data to compute pseudo-observations
alt_max = 15000.            # maximum height of the radar data used in the computation of the pseudo-observations
radarloc="center"           # radar location: center or latlon (if latlon ==> to be defined below buy user)


# ========== Directories / file's name ========= #
# Tmatrix directory
table_ind = "" # number of the selected Tmatrix table 
repTmat   = "/cnrm/precip/SAVE/davidcl/THESE/DPOLSIMUL"

# Output directory
pathTmat = repTmat+"/OUTPUT/"

# Model files paths
arome_commonPath  = f"/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/"  # for recurrent nomenclature
arome_commonFilename = "historic.arome.franmg-01km30+00"                    # idem
arome_experience_name = None                                                # set None if multiple
outPath = f"/cnrm/precip/SAVE/davidcl/THESE/operadar_files"
csvPath = "/cnrm/precip/SAVE/davidcl/THESE/studyCases_from_ods.csv"

# ========== Constants =======================
RT = 6371.229*10**3 # Earth radius
CCIconst=800.


# =========== CSV SETUP ==============
""" Please explicit the name of your csv columns based on the given example hereafter

    - Domain :
        * Either your domain is written within one column with the format [lonmin,lonmax,latmin,latmax]
          --> provide the name of the column in the following format :  "name_of_my_column"
        * Or the latitude/longitude min/max are written in 4 distinct columns
          --> provide the name of the columns in the following format : ["lonmin","lonmax","latmin","latmax"]
              (the names can differ, but needs to be written in the same order as the example)
"""
csv_time_columns = ["model_start_time", "model_end_time"]
csv_model_run_column = "model_run" #set to None in case of a MesoNH file
csv_datetime_format = "%Y-%m-%d %H:%M"
csv_delimiter = ";"
csv_domain_columns = "domain"


liste_var_pol = ["Zhh", "Zdr", "Kdp","Rhohv"]
liste_var_calc=["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]