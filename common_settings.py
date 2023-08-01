#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 7 2023
@author: davidcl

"""
import datetime as dt

step = dt.timedelta(minutes=5)
save_npz    = False
save_netcdf = True

# Model files paths
commonPath_fa  = f"/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/"
commonFilename = "historic.arome.franmg-01km30+00" #08:00.fa"
outPath = f"/cnrm/precip/SAVE/davidcl/THESE/operadar_files"
csvPath = "./expe_olive.csv"

# Radar options
distmax_rad = 255.*1000 #150*1000 # Distance max des données radar dont on calcule les pseudo-observations
alt_max = 15000. # Altitude max des données radar utilisées pour le calcul des pseudo-observations
radarloc="center" # radar location: center or latlon (if latlon ==> to be defined below buy user)
#latrad=
#lonrad=