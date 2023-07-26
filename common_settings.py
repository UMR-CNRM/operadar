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