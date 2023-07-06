#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 7 2023
@author: davidcl

"""
import pandas as pd
import datetime as dt

radar_ids = 'OPOU-MCLA-NIME'
run = "00"
deb = pd.Timestamp('2022-08-16 14:00')
fin = pd.Timestamp('2022-08-16 23:45')
step = dt.timedelta(minutes=5)

lat_min,lat_max=42,45
lon_min,lon_max=1,5