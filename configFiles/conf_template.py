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

# Model name : can be 'Arome' or 'MesoNH'
model = 'Arome'

# Microphysics scheme name : can be ICE3, ICE4 or LIMA + an accepted name extension (e.g. LIMA_noHail or ICE3_CIBU_moins)
micro_scheme = 'ICE3'

# Number of moments for each hydrometeor of the microphysics scheme
hydrometeors_moments = {'cc':1,'rr':1,'ss':1,'gg':1,'ii':1,'wg':1}

# Subdomain : written as [lon_min,lon_max,lat_min,lat_max] or None (will use all the points in the file)
subDomain = [-2.61,2.1,42.68,46.5]

# Mixed phase simulation : can be 'T_pos' or 'Fw_pos' or 'Fw_posg'. Please, have a look at the README beforehand.
MixedPhase = 'Fw_posg'

# Additional output : if True, compute the dual-pol variables for each hydrometeor class and save the resultant netcdf (1 file/hydrometeor class)
singletype = False

# Tmatrix directory path
path_Tmatrix = "/home/davidcl/Programmation/data/DPOLSIMUL/OUTPUT/"  #"/cnrm/precip/SAVE/davidcl/THESE/DPOLSIMUL/OUTPUT/"

# INPUT file(s) directory
input_directory  = "/home/davidcl/Programmation/data/expeOLIVE/arome/3dvarfr/"  #"/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/" 

# If applicable, to use if experiences are separated in different sub-directories. Else set to ''.
experience_name = 'GOV2'

# OUPUT file(s) directory
outPath = "/home/davidcl/Programmation/output_test/operadar_test/" #"/cnrm/precip/users/davidcl/operadar_test/testing/"


# ----- Radar simulation options ----- #
radar_band = 'C'                       # radar band (C, X, S, W or K)
distmax_rad = 255.*1000                # maximum radius of the radar data to compute pseudo-observations
#alt_max = 15000.      NOT USED ?      # maximum height of the radar data used in the computation of the pseudo-observations
radarloc="center"                      # radar location: 'center' or [lat_radar,lon_radar]


# ===== WILL BE REMOVED LATER ===== #
n_interpol = 32 # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)
RT = 6371.229*10**3 # Earth radius constant
CCIconst=800. # Ice concentration constant for one moment