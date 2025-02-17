#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 *    This is a template of the configuration file
 *   ----------------------------------------------
 *
 *   You may want to adjust some parameters.
 *   Please make a copy of this file before any changes
 *
"""

# ----- INPUT file path (filename will be used to name the output file)
input_filePath  = "home/my_input_folder/" 

# ----- OUPUT file(s) directory
outPath = f"/home/my_output_folder/"

# ----- Tmatrix directory (TmatCoefInt_SCXW) path
path_Tmatrix = "/my_path/TmatCoefInt_SCXW/"

# ----- Model name : can be 'Arome' or 'MesoNH'
model = 'Arome'
   
# ----- Microphysics scheme name : can be 'ICE3', 'ICE4' or 'LIMA'
#       + a name extension (e.g. 'LIMA_noHail' or 'ICE3_CIBU_moins', optional)
micro_scheme = 'LIMA_exp_GNZR'

# ----- Number of moments for each hydrometeor of the microphysics scheme
hydrometeors_moments = {'cc':2,'rr':2,'ss':1,'gg':1,'ii':2,'wg':1}

# ----- Subdomain : written as [lon_min,lon_max,lat_min,lat_max]
#                           or None (will use all the points in the file)
subDomain = None

# ----- Mixed phase simulation : can be 'T_pos' or 'Fw_pos' or 'Fw_posg'.
#                                Please, have a look at the README beforehand.
MixedPhase = 'Fw_posg'

# ----- Additional output : if True, compute the dual-pol variables for each hydrometeor class
#                           and save the resultant netcdf (1 file/hydrometeor class)
save_netcdf_single_hydrometeor = False

# ----- Radar simulation options 
radar_band = 'C'                    # radar band (C, X, S, W or K)
distmax_rad = 255.*1000             # maximum radius of the radar data to compute pseudo-observations
#alt_max = 15000.      NOT USED ?   # maximum height of the radar data used in the computation of the pseudo-observations
radarloc="center"                   # radar location: 'center' or [lat_radar,lon_radar]


# ===== WILL BE REMOVED LATER ===== #
n_interpol = 32 # nb bornes to interpol (2**5: min et max pour LAM, ELEV, T, M, Fw)
RT = 6371.229*10**3 # Earth radius constant
CCIconst=800. # Ice concentration constant for one moment
LIMToption=""