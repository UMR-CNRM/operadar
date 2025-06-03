#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 *    This is a template of the configuration file
 *   ----------------------------------------------
 *
 *   You may want to adjust some parameters.
 *   /!\ Possible options are displayed. Be careful and respect if quotes are needed or not.
 *   Please make a copy of this file before any changes
 *
"""

# ----- INPUT file path (path to folder containing the files)
input_filePath  = "home/my_input_folder/" 

# ----- OUPUT file(s) directory
outPath = f"/home/my_output_folder/"

# ----- Lookup tables directory path
path_tables = "./tmatrix_generator/tables/default/"

# ----- Model name : can be 'Arome' or 'MesoNH'
model = 'Arome'
real_case=False   # for MesoNH only

# ----- Microphysics scheme name : can be 'ICE3', 'ICE4' or 'LIMA'
#       + a name extension (e.g. 'LIMA_noHail' or 'ICE3_CIBU_moins', optional)
#       Note : only the three first characters are used to select the right scheme
micro_scheme = 'LIMA_exp_GNZR'

# ----- Number of moments for each hydrometeor of the microphysics scheme
hydrometeors_moments = {'cc':2,'rr':2,'ss':1,'gg':1,'ii':2,'wg':1}

# ----- Subdomain : written as [lon_min,lon_max,lat_min,lat_max] for a real case
#                           or as [i_min,i_max,j_min,j_max] for an idealized case
#                           or None (will use all the points in the file)
subDomain = None

# ----- Mixed phase simulation : can be 'T_pos' or 'Fw_pos' or 'Fw_posg'.
#                                Please, have a look at the README beforehand.
MixedPhase = 'Fw_posg'

# ----- Additional output : if True, compute the dual-pol variables for each hydrometeor class
#                           and save the resultant netcdf (1 file/hydrometeor class)
save_netcdf_single_hydrometeor = False

# ----- Dual-pol variables to add in the output file : provide a list of at least one element 
#                                                      in ['Zh','Zdr','Kdp','Rhohv', 'Ah', 'Av'] 
dpol2add = ['Zh','Zdr','Kdp']

# ----- Scattering method (TO COME)
#       The user can specify if one method ('Tmatrix' or 'Rayleigh') or 'both' methods are employed to 
#       compute the polarimetric variables, or, specify the method for each hydrometeor individually.
#       See the examples below.
#       * scattering_method = 'both' will create two output files instead of one, respectively named *_Tmatrix.nc and *_Rayleigh.nc
#                                    + option --append --> will instead append fieldName_Tmatrix and fieldName_Rayleigh in the input file
#       * scattering_method = 'Tmatrix' or 'Rayleigh' will create one output file named *_{scattering_method}.nc
#                                                     + option --append --> will instead append fieldName_{scattering_method} in the input file
#       * scattering_method = {'cc':'Rayleigh','rr':'Rayleigh','ss':'Tmatrix','gg':'Tmatrix','ii':'Tmatrix','wg':'Tmatrix'}
#                             --> will create a unique file with the polarimetric fields resulting from the combinations
#                                 of the chosen methods in this dictionnary
scattering_method = "Tmatrix"

# ----- Radar simulation options 
radar_band = 'C'                    # radar band (C, X, S, W or K)
distmax_rad = 255.*1000             # maximum radius of the radar data to compute pseudo-observations
radarloc=None                       # radar location: 'center' or [lat_radar,lon_radar] or None (None will not simulate radar beams)
cnst_angle=90                       # will be used when radarloc=None to simulate either a horizontal (cnst_angle=0°)
                                    #  or a vertical (cnst_angle=90°) pointing radar at each grid point
