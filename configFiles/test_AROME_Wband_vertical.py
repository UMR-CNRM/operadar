#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 *    Example configuration: AROME, W band, vertically pointing radar at the ground at each grid point
 *   -------------------------------------------------------------------------------------------------
 *   (cnst_angle = 90, cloud radar looking up), with the two-way attenuation from the ground (Zh_att).
 *   File: modelFiles/AROME/historic.arome.franmg-01km30+0018:00.fa (31 Aug. 2025, 18 UTC), subdomain over the
 *   south-east of France. Needs epygram and the FA library to read the .fa file.
 *       python -m operadar "historic.arome.franmg-01km30+0018:00.fa" test_AROME_Wband_vertical
 *   Sensitivity: -50 dBZ at 1 km (W band default), increasing with the height above the radar.
"""

# ----- INPUT file path (path to folder containing the files), e.g. "./modelFiles/AROME/" or "./modelFiles/MesoNH/"
input_filePath  = f"./modelFiles/AROME/"
# ----- OUPUT file(s) directory, e.g. "./operadarFiles/AROME/" (keeps the outputs apart from the model files)
output_filePath = f"./operadarFiles/AROME/"

# ----- Lookup tables directory path
path_tables = f"./tables_generator/tables/David2026PhD/"   # default configuration, all bands (David2025AMT: C band only)

# ----- Model name : can be 'Arome' or 'MesoNH'
model = 'Arome'
real_case=True

# ----- Microphysics scheme name : can be 'ICE3', 'ICJW', 'LIMA' or 'LIMC'
#       + a name extension (e.g. 'LIMA_noHail' or 'ICE3_CIBU_moins', optional)
#       Note : only the four first characters are used to select the right table
microphysics_scheme = 'ICE3'

# ----- Number of moments for each hydrometeor of the microphysics scheme
#       /!\ Depending on your simulation, please chose cloud water (cc) over land or sea
hydrometeors_moments = {'cc':1,'rr':1,'ss':1,'gg':1,'ii':1,'wg':1}   # add 'hh':1,'wh':1 for hail and wet hail (L, S, C, X bands)
cloud_water_over='land'     # 'land' or 'sea' --> will read either the cl (land) or cs (sea) lookup table

# ----- Subdomain : written as [lon_min,lon_max,lat_min,lat_max] for a real case
#                           or as [i_min,i_max,j_min,j_max] for an idealized case
#                           or None (will use all the points in the file)
#       e.g. [2.0,4.0,45.0,46.5] (AROME, south-east of France), [53.,60.,-24.,-17.] (MesoNH 3BL13, cyclone near
#       La Reunion). Large MesoNH domains (e.g. 1202 x 1202 x 82) do not fit in ~23 GB of memory: use a subdomain.
subDomain = [2.0,4.0,45.0,46.5]

# ----- Mixed phase simulation : can be 'T_pos' or 'Fw_pos' or 'Fw_posg'.
#                                Please, have a look at the README beforehand.
mixed_phase_parametrization = 'Fw_posg'

# ----- Additional output : if True, compute the dual-pol variables for each hydrometeor class
#                           and save the resultant netcdf (1 file/hydrometeor class)
save_netcdf_single_hydrometeor = False

# ----- Dual-pol variables to add in the output file : provide a list of at least one element 
#                                                      in ['Zh','Zdr','Kdp','Rhohv', 'Ah', 'Av'] 
dpol2add = ['Zh','Zdr','Kdp','Ah']

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

# ---- Compute attenuation for vertical pointing radar
# If compute_attenuation = True, Ah must be included in dpol2add (e.g. dpol2add = ['Zh','Zdr','Kdp','Ah']);
# Zh_att is then added to the output. The attenuation is integrated along the vertical from radar_altitude
# (radar at the ground looking up, or airborne radar), whatever cnst_angle: use cnst_angle = 90 below for
# Zh, Zdr, Kdp consistent with that geometry (with cnst_angle = 0, Zdr and Kdp are those of a horizontal beam).

compute_attenuation = True
radar_altitude = 0 # Radar altitude in meters (can be the aircraft altitude for an onboard radar)


# ----- Radar sensitivity : Zh, Zdr, Kdp and Rhohv are set to NaN where Zh is below the minimum detectable
#       reflectivity Zmin(r) = Zmin_1km + 20 log10(r / 1 km), r = range from the radar (when radarloc is set, or
#       for a vertical pointing radar, cnst_angle=90); constant Zmin_1km when the range is not known.
#       Zmin_1km = None uses a default per band (dBZ at 1 km): L, S, C, X: -30 ; Ku, K, Ka: -40 ; W: -50
sensitivity_threshold = True
Zmin_1km = None

# ----- Radar simulation options 
radar_band = 'W'                    # radar band (C, X, S, W or K)
distmax_rad = 1000.*1000             # maximum radius of the radar data to compute pseudo-observations
radarloc=None                       # radar location: 'center' or [lat_radar,lon_radar] or None (None will not simulate radar beams)
cnst_angle=90                      # will be used when radarloc=None to simulate either a horizontal (cnst_angle=0°)
                                    #  or a vertical (cnst_angle=90°) pointing radar at each grid point
