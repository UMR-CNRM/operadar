#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:56:47 2023

Create single particule scattering table
using pytmatrix or hydroscatt (hail)

input: pytmat/param/TmatParam
output:

@author: augros (based on Jordi Figueras scattering_rain in hydroscatt )
"""

#import datetime
#import argparse
#import atexit
from warnings import warn

import numpy as np
import pandas as pd
import xarray as xr
import cmath

# pytmatrix libs
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, tmatrix_aux, refractive

# hydroscatt libs
import sys
sys.path.insert(0,"./hydroscatt/")
from scattering import compute_scattering_sp_tm
from scattering import compute_scattering_canting_sp
from scattering import compute_angular_moments_analytical
from scattering_io import get_save_dir
from graph import plot_sp_scatt_quantities

# pyscattering libs
import lib_dielectric as diel
import util

#def main():
#    """
#    main
#    """
#    # parse the arguments
#    parser = argparse.ArgumentParser(
#        description='Entry to rain scattering simulations framework')
#
#    # keyword arguments
#    parser.add_argument(
#        '--path', type=str,
#        default='../output/',
#        help='output data path')
#
#    parser.add_argument(
#        '--band', type=str,
#        default='C',
#        help='frequency band. Default C')
#
#    parser.add_argument(
#        '--hydro_type', type=str,
#        default="rain",
#        help='hydrometeor type (rain, snow, graupel, hail. Default rain')
#        
#    parser.add_argument(
#        '--analytical_cant_angl', type=int,
#        default=1,
#        help='If 1 the canting angle will be computed analytically. Default 1')
#
#    args = parser.parse_args()
#
#    print(f'======  single particle scattering simulation started: '
#          f'{datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")}')
#    atexit.register(_print_end_msg,
#     "====== single particle scattering simulation finished: ")

#=============== Parameters =====================

# == Scattering options   
band,wavelength="C",53.2
hydro_type="rain"
analytical_cant_angl=1. #1.
canting_angle = 10. #10.

# == Table options
diam_min,diam_max,diam_step= 0.05,8.,0.05
elev_min, elev_max,elev_step=0.,20.,4.
temp_min, temp_max,temp_step=-20.,40.,1.

# == Plot options
temp_choice=20.
elev_choice=0.

# == Directories
path='../output/'
with_subdirs = False

savedir = get_save_dir(
    path, hydro_type, band, create_dir=True,
    with_subdirs=with_subdirs)

# == Output variables
sp_var_list = ['sca_xsect_h', 'sca_xsect_v', 'im_delta_co','re_delta_co','refl_h', 'zdr', 'kdp', 'rho_hv']
sp_x_var_list = ['diameter']
sp_y_var_list = ['refl', 'zdr','kdp', 'rho_hv']


# ============== Programm ===============================

print(hydro_type)
print('band:', band)
radar_freq=util.freq(wavelength)

diam_list = np.arange(diam_min, diam_max+diam_step, diam_step)
nb_diam = diam_list.size
elev_list = np.arange(elev_min, elev_max+elev_step, elev_step)
nb_elev = elev_list.size
temp_list=np.arange(temp_min, temp_max+temp_step, temp_step)
nb_temp = temp_list.size


# ==== Tables initialization
# Canting option
if analytical_cant_angl:
    ang_moments_dict = compute_angular_moments_analytical(canting_angle)
    fv180 = np.empty([nb_temp,nb_elev,nb_diam], dtype=complex)
    fh180 = np.empty([nb_temp, nb_elev,nb_diam], dtype=complex)
    fv0 = np.empty([nb_temp, nb_elev,nb_diam], dtype=complex)
    fh0 = np.empty([nb_temp, nb_elev,nb_diam], dtype=complex)
else:
    single_part_dict = {'temp': temp_list,'elev': elev_list,'d': diam_list}
    for var in sp_var_list:
        single_part_dict.update({var: np.zeros(nb_temp,nb_elev,nb_diam)})
# end canting option


# ==== Loop over temperatures
for itemp, temp in enumerate(temp_list):

    # Diel constant
    tempK=temp+273.15 #!T (°K)
    
    
    if (hydro_type=="rain"):
        m=cmath.sqrt(diel.QEPSW(tempK, radar_freq))
    
    scatterer = Scatterer(wavelength=wavelength, m=m) 
    
    # Canting option
    if analytical_cant_angl:
        scatterer.orient = orientation.orient_single
    else:
        scatterer.orient = orientation.orient_averaged_fixed
        scatterer.or_pdf = orientation.gaussian_pdf(canting_angle)
    # end canting option

    # === Loop over elevations
    for ielev, elev in enumerate(elev_list):
        # Scatterer geometry initialization
        geom_back = (90.0-elev, 90.0+elev, 0.0, 180.0, 0.0, 0.0)
        geom_forw = (90.0-elev, 90.0-elev, 0.0, 0.0, 0.0, 0.0)


        # ==== Loop over diameters
        for idiam, diam in enumerate(diam_list):
            #print(f'Computing point {idiam} at D={diam}')             
            scatterer.radius = diam/2.
            
            # === Axis ratio
            scatterer.axis_ratio = 1.0/tmatrix_aux.dsr_thurai_2007(diam)
            
            
            if analytical_cant_angl: # only scattering matrix
                scatterer.set_geometry(geom_back)
                s_mat = scatterer.get_S()
                fv180[itemp,ielev,idiam] = 1e-3*s_mat[0][0]
                fh180[itemp,ielev,idiam] = -1e-3*s_mat[1][1]
                scatterer.set_geometry(geom_forw)
                s_mat = scatterer.get_S()
                fv0[itemp,ielev,idiam] = 1e-3*s_mat[0][0]
                fh0[itemp,ielev,idiam] = 1e-3*s_mat[1][1]
            else: # scattering coefficients and variables using pytmatrix
                scatt_sp_dict = compute_scattering_sp_tm(
                    scatterer, geom_back=geom_back, geom_forw=geom_forw,
                    var_list=sp_var_list)
        
                for var in sp_var_list:
                    single_part_dict[var][itemp,ielev,idiam] = scatt_sp_dict[var]
        
        # ====== end loop over diameters
    # ==== end loop over elevations
# ==== end loop over temperatures
        
# Compute scattering coefficients and dual pol variables
if analytical_cant_angl:
    ds_single_part = compute_scattering_canting_sp(
        wavelength,temp_list,elev_list,diam_list,
        fv180, fh180, fv0, fh0, ang_moments_dict,
        var_list=sp_var_list)

else:
    df_single_part = pd.DataFrame.from_dict(single_part_dict)

# save single particle results
fname = (
    f'{savedir}sp_{hydro_type}_{band}_scattering.nc')
ds_single_part.to_netcdf(fname)

print(f'saved {fname}')

# plot single particle results for selected temperature temp_choice, elev_choice
ds_plot=ds_single_part.sel(elevation=elev_choice,temperature=temp_choice)
df_plot=ds_plot.to_pandas()
df_plot['diameter']=df_plot.index

plot_sp_scatt_quantities(
    df_plot, savedir, band, hydro_type, temp_choice, elev_choice, 
     x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)



# end main
# ===========================================================================    
    
#def _print_end_msg(text):
#    """
#    prints end message
#
#    Parameters
#    ----------
#    text : str
#        the text to be printed
#
#    Returns
#    -------
#    Nothing
#
#    """
#    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
#if __name__ == "__main__":
#    main()