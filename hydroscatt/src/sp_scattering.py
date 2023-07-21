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

import datetime
import argparse
import atexit
from warnings import warn

import numpy as np
import pandas as pd

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, tmatrix_aux, refractive

#from scattering_io import get_save_dir
#from graph import plot_sp_scatt_quantities, plot_psd_scatt_quantities
#from part_descrip import compute_velocity_rain
#from precip import compute_lwc, compute_rainfall_rate
from scattering import compute_scattering_psd_tm, compute_scattering_sp_tm
from scattering import compute_scattering_canting_sp
from scattering import compute_scattering_canting_psd
from scattering import compute_angular_moments_analytical

# parameters
sp_var_list = [
    'sca_xsect_h', 'sca_xsect_v', 'refl_h', 'refl_v', 'ldr_h', 'ldr_v',
    'zdr', 'rho_hv', 'delta_hv', 'ext_xsect_h', 'ext_xsect_v', 'kdp',
    'A_h', 'A_v', 'Adp']

sp_x_var_list = ['d']
sp_y_var_list = [
    'sca_xsect', 'ext_xsect', 'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv',
    'kdp', 'A', 'Adp']

psd_var_list = [
    'refl_h', 'refl_v', 'ldr_h', 'ldr_v', 'zdr', 'rho_hv', 'delta_hv',
    'kdp', 'A_h', 'A_v', 'Adp']

psd_x_var_list = ['refl_h', 'lwc', 'rr', 'D0']
psd_y_var_list = [
    'refl', 'ldr', 'zdr', 'rho_hv', 'delta_hv', 'kdp', 'A', 'Adp']

diam_min = 0.1
diam_max = 7.
step = 0.1
diam = np.arange(diam_min, diam_max+step, step)
num_points = diam.size

canting_angle = 10.

# geometry = (theta0, theta, phi0, phi, alpha, beta)
# geom_horiz_back = (90.0, 90.0, 0.0, 180.0, 0.0, 0.0) #horiz. backscatter
# geom_horiz_forw = (90.0, 90.0, 0.0, 0.0, 0.0, 0.0) #horiz. forward scatter
# geom_vert_back = (0.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. backscatter
# geom_vert_forw = (180.0, 180.0, 0.0, 0.0, 0.0, 0.0) #vert. forward scatter
geom_back = (90.0-args.ele, 90.0+args.ele, 0.0, 180.0, 0.0, 0.0)
geom_forw = (90.0-args.ele, 90.0-args.ele, 0.0, 0.0, 0.0, 0.0)

print(hydro_type)
print('band:', args.band)
print('temp:', args.temp)
print('elevation angle', args.ele)

if args.band == 'S':
    wavelength = tmatrix_aux.wl_S
elif args.band == 'C':
    wavelength = tmatrix_aux.wl_C
elif args.band == 'X':
    wavelength = tmatrix_aux.wl_X

if args.temp == 0.:
    m = refractive.m_w_0C[wavelength]
elif args.temp == 10.:
    m = refractive.m_w_10C[wavelength]
elif args.temp == 20.:
    m = refractive.m_w_20C[wavelength]

scatterer = Scatterer(wavelength=wavelength, m=m)
if args.analytical_cant_angl:
    scatterer.orient = orientation.orient_single
    ang_moments_dict = compute_angular_moments_analytical(canting_angle)
else:
    scatterer.orient = orientation.orient_averaged_fixed
    scatterer.or_pdf = orientation.gaussian_pdf(canting_angle)

savedir = get_save_dir(
    args.path, hydro_type, args.band, args.temp, create_dir=True,
    with_subdirs=with_subdirs)

if args.compute_sp:
    # single particle scattering
    if args.analytical_cant_angl:
        fv180 = np.empty(num_points, dtype=complex)
        fh180 = np.empty(num_points, dtype=complex)
        fv0 = np.empty(num_points, dtype=complex)
        fh0 = np.empty(num_points, dtype=complex)
    else:
        single_part_dict = {'d': diam}
        for var in sp_var_list:
            single_part_dict.update({var: np.zeros(num_points)})

    for ind, d_part in enumerate(diam):
        print(f'Computing point {ind} at D={d_part}')

        scatterer.radius = d_part/2.
        scatterer.axis_ratio = 1.0/tmatrix_aux.dsr_thurai_2007(d_part)
        if args.analytical_cant_angl:
            scatterer.set_geometry(geom_back)
            s_mat = scatterer.get_S()
            fv180[ind] = 1e-3*s_mat[0][0]
            fh180[ind] = -1e-3*s_mat[1][1]
            scatterer.set_geometry(geom_forw)
            s_mat = scatterer.get_S()
            fv0[ind] = 1e-3*s_mat[0][0]
            fh0[ind] = 1e-3*s_mat[1][1]
        else:
            scatt_sp_dict = compute_scattering_sp_tm(
                scatterer, geom_back=geom_back, geom_forw=geom_forw,
                var_list=sp_var_list)

            for var in sp_var_list:
                single_part_dict[var][ind] = scatt_sp_dict[var]

    if args.analytical_cant_angl:
        df_single_part = compute_scattering_canting_sp(
            wavelength, fv180, fh180, fv0, fh0, ang_moments_dict,
            var_list=sp_var_list)
        df_single_part['d'] = diam
    else:
        df_single_part = pd.DataFrame.from_dict(single_part_dict)

    # save single particle results
    fname = (
        f'{savedir}sp_{hydro_type}_{args.band}_{int(args.temp*100):04d}'
        f'_ele{int(args.ele*100.):05d}_scattering.csv')
    df_single_part.to_csv(fname, index=False)
    print(f'saved {fname}')

    plot_sp_scatt_quantities(
        df_single_part, savedir, args.band, args.temp, hydro_type,
        ele=args.ele, x_var_list=sp_x_var_list, y_var_list=sp_y_var_list)