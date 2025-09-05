# -*- coding: utf-8 -*-
"""
Creation 03 March 2025
@author: davidcl
"""

import argparse
from plot_tools.utils.sensitivity_test import *


Path_tables = "/home/davidcl/Programmation/operadar/tables_generator/tables/"
dir_fig="/home/davidcl/Programmation/operadar/plot_tools/sensitivity_study/"

Fw_list,Fw_ls=[0.0,0.1,0.6,1.0],['-.',':','--','-']

Fwchoix=0
ELEVchoix=0 #0 pour radars sol, 90 pour rasta
Nii=800 #selected number concentration for primary ice



def main(h:str,band:str,method:str,micro:str,nmoment:int,axe:str,dictParam:dict,add_ref:bool,combine:list):
    
    sensitivity_test(hydrometeor=h,
                     axeX=axe,
                     dictParam=dictParam,
                     which_dpolVar=['Zh','Zdr','Kdp'],#,'Rhohv'],
                     method=method,
                     microphysics=micro,
                     moment=nmoment,
                     band=band,
                     folder_tables=Path_tables,
                     folder_figures=dir_fig,
                     ref=add_ref,
                     combine=combine,
                     )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="---------------- LOOKUP TABLES PLOT TOOL ----------------")
    parser.add_argument("hydro", type=str, default='rr',
                        help='Hydrometeor type : rr, ss, gg, wg, wh, hh, cc, ii')
    parser.add_argument("band", type=str, default='C',
                        help='Band type : C, S, X, W, K')
    parser.add_argument("--method", type=str, default='Both',
                        help='Method : Tmatrix, Rayleigh, Both')
    parser.add_argument("--micro", type=str, default='ICE3',
                        help='Microphysics scheme : ICE3, LIMA')
    parser.add_argument("--axe", type=str, default='D',
                        help='Horizontal axis : D (diameter) or M (content)')
    parser.add_argument("--ref", action="store_true",default=False,
                        help='Wether to add or not the base configuration on every subfigure.')
    parser.add_argument("--ARvalue", type=float, nargs='*', default=[],
                        help='Axis ratio : provide as many values as wanted.')
    parser.add_argument("--ARfunc", type=str, nargs='*', default=[],
                        help='Axis ratio function : AUds, CNST, BR02, RYdg, RYwg')
    parser.add_argument("--Fw", type=float, nargs='*', default=[],
                        help='Liquid water fraction (only for wet hydrometeors) : provide as many values as wanted.')
    parser.add_argument("--Nc", type=float, nargs='*', default=[],
                        help='Number concentration (only if --axe M) : provide any number of particles per m3.')
    parser.add_argument("--DSTYfunc", type=str, nargs='*', default=[],
                        help='Density function : BR07, RHOX')
    parser.add_argument("--FRIM", type=float, nargs='*', default=[],
                        help='Riming fraction : (only applies when DSTYfunc=BR07) provide as many values as wanted.')
    parser.add_argument("--CANTING", type=float, nargs='*', default=[],
                        help='Canting angle : provide as many values as wanted.')
    parser.add_argument("--DIEL", type=str, nargs='*', default=[],
                        help='Dielectric function : Liebe91, RY19dry, LBwetgr, MGwMA08')
    parser.add_argument("--combine", type=str, nargs='*', default=[],
                        help='Combine at least 2 parameters together so it appears on plot as one curve : axis_ratio, axis_ratio_func, canting_angle, density_func, riming_fraction, diel_func')
    
    parser_args = parser.parse_args()
    
    dictArgs = dict(axis_ratio = parser_args.ARvalue,
                    axis_ratio_func = parser_args.ARfunc,
                    canting_angle = parser_args.CANTING,
                    density_func = parser_args.DSTYfunc,
                    riming_fraction = parser_args.FRIM,
                    diel_func = parser_args.DIEL,
                    liquid_water_fraction = parser_args.Fw,
                    number_concentration = parser_args.Nc,
                    )
    if len(parser_args.Nc)>=1 : nMoment=2
    else : nMoment = 1
    main(h=parser_args.hydro,
         band=parser_args.band,
         method=parser_args.method,
         micro=parser_args.micro,
         nmoment=nMoment,
         axe=parser_args.axe,
         dictParam = dictArgs,
         add_ref=parser_args.ref,
         combine=parser_args.combine,
         )