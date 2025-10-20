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



def main(h:str,band:str,method:str,micro:str,nmoment:int,axe:str,dictParam:dict,
         add_ref:bool,
         combine:list,
         invert_column_and_legend:bool,
         ):
    
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
                     invertColLegend=invert_column_and_legend,
                     )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="---------------- LOOKUP TABLES PLOT TOOL ----------------")
    parser.add_argument("hydro", type=str, default='rr',
                        help='Chose one hydrometeor type : rr, ss, gg, wg, wh, hh, cs, cl, ii')
    parser.add_argument("band", type=str, default='C',
                        help='Chose one frequency band : C, S, X, W, K')
    parser.add_argument("--method", type=str, default='Both',
                        help='Which scattering computation method to plot : Tmatrix, Rayleigh, Both. Default : Both')
    parser.add_argument("--micro", type=str, default='ICE3',
                        help='Microphysics scheme name (has to correspond with the one used to name the lookup table). Default : ICE3')
    parser.add_argument("--axe", type=str, default='D',
                        help='The plot horizontal axis : D (diameter) or M (content). Default : D')
    parser.add_argument("--ref", action="store_true",default=False,
                        help='Wether to add or not the base configuration on every subfigure.')
    parser.add_argument("--invert", action="store_true",default=False,
                        help='Invert what is shown in the legend with what is displayed in the column.')
    parser.add_argument("--ARvalue", type=float, nargs='*', default=[],
                        help='Axis ratio : provide as many values as wanted. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--ARfunc", type=str, nargs='*', default=[],
                        help='Axis ratio function : AUds, CNST, BR02, RYdg, RYwg. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--Fw", type=float, nargs='*', default=[],
                        help='Liquid water fraction (only applies to wet hydrometeors) : provide as many values as wanted. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--Nc", type=float, nargs='*', default=[],
                        help='Number concentration (only if --axe M) : provide any number of particles per m3. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--DSTYfunc", type=str, nargs='*', default=[],
                        help='Density function : BR07, RHOX. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--FRIM", type=float, nargs='*', default=[],
                        help='Riming fraction : (only applies when hydrometeor=ss) provide as many values as wanted. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--CANTING", type=float, nargs='*', default=[],
                        help='Canting angle : provide as many values as wanted. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--DIEL", type=str, nargs='*', default=[],
                        help='Dielectric function : Liebe91, RY19dry, LBwetgr, MGwMA08. If not provided, will use the value stored in the base configuration.')
    parser.add_argument("--combine", type=str, nargs='*', default=[],
                        help='Combine at least 2 parameters together so it appears on plot as one curve : axis_ratio, axis_ratio_func, canting_angle, density_func, riming_fraction, diel_func, liquid_water_fraction')
    
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
         invert_column_and_legend=parser_args.invert,
         )