#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import gamma
import numpy as np



def get_PSD_parameters_dict(micro_scheme:str,hydrometeor:str,moment:int)-> dict:
    parameters = dict(
        ii = dict(alpha=3 , nu=3 ,a=0.82 , b=2.5 , C=None , x=0),   # pristine ice
        ss = dict(alpha=1 , nu=1 ,a=0.02 , b=1.9 , C=5.0 , x=1),    # snow
        gg = dict(alpha=1 , nu=1 ,a=19.6 , b=2.8 , C=5e5 , x=-0.5), # graupel
        hh = dict(alpha=1 , nu=8 ,a=470. , b=3.0 , C=4e4 , x=-1),   # hail
        rr = dict(alpha=1 , nu=1 ,a=524. , b=3.0 , C=8e6 , x=-1),   # rain
        cc = dict(alpha=3 , nu=1 ,a=524. , b=3.0 , C=300e6 , x=0),  # cloud water (300/cm3 over land 100 over sea)
    )
    if micro_scheme[0:3] == 'LIM' and moment==2 :
        parameters['rr']['nu'] = 2
    
    return parameters[hydrometeor]



def compute_slope_parameter(micro_scheme:str,hydrometeor:str,moment:int,contents:dict[np.ndarray],concentrations:dict[np.ndarray]=None):
    """See p.156 of http://mesonh.aero.obs-mip.fr/mesonh57/BooksAndGuides?action=AttachFile&do=view&target=scidoc_p3.pdf"""
    
    alpha, nu, a, b, C, x = get_PSD_parameters_dict(micro_scheme,hydrometeor,moment).values()
    content = contents[hydrometeor]
    
    if moment == 1 :
        return ( (a*C*gamma(nu + b/alpha)) / (gamma(nu)*content) )^(1/(b-x)) # contents in kg/m3 (or in kg/kg and need to multiply by density)
    elif moment == 2 :
        concentration = concentrations[hydrometeor]
        return ( (a*concentration*gamma(nu + b/alpha)) / (gamma(nu)*content) )^(1/b) # contents in kg/m3 if concentration in m-3
    


def compute_interceipt(micro_scheme:str,hydrometeor:str,moment:int,contents:dict[np.ndarray],concentrations:dict[np.ndarray],temperature:np.ndarray) :
    """Compute the number concentration (N) of 1-moment species according to the MesoNH doc (http://mesonh.aero.obs-mip.fr/mesonh57/BooksAndGuides?action=AttachFile&do=view&target=scidoc_p3.pdf)
    - Cloud water N"cc" : see section 7.2.2 (p.128 of the doc)
    - Pristine ice N"ii" : see section 7.2.3 (p.128 of the doc)
    - Other hydrometeors : N = C*lambda^x
    
    TODO : link doc MNH formula
    """
    alpha, nu, a, b, C, x = get_PSD_parameters_dict(micro_scheme,hydrometeor,moment).values()
    slope_parameter = compute_slope_parameter(micro_scheme, hydrometeor, moment, contents)
    
    if hydrometeor == 'ii':
        Tt = 0 # dans la doc Tt = temperature top atmosphere (p.7) mais aussi Tt = point triple de l'eau = 273.16K (p.86) ? pas clair
        N = np.zeros(temperature.shape)
        if temperature-Tt >= -5 :
            return None
        elif temperature-Tt < -5 :
            return None
    elif hydrometeor == "cc":
        N = 300e6 
    else :    
        N = C * (slope_parameter)^x    
    return N