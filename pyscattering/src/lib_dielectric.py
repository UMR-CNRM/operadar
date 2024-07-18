#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:59:07 2023

Dielectric functions (python conversion from Tmatrix.f and MesoNH subroutines)

@author: augros
"""




# ================ Dielectric functions ===================

def QEPSW(PTEMP, PFREQ):
    # water complex dielectric function (Liebe et al., 1991)
    # electromagnetic fields in exp(-i*omega*t), i.e. Im(epsw)>=0
    # in  : ptemp=temperature in K
    #       pfreq=frequency in Hz
    # out : eps=epsilon
    ZTHETA = 1 - 300. / PTEMP
    ZEZ = 77.66 - 103.3*ZTHETA
    ZEINF = 0.066*ZEZ
    ZF = (20.27 + 146.5*ZTHETA + 314.*ZTHETA**2)*1.E9
    EPSW = ZEINF + (ZEZ - ZEINF) / (1. - complex(0., 1.)*PFREQ/ZF)
    return EPSW