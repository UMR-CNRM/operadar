#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:20:25 2023

@author: augros
"""

# ================ Constants ==============================
LIGHTSPEED = 299792458. # m/s

# ================ Functions ==============================
def freq(wavelength):
    """
    input = wavelength (mm)
    output = radar frequency (Hz)
    """
    
    radar_freq=LIGHTSPEED*1E3/wavelength
    
    return radar_freq