#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:00:34 2023

@author: augros
"""


import numpy as np
import xarray as xr

# ================== Save model (Arome or MesoNH) + simulated radar variables ===================================

"""
Save dpol variables in npz or netcdf file
input: M, CC, CCI,Vm_k, Tc, Z, X, Y, lat, lon (Arome), fick
output: file saved
"""
def save_dpolvar(M, Nc, Vm_k, Tc, Z, X, Y,lat,lon,datetime,outfile):
    
    # M dict formatting for dataset
    hydromet_list = list(M.keys())
    Nc_list = list(Nc.keys())
    contents = np.array([M[hydromet]*1000 for hydromet in hydromet_list]).astype('f4') # from kg to g/m3
    concentrations = np.array([Nc[hydromet] for hydromet in Nc_list]).astype('f4') # from kg to g/m3
    
    ds=xr.Dataset(
            data_vars=dict(
                Zh     = (["level","y","x"],Vm_k["Zhh"].astype('f4'), {"units": "dBZ"}),
                Zdr    = (["level","y","x"],Vm_k["Zdr"].astype('f4'), {"units": "dB"}),
                Kdp    = (["level","y","x"],Vm_k["Kdp"].astype('f4'), {"units": "°/km"}),
                Rhohv  = (["level","y","x"],Vm_k["Rhohv"].astype('f4'), {"units": "1"}),
                M      = (["hydrometeor","level","y","x"],contents, {"units": "g/m3"}),
                Nc     = (["hydrometeor","level","y","x"],concentrations, {"units": "kg^-1"}),
                T      = (["level","y","x"],Tc.astype('f4'), {"units": "°C"}),
                Alt    = (["level","y","x"],Z.astype('i4'), {"units": "m"}),
                ),
                coords=dict(
                   y   = (["y"], Y.astype('f4')),
                   x   = (["x"], X.astype('f4')),
                   lon = (["y","x"], lon.astype('f4')),
                   lat = (["y","x"], lat.astype('f4')),
                   level=(["level"], np.arange(Z.shape[0]).astype('i4')),
                   hydrometeor = (["hydrometeor"],hydromet_list),
		           time = (datetime),
		           #Radloc = (["radpos"],Radpos),
	 ),
    )    
    ds.to_netcdf(outfile+".nc")
    ds.close() ;
    del ds
    print("Saving model and dpol variables in: ",outfile+".nc")

    
    # # =========== Plot Zh at first level to test (level 2 = 10 m in MesoNH)=======
    # ds.Zh.sel(level=2).plot(x="lon",y="lat",cmap="viridis",vmin=0)
  

# ============================================================================         
