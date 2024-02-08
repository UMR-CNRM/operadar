#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:00:34 2023

@author: augros
"""


import numpy as np
import xarray as xr


# ================== Save Arome variables ===================================
"""
Save Arome dpol variables in npz or netcdf file
input: Vm_k, Tc, Z, lat, lon , fick
output: file saved
"""
def save_dpolvar_arome(liste_var_pol,M, CC, CCI, Vm_k, Tc, Z,lat,lon,fick,time,
                       save_npz=True, save_netcdf=True):
    
    IKE=Tc.shape[0] 
    IJE=Tc.shape[1]
    IIE=Tc.shape[2]
    
    # ============= Save in npz file : tab_fick 
    if save_npz :
        nb_lin=IKE*IJE*IIE

        # dtype = "[('i',int),('j',int),('k',int),('Zhh', float), ('Zdr', float), ('Kdp', float), ('Rhohv', float)]"
        # cmd = "tab_fick = np.zeros(nb_lin, dtype = "+dtype+")"
        # exec(cmd)
        dtype = [('i',int),('j',int),('k',int),('Zhh', float), ('Zdr', float), ('Kdp', float), ('Rhohv', float)]
        tab_fick = np.zeros(nb_lin, dtype = dtype)
        
        index=np.array([[[(k,j,i) for i in range(IIE)] for j in range(IJE)] for k in range(IKE)])
            
        tab_fick['i'] = (np.rollaxis(np.rollaxis(index[:,:,:,2], 2, 0), 2, 1)).ravel()
        tab_fick['j'] = (np.rollaxis(np.rollaxis(index[:,:,:,1], 2, 0), 2, 1)).ravel()
        tab_fick['k'] = (np.rollaxis(np.rollaxis(index[:,:,:,0], 2, 0), 2, 1)).ravel()  
        for var in liste_var_pol:
            tab_fick[var] = (np.rollaxis(np.rollaxis(Vm_k[var], 2, 0), 2, 1)).ravel()
        
        print ("Writing varpol in fick :")
        print (fick)
        np.savez_compressed(fick+'.npz', tab_fick)
  
    # ============ Save var pol in netcdf file with xarray
    
    # M dict formatting for dataset
    hydromet_list = list(M.keys())
    contents = np.array([M[hydromet]*1000 for hydromet in hydromet_list]).astype('f4') # from kg to g/kg
    
    if save_netcdf :
        ds=xr.Dataset(
                data_vars=dict(
                        Zh     = (["level","y","x"],Vm_k["Zhh"].astype('f4'), {"units": "dBZ"}),
                        Zdr    = (["level","y","x"],Vm_k["Zdr"].astype('f4'), {"units": "dB"}),
                        Kdp    = (["level","y","x"],Vm_k["Kdp"].astype('f4'), {"units": "°/km"}),
                        Rhohv  = (["level","y","x"],Vm_k["Rhohv"].astype('f4'), {"units": "1"}),
                        M      = (["hydrometeor","level","y","x"],contents, {"units": "g/kg of dry air"}),
                        CCrain = (["level","y","x"],CC.astype('f4'), {"units": "kg^-1"}),
                        CCice  = (["level","y","x"],CCI.astype('f4'), {"units": "kg^-1"}),
                        T      = (["level","y","x"],Tc.astype('f4'), {"units": "°C"}),
                        Alt    = (["level","y","x"],Z.astype('i4'), {"units": "m"}),
                ),
                coords=dict(
                    y   = (["y"], np.arange(Tc.shape[1]).astype('i4')),
                    x   = (["x"], np.arange(Tc.shape[2]).astype('i4')),
                    lon = (["y","x"], lon.astype('f4')),
                    lat = (["y","x"], lat.astype('f4')),
                    level = (["level"], np.arange(Tc.shape[0]).astype('i4')),
                    time  = (time),             
                    hydrometeor = (["hydrometeor"],hydromet_list),
        ),
        )    
        ds.to_netcdf(fick+".nc")
        ds.close() ; del ds
    
    # =========== Plot Zh at first level to test ==================
    #ds.Zh.sel(level=89).plot(x="lon",y="lat",cmap="viridis",vmin=0)   

# ============================================================================    


# ================== Save MesoNH variables ===================================

"""
Save dpol variables in npz or netcdf file
input: Vm_k, Tc, Z, lat, lon (Arome), fick
output: file saved
"""
def save_dpolvar_mesonh(liste_var_pol, Vm_k, Tc, Z, X, Y,lat,lon,time,fick,
                        save_npz=True, save_netcdf=True):
    
    IKE=Tc.shape[0]
    IJE=Tc.shape[1]
    IIE=Tc.shape[2]
    
    # ============= Save in npz file : tab_fick 
    if save_npz :
        nb_lin=IIE*IJE*IKE
        dtype = [('i',int),('j',int),('k',int),('Zhh', float), ('Zdr', float), ('Kdp', float), ('Rhohv', float)]
        tab_fick = np.zeros(nb_lin, dtype = dtype)
        
        index=np.array([[[(k,j,i) for i in range(IIE)] for j in range(IJE)] for k in range(IKE)])
            
        tab_fick['i'] = (np.rollaxis(np.rollaxis(index[:,:,:,2], 2, 0), 2, 1)).ravel()
        tab_fick['j'] = (np.rollaxis(np.rollaxis(index[:,:,:,1], 2, 0), 2, 1)).ravel()
        tab_fick['k'] = (np.rollaxis(np.rollaxis(index[:,:,:,0], 2, 0), 2, 1)).ravel()  
        for var in liste_var_pol:
            tab_fick[var] = (np.rollaxis(np.rollaxis(Vm_k[var], 2, 0), 2, 1)).ravel() 
        
        print ("Writing varpol in fick :")
        print (fick)
        np.savez_compressed(fick+'.npz', tab_fick)
      
    # ============ Save var pol in netcdf file with xarray
    if save_netcdf :
        ds=xr.Dataset(
                data_vars=dict(
                        Zh=(["level","y","x"],Vm_k["Zhh"]),
                        Zdr=(["level","y","x"],Vm_k["Zdr"]),
                        Kdp=(["level","y","x"],Vm_k["Kdp"]),
                        Rhohv=(["level","y","x"],Vm_k["Rhohv"]),
                        T=(["level","y","x"],Tc),
                        Alt=(["level"],Z),                
                            #             "Kdp":(["lat","lon","alt"],Vm_k["Kdp"]),
        #             "Rhohv":(["lat","lon","alt"],Vm_k["Rhohv"]),},
                ),
		 coords=dict(
		           lon=(["x"],lon),
		           lat=(["y"],lat),
		           Time = (["time"], time),
		           level=(["level"], np.arange(Z.shape[0])),
		           #Radloc = (["radpos"],Radpos),
	 ),
        )    
        ds.to_netcdf(fick+".nc")
        ds.close() ;
        del ds

    
    # # =========== Plot Zh at first level to test (level 2 = 10 m in MesoNH)=======
    # ds.Zh.sel(level=2).plot(x="lon",y="lat",cmap="viridis",vmin=0)
  

# ============================================================================         
