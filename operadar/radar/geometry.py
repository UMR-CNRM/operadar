#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import time as tm
import operadar.operad_conf as cf



def compute_radar_geometry(X:np.ndarray, Y:np.ndarray, Z:np.ndarray,
                           Tc:np.ndarray, elev_max:float,
                           distmax_rad:float=cf.distmax_rad)-> tuple[np.ndarray,np.ndarray]:
    """Mimic radar geometry in the model grid.

    Args:
        X (np.ndarray): 1D horizontal coordinates
        Y (np.ndarray): 1D horizontal coordinates
        Z (np.ndarray): 1D vertical coordinates
        Tc (np.ndarray): temperature field (3D)
        elev_max (np.ndarray): maximum elevation angle
        distmax_rad (float, optional): maximum radius of the radar data to compute pseudo-observations. Defaults to cf.distmax_rad.

    Returns:
        distance_mask (np.ndarray) : mask of the grid points within the maximum radar radius (3D).
        elev (np.ndarray) : corresponding elevation value for each grid point (3D).
    """

    print("Compute radar geometry") ; deb_timer = tm.time()
    np.seterr(invalid='ignore') # silence warning of invalid division 0 by 0 (result in a nan)
    
    if (cf.model=="Arome"):
        print('\tNo radar geometry computation implemented yet in Arome.')
        elev = np.zeros(Tc.shape)
        distance_mask = (elev >= 0.)
        
    elif (cf.model=="MesoNH"):
        if type(cf.radarloc)==str and cf.radarloc=="center" :
            print('\tComputing radar geometry at the center of the domain.')
            X0, Y0, Z0 = np.nanmean(X), np.nanmean(Y), 0.
        elif type(cf.radarloc)==list :
            radar_lat = cf.radarloc[0]
            radar_lon = cf.radarloc[1]
            print(f'\tComputing radar geometry at latitude {radar_lat}° and longitude {radar_lon}°.')
            
        distance_mask, radar_dist_3D = compute_distance_mask(X, Y, Z, X0, Y0, Z0)
        elev = compute_radar_elevation(radar_dist_3D, Z, elev_max)
        
    return distance_mask, elev
    
    
    
def compute_distance_mask(X:np.ndarray, Y:np.ndarray, Z:np.ndarray,
                          X0:float, Y0:float, Z0:float,
                          distmax_rad:float=cf.distmax_rad)-> tuple[np.ndarray,np.ndarray] :
    """Mask grid points that are beyond `distmax_rad` in config file."""
    XX,YY=np.meshgrid(X-X0,Y-Y0)
    radar_dist=(XX**2+YY**2)**0.5 # 2D array with distance from radar for each grid point
    radar_dist_3D= np.stack([radar_dist]*Z.shape[0],axis=0) 
    mask_distmax=(radar_dist_3D < distmax_rad)
    return mask_distmax, radar_dist_3D



def compute_radar_elevation(radardist3D:np.ndarray, Z:np.ndarray, elev_max:float)-> np.ndarray:
    """Compute the radar elevation angle value at each grid point"""
    earth_radius = 6371.229e3
    tanel = Z/radardist3D - 3.*radardist3D/(8.*earth_radius)
    elev = np.arctan(tanel)*180./math.pi
    elev[elev<0] = 0.
    elev[elev>elev_max] = elev_max
    elev[:,:,:]=0.
    return elev