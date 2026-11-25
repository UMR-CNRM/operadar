#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: davidcl
"""

import sys
import numpy as np
from netCDF4 import Dataset
from typing import Sequence
from operadar.utils.formats_data import get_lat_lon_from_subdomain



def check_variable_is_in_dataset(mnhFile:Dataset):
    necessary_variables = ['XHAT','YHAT','ZHAT','PABST','THT','RHOREFZ']
    variables_list = mnhFile.variables.keys()
    missing_var = []
    for var in necessary_variables :
        if var not in variables_list :
            missing_var += [var]

    if len(missing_var)>0 :
        print('_____________')
        print('/!\\ ERROR /!\\ Missing variables :',missing_var)
        sys.exit()



def get_geometry(mnhFile:Dataset,
                 real_case:bool,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 ):
    X = mnhFile.variables['XHAT'][:][j_min:j_max]
    Y = mnhFile.variables['YHAT'][:][i_min:i_max]
    ZHAT = mnhFile.variables['ZHAT'][:] # height above ground (m)
    if real_case :
        ZS = mnhFile.variables['ZS'][:][i_min:i_max,j_min:j_max] # altitude (m) of model surface (ground)
        LAT = mnhFile.variables['latitude'][:][i_min:i_max,j_min:j_max]
        LON = mnhFile.variables['longitude'][:][i_min:i_max,j_min:j_max]
    else :
        ZS = np.zeros((X.shape[0],Y.shape[0])) # is null for idealized cases
        LAT =float('nan')
        LON = float('nan')
    Z_3D = np.empty((ZHAT.shape[0],ZS.shape[0],ZS.shape[1]))
    for level in range(ZHAT.shape[0]):
        Z_3D[level,:,:]=ZS[:,:]+ZHAT[level]
    return X, Y, Z_3D, LAT, LON



def get_subdomain_indices(mnhFile:Dataset,
                          subDomain:Sequence[float]|Sequence[int],
                          real_case:bool,
                          ) -> tuple[int,int,int,int]:
    if real_case:
        LAT = mnhFile.variables['latitude'][:]
        LON = mnhFile.variables['longitude'][:]
        lon_min, lon_max, lat_min, lat_max = get_lat_lon_from_subdomain(subDomain)
        mask_zoom = ((LON>lon_min) &(LON<lon_max) & (LAT>lat_min) & (LAT<lat_max) )
        [ilon,jlat]=np.where(mask_zoom)
        # +1: the indices are used as slice ends (exclusive), the last row/column inside the box is kept
        i_min, i_max = np.nanmin(ilon), np.nanmax(ilon) + 1
        j_min, j_max = np.nanmin(jlat), np.nanmax(jlat) + 1
    else :
        i_min, i_max = subDomain[0], subDomain[1]
        j_min, j_max = subDomain[2], subDomain[3]
    return i_min, i_max, j_min, j_max



def get_pressure_temperature_density(mnhFile:Dataset,
                                i_min:int, i_max:int,
                                j_min:int, j_max:int,
                                ):
    # Pressure
    p=mnhFile.variables['PABST'][0,:,:,:][:,i_min:i_max,j_min:j_max]
    p[np.where(p==999.)] = float('nan')
    # Temperature
    Th=mnhFile.variables['THT'][0,:,:,:][:,i_min:i_max,j_min:j_max]
    Th[np.where(Th==999.)] = float('nan')
    temperature_celsius = Th*((p/100000)**(0.4/1.4))-273.15
    # Density
    rhodref = mnhFile.variables['RHOREFZ'][:]
    density_3D = np.ones(temperature_celsius.shape)
    IKE = temperature_celsius.shape[0]
    for k in range(IKE):    
        density_3D[k,:,:]=rhodref[k]
    return p, temperature_celsius, density_3D



def get_contents(mnhFile:Dataset,
                 hydrometeors:list,
                 name_var_hydro:dict,
                 temperature:np.ndarray,
                 rho3D:np.ndarray,
                 i_min:int, i_max:int,
                 j_min:int, j_max:int,
                 )->tuple[dict[str,np.ndarray],np.ndarray]:
    contents = {}
    for key in hydrometeors:
        contents[key] = np.empty(temperature.shape)
        contents[key] = mnhFile.variables[name_var_hydro[key]][0,:,:,:][:,i_min:i_max,j_min:j_max]*rho3D[:,:,:] # kg/kg of dry air
        contents[key][contents[key]==999.] = float('nan')
    # water vapour mixing ratio (kg/kg), used for the gas attenuation (it was read but not kept: qv was undefined)
    qv = np.array(mnhFile.variables[name_var_hydro['vv']][0,:,:,:][:,i_min:i_max,j_min:j_max], dtype=float)
    qv[qv==999.] = float('nan')
    return contents,qv



# Number concentration variables of the 2-moment species (current name, former name)
MESONH_CONCENTRATIONS = {'rr': ('CRAIN', 'CRAINT'), 'ii': ('CICE', 'CICET'), 'cc': ('CCLOUD', 'CCLOUDT')}
# Pristine ice concentration of the 1-moment schemes when it is not in the file (as for AROME)
DEFAULT_ICE_CONCENTRATION = 800.


def read_concentration(mnhFile:Dataset, names:tuple, rho3D:np.ndarray,
                       i_min:int, i_max:int, j_min:int, j_max:int) -> np.ndarray|None:
    """Number concentration (m-3) of the first variable of `names` found in the file, None if none is found.
    The variables in kg-1 (LIMA prognostic concentrations) are multiplied by the dry air density; the
    variables in m-3 (CIT of the 1-moment schemes) are kept."""
    for name in names:
        if name in mnhFile.variables:
            var = mnhFile.variables[name]
            field = np.array(var[0,:,:,:][:,i_min:i_max,j_min:j_max], dtype=float)
            field[field==999.] = float('nan')
            units = getattr(var, 'units', '').replace(' ', '')
            if 'kg-1' in units or 'kg^-1' in units or '/kg' in units:
                field *= rho3D
            elif units not in ('m-3', 'm^-3', '/m3'):
                print(f"\tWarning: unknown units '{units}' for {name}: taken as m-3")
            return field
    return None


def get_concentrations(mnhFile:Dataset,
                       hydrometeorMoments:dict[str,int],
                       hydrometeors:list,
                       temperature:np.ndarray,
                       rho3D:np.ndarray,
                       i_min:int, i_max:int,
                       j_min:int, j_max:int,):
    """Number concentrations (m-3), depending on the number of moments of each species (and not on the name
    of the scheme: LIMC is 2-moment for cloud water and rain only, its ice is 1-moment as in ICE3/ICJW):
    * 2-moment species: prognostic concentration (CRAIN, CICE, CCLOUD), in kg-1 -> times the dry air density
    * 1-moment pristine ice: diagnosed concentration CIT (m-3), or 800 m-3 if not in the file (as for AROME)
    * other 1-moment species: not used (0)."""
    concentrations = {key: np.zeros(temperature.shape) for key in hydrometeors}
    for h in hydrometeors:
        if hydrometeorMoments.get(h, 1) == 2:
            names = MESONH_CONCENTRATIONS.get(h)
            field = None if names is None else read_concentration(mnhFile, names, rho3D, i_min, i_max, j_min, j_max)
            if field is None:
                raise KeyError(f"{h} is 2-moment but its concentration ({names}) is not in the file")
            concentrations[h] = field
        elif h == 'ii':
            field = read_concentration(mnhFile, ('CIT',), rho3D, i_min, i_max, j_min, j_max)
            if field is None:
                print(f"\tCIT not in the file: pristine ice concentration = {DEFAULT_ICE_CONCENTRATION:g} m-3")
                field = np.full(temperature.shape, DEFAULT_ICE_CONCENTRATION)
            concentrations[h] = field
    return concentrations
