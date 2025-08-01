#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import epygram
import numpy as np
from operadar.utils.make_links import link_varname_with_arome_name


def get_2D_lat_lon_epygram(epygram_file) -> tuple[np.ndarray,np.ndarray]:
    """Extract the 2D latitude and longitude field of an epygram file."""
    surface_pressure = epygram_file.readfield('SURFPRESSION')
    surface_pressure.sp2gp() # spectral to grid points
    (lon,lat) = surface_pressure.geometry.get_lonlat_grid()
    return lon, lat



def get_geometry(epygram_file,
                 hybrid_pressure_coefA:np.ndarray,
                 hybrid_pressure_coefB:np.ndarray,
                 )-> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """Retrieve the 3D pressure field (p), the surface pressure (psurf), the 2D geopotential at
    the surface (geosurf) and the pressure departure (pdep).
    
    NOTE : Pressure departure = difference at z_level between pressure and hydrostatic state.
    """
    # Horizontal coordinates 
    surface_pressure = epygram_file.readfield('SURFPRESSION')
    surface_pressure.sp2gp() # spectral to grid points
    surface_pressure = np.exp(surface_pressure.getdata())
    
    # 3D Pressure (Pa)
    pressure = epygram.profiles.hybridP2masspressure(A=hybrid_pressure_coefA,
                                                     B=hybrid_pressure_coefB,
                                                     Psurf=surface_pressure,
                                                     vertical_mean='geometric')

    # 2D Geopotential at surface
    surface_geopotential=epygram_file.readfield('SPECSURFGEOPOTEN')
    surface_geopotential.sp2gp()
    
    # Pressure departure
    pressure_departure = np.zeros(pressure.shape)   
    field_all_levels = epygram_file.readfields('S0*PRESS.DEPART')
    for level,field in enumerate(field_all_levels):
        field.sp2gp()
        pressure_departure[level,:,:] = field.getdata()
    del field_all_levels
    
    # Orography
    #oro = surface_geopotential.getdata(subzone='C')/epygram.profiles.g0

    return pressure, surface_pressure, pressure_departure, surface_geopotential



def get_contents_T_and_R(epygram_file,
                         pressure:np.ndarray,
                         hydrometeors: list,
                         )-> tuple[dict[str,np.ndarray],np.ndarray,np.ndarray]:
    """Retrieve the 3D temperature field, the 3D content fields and compute the gas constant."""
    # 3D kelvin temperature T
    T = np.zeros(pressure.shape)
    temperature_all_levels = epygram_file.readfields('S0*TEMPERATURE')
    for level,field in enumerate(temperature_all_levels):
        field.sp2gp()
        T[level,:,:] = field.getdata()
    del temperature_all_levels 

    # 3D specific content (q in kg/total kg of air)
    name_hydro = link_varname_with_arome_name()
    q = {}
    for key in hydrometeors+['vv'] :
        q[key]=np.zeros(pressure.shape)
        field_all_levels = epygram_file.readfields('S0[0-9][0-9]'+name_hydro[key])
        for level,field in enumerate(field_all_levels):
            q[key][level,:,:] = field.getdata()
        del field_all_levels
    
    # Gas constant computation for the mix dry air/vapour 
    Rd = epygram.profiles.Rd # dry air
    Rv = epygram.profiles.Rv # water vapour
    R = Rd + q["vv"]*(Rv-Rd) - Rd*np.sum(q[x] for x in hydrometeors)
    
    # Transformation of specific content to contents (M in kg/m3)
    M  = {h:q[h]*pressure/(R*T) for h in hydrometeors}

    return M,T,R



def get_concentrations(epygram_file,
                       hydrometeorsConfig: dict,
                       content:dict,
                       temperature:np.ndarray,
                       )-> dict[str,np.ndarray]:
    """Retrieve concentration fields for all hydrometeor classes, depending on their moments."""
    Nc={}
    
    arome_hydrometeors=link_varname_with_arome_name()
    for key,val in arome_hydrometeors.items():
        arome_hydrometeors[key]=val.split('_')[0]
    
    for hydrometeor,moment in hydrometeorsConfig.items():
        Nc[hydrometeor]=np.zeros(temperature.shape)
        if moment == 2:
            extract_N = epygram_file.readfields('S0*N_'+arome_hydrometeors[hydrometeor])
            for level in range(len(extract_N)):
                Nc[hydrometeor][level,:,:] = extract_N[level].getdata()
        elif moment == 1 :
            if hydrometeor == 'ii' :
                Nc["ii"]=800*np.ones(temperature.shape)
            # TODO : explicitely compute Nc for all one moments species   
            # elif hydrometeor == 'ii' :
            #     Nc["ii"]=800*np.ones(temperature.shape)        
    return Nc



def get_altitude(hybrid_pressure_coefA,
                 hybrid_pressure_coefB,
                 temperature,
                 pressure_departure,
                 surface_pressure,
                 surface_geopotential,
                 specific_gas_constant,
                 )-> np.ndarray:
    """Computes the altitude of mass levels defined by hybrid-pressure coefficients."""
    altitude3D = epygram.profiles.hybridP2altitude(A=hybrid_pressure_coefA,
                                                   B=hybrid_pressure_coefB,
                                                   R=specific_gas_constant,
                                                   T=temperature,
                                                   Psurf=surface_pressure,
                                                   vertical_mean='geometric',
                                                   Pdep=pressure_departure,
                                                   Phi_surf=surface_geopotential.getdata(),
                                                   Ptop=np.zeros(surface_pressure.shape)
                                                   )
    return altitude3D