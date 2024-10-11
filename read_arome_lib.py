import numpy as np
import epygram
import math
import time as tm
# ======== Horizontal, vertical coordinates, pressure ===========================
"""
Horizontal, vertical coordinates, pressure 
input: ficsubdo
output: p, A, B, nkA
"""

def get_lat_lon_epygram(ficsubdo):
    ps = ficsubdo.readfield('SURFPRESSION')
    ps.sp2gp() # spectral to grid points
    (lon,lat) = ps.geometry.get_lonlat_grid(subzone='C')
    return lon, lat


def get_geometry(ficsubdo, A, B):
    # === Horizontal coordinates 
    ps = ficsubdo.readfield('SURFPRESSION')
    ps.sp2gp() # spectral to grid points
    psurf = np.exp(ps.getdata())
    
    # 3D Pressure (Pa)
    p = epygram.profiles.hybridP2masspressure(A, B, psurf, 'geometric')

    # 2D Geopotential at surface
    phis=ficsubdo.readfield('SPECSURFGEOPOTEN')
    phis.sp2gp()
    
    # Pressure depart: difference at z_level: pressure - hydrostatic state
    pdep = np.zeros(p.shape)   
    field_all_levels = ficsubdo.readfields('S0*PRESS.DEPART')
    for k,field in enumerate(field_all_levels):
        field.sp2gp()
        pdep[k,:,:] = field.getdata()
    del field_all_levels
    
    # Orography
    #oro = phis.getdata(subzone='C')/epygram.profiles.g0

    return p, psurf, pdep, phis

# ==============================================================================

def link_varname_with_realname ():
    list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    list_hydro=['HUMI.SPECIFI','CLOUD_WATER','RAIN','ICE_CRYSTAL','SNOW','GRAUPEL','HAIL']
    name_hydro_linked={t:list_hydro[it] for it,t in enumerate(list_t_full)}
    return name_hydro_linked, list_t_full


# ========== Get hydrometeor contents and temperature =============================
"""
   input: ficsubdo, p, list of hydrometeors type
   output: M, T, R 
"""
def get_contents_and_T(ficsubdo, p, hydrometeor_type: list):
    
    name_hydro , list_t_full = link_varname_with_realname()

    # Arrays initialisation
    q = {t:np.zeros(p.shape) for t in list_t_full}
    T = np.zeros(p.shape)
    
    # 3D temperature T
    temperature_all_levels = ficsubdo.readfields('S0*TEMPERATURE')
    for k,field in enumerate(temperature_all_levels):
        field.sp2gp()
        T[k,:,:] = field.getdata()
    del temperature_all_levels 
    
    # 3D specific content q 
    for htype in hydrometeor_type :
        field_all_levels = ficsubdo.readfields('S0[0-9][0-9]'+name_hydro[htype])
        for k,field in enumerate(field_all_levels):
            q[htype][k,:,:] = field.getdata()
        del field_all_levels
    
    # Calcul de la "constante" des gaz parfaits du mélange air sec/vapeur 
    Rd = epygram.profiles.Rd # air sec
    Rv = epygram.profiles.Rv # vapeur d'eau
    R  = Rd + q["vv"]*(Rv-Rd) - Rd*(q["cc"]+q["ii"]+q["rr"]+q["ss"]+q["gg"])
    
    M  = {t:q[t]*p/(R*T) for t in hydrometeor_type}

    return M,T,R


# ========== Hydrometeor concentrations =============================
def get_concentrations(ficsubdo, p, microphysics: str, hydrometeor_type: list, iceCC_cst:float = 800.,hydrometeor_content=None):

    cc_rain = np.empty(p.shape)
    cc_ice  = iceCC_cst*np.ones(p.shape)
    if microphysics == "ICE3":
        #print(hydrometeor_content["rr"])
        cc_rain = 8e6 * SlopeParameterICE3(M= hydrometeor_content["rr"],
                                           alpha= 1, nu=1, a=524, b=3, C=8e6, X=-1,
                                           ) # A TESTER
    if microphysics[0:4] == "LIMA":
        extract_rr_cc = ficsubdo.readfields('S0*N_RAIN')
        extract_ii_cc = ficsubdo.readfields('S0*N_ICE')
        for k in range(len(extract_rr_cc)):
            cc_rain[k,:,:] = extract_rr_cc[k].getdata()
            cc_ice[k,:,:] = extract_ii_cc[k].getdata()
         
    return cc_rain, cc_ice

# =============== Altitude ======================================================
"""
   Altitude z of each level
   input: A, B, nkA, T, p, R
   output: z [i,j,k]
"""
#def get_altitude(A, B, niA, njA, nkA, T, p, pdep, Psurf, phis, R):
def get_altitude(A, B, T, p, pdep, Psurf, phis, R):
    z = epygram.profiles.hybridP2altitude(A, B, R, T, Psurf, 'geometric', Pdep=pdep, Phi_surf=phis.getdata(), Ptop=np.zeros(Psurf.shape))

    return z
#===========Slope parameter (only for ICE3)==========================================
"""
   INPUTS :
    - M          = hydrometeor content
    - nu & alpha = dispersion parameters
    - a & b      = prognostic variables
    - C          = densite volumique de particules ?
    - X
"""
def SlopeParameterICE3(M,alpha,nu,a,b,C,X):
    lambda_hydrometeor =((M*math.gamma(nu))/(a*C*math.gamma(nu+b/alpha)))**(1.0/(X-b))      
    return lambda_hydrometeor 