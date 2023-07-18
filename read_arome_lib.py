import numpy as np
import epygram


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


def get_geometry(ficsubdo, A, B, IKE):
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
    for k in range(IKE): #going downward
        pdep_cur=ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+'PRESS.DEPART')
        pdep_cur.sp2gp()
        pdep[k,:,:] = pdep_cur.getdata()
        del pdep_cur

    # Orography
    #oro = phis.getdata(subzone='C')/epygram.profiles.g0

    return p, psurf, pdep, phis

# ==============================================================================

def link_varname_with_realname ():
    list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    list_hydro=['HUMI.SPECIFI','CLOUD_WATER','RAIN','ICE_CRYSTAL','SNOW','GRAUPEL','HAIL']
    name_hydro_linked={t:list_hydro[it] for it,t in enumerate(list_t_full)}
    return name_hydro_linked, list_t_full


# ========== Hydrometeor contents and temperature =============================
"""
   Hydrometeor contents and temperature
   get_contents_and_T
   input: ficsubdo, p
   output: M, T, R 
"""
def get_contents_and_T(ficsubdo, p, hydrometeor_type: list):
    
    name_hydro , list_t_full = link_varname_with_realname()

    # Arrays initialisation
    q = {t:np.zeros(p.shape) for t in list_t_full}
    T = np.zeros(p.shape)
    
    # 3D specific content q and temperature T 
    IKE=p.shape[0]
    for k in range(IKE): #going downward
        # Temperature
        T_cur = ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+'TEMPERATURE')
        T_cur.sp2gp()
        T[k,:,:] = T_cur.getdata()
        del T_cur    
        # Hydrometeors
        q_cur={}
        for htype in hydrometeor_type :
            q_cur[htype]=ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+name_hydro[htype])
            q[htype][k,:,:] = q_cur[htype].getdata()
            del q_cur[htype]
        
    # Calcul de la "constante" des gaz parfaits du mélange air sec/vapeur 
    Rd = epygram.profiles.Rd # air sec
    Rv = epygram.profiles.Rv # vapeur d'eau
    R  = Rd + q["vv"]*(Rv-Rd) - Rd*(q["cc"]+q["ii"]+q["rr"]+q["ss"]+q["gg"])
    
    M  = {t:q[t]*p/(R*T) for t in hydrometeor_type}

    return M,T,R


# ========== Hydrometeor concentrations =============================
def get_concentrations(ficsubdo, p, microphysics: str, hydrometeor_type: list, iceCC_cst:float = 800.):

    cc_rain = np.empty(p.shape)
    cc_ice  = iceCC_cst*np.ones(p.shape)
               
    if microphysics[0:4] == "LIMA":

        for k in range(p.shape[0]):
            extract_rr_cc = ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+'N_RAIN')
            cc_rain[k,:,:] = extract_rr_cc.getdata()
            extract_ii_cc = ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+'N_ICE')
            cc_ice[k,:,:] = extract_ii_cc.getdata()
         
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
#===============================================================================