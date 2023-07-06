import numpy as np
import epygram


# ======== Horizontal, vertical coordinates, pressure ===========================
"""
Horizontal, vertical coordinates, pressure 
input: ficA
output: p, A, B, nkA, lon, lat 
"""
def get_geometry(ficA,ficsubdo):
    # === Horizontal coordinates 
    ps = ficsubdo.readfield('SURFPRESSION')
    ps.sp2gp() # spectral to grid points
    psurf = np.exp(ps.getdata())

    (lon,lat) = ps.geometry.get_lonlat_grid(subzone='C')


    # Vertical levels values
    A = [level[1]['Ai'] for level in ficA.geometry.vcoordinate.grid['gridlevels']][1:]
    B = [level[1]['Bi'] for level in ficA.geometry.vcoordinate.grid['gridlevels']][1:]
        
    # Number of vertical levels
    IKE=len(ficA.geometry.vcoordinate.levels)

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

    return p, psurf, pdep, phis, A, B, lon, lat 

# ==============================================================================

def link_varname_with_realname ():
    list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    list_hydro=['HUMI.SPECIFI','CLOUD_WATER','RAIN','ICE_CRYSTAL','SNOW','GRAUPEL','HAIL']
    
    name_hydro_linked={}
    for it,t in enumerate(list_t_full):
        name_hydro_linked[t]=list_hydro[it]
        
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
    q={}
    for it,t in enumerate(list_t_full):
        q[t]=np.zeros(p.shape)
        
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
    # Constante des gaz parfait pour l'air sec
    Rd = epygram.profiles.Rd
    
    # Constante des gaz parfait pour la vapeur d'eau
    Rv = epygram.profiles.Rv        
    R = Rd + q["vv"]*(Rv-Rd) - Rd*(q["cc"]+q["ii"]+q["rr"]+q["ss"]+q["gg"])
    
    M={}
    for t in hydrometeor_type:
        M[t]=q[t]*p/(R*T)

    return M,T,R
# ===============================================================================

def get_concentrations(ficsubdo, p, Tc, microphysics: str, hydrometeor_type: list, iceCC_cst:float = 800):
    
    cc_rain = np.empty(Tc.shape)
    cc_ice = iceCC_cst*np.ones(Tc.shape)
               
    if microphysics[0:4] == "LIMA":
        name_hydro , list_t_full = link_varname_with_realname()
        """
        # Arrays initialisation
        q={}
        for it,t in enumerate(list_t_full):
            q[t]=np.zeros(p.shape)
        IKE=p.shape[0]
        for k in range(IKE):
            q_cur={}
            for htype in hydrometeor_type:
                q_cur[htype]=ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+name_hydro[htype])
                q[htype][k,:,:] = q_cur[htype].getdata()
                del q_cur[htype]
        """   
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