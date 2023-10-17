import numpy as np
import epygram


import operad_conf as cf

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


# ========== Hydrometeor contents and temperature =============================
"""
   Hydrometeor contents and temperature
   get_contents_t
   input: ficsubdo, p
   output: M, T, R 
"""
def get_contents_t(ficsubdo, p):
    list_t_full=['vv','cc','rr','ii','ss','gg','hh']
    list_hydro=['HUMI.SPECIFI','CLOUD_WATER','RAIN','ICE_CRYSTAL','SNOW','GRAUPEL','HAIL']
    q_cur={}
    q={}
    name_hydro={}
    

    # Arrays initialisation
    for it,t in enumerate(list_t_full):
        name_hydro[t]=list_hydro[it]
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
    
        for t in cf.list_types:
            q_cur[t]=ficsubdo.readfield('S'+'{0:03d}'.format(k+1)+name_hydro[t])
            q[t][k,:,:] = q_cur[t].getdata()
            del q_cur[t]
        
    # Calcul de la "constante" des gaz parfaits du mélange air sec/vapeur 
    # Constante des gaz parfait pour l'air sec
    Rd = epygram.profiles.Rd
    
    # Constante des gaz parfait pour la vapeur d'eau
    Rv = epygram.profiles.Rv        
    R = Rd + q["vv"]*(Rv-Rd) - Rd*(q["cc"]+q["ii"]+q["rr"]+q["ss"]+q["gg"])
    
    M={}
    for t in cf.list_types:
        M[t]=q[t]*p/(R*T)

    return M,T,R
# ===============================================================================




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