#!/usr/bin/env python
# -*- coding: utf-8 -*-

#===========================================================================
# C. Augros avril 2015
# Adaptation a python de la fonction CALC_KTMAT du simulateur radar
# mode_readtmat.f90 dans MesoNH

# C. Augros 27/02/2020
# Function CALC_KTMAT : add argument NMOMENTS: (1 moment P3=Fw, 2 moments: P3=concentration) 

# C. Augros 12/05/2020
# removal old CALC_KTMAT

# C. Augros 06/2020
# CALC_KTMAT:
# If LAM, ELEV, Tc, P3 or M are outside min and max ranges:
# => warning and the values are set to the min (if below min) or max (if over max)
# Addition Z2dBZ

#===========================================================================

import numpy as np
import math


# ========== Define P3 : CC (2 moments) or Fw (1 moment) ===================
def defineP3(t,moments,mask_tot,Nc,expCCmin, expCCmax, expCCstep,Fw_temp,Fwmin, Fwmax, Fwstep):
    """ Compute 3d parameter in Tmatrix table 
    @todo : change the reading of Tmatrix table so that the parameters are kept the same !

    * pristine ice (treated as 2-moment => the concentration can vary in MesoNH and is registered)
    * rainwater : 2-moment for all LIMA versions / 1-moment for ICE3
    * all other species (graupel / snow / wet graupel / wet snow ? )=> treated as 1-moment but can have a variable wet fraction
    => 3d parameter = Fw
    """    
    if((t =='rr' and moments==2) or (t=='ii')):
        P3name="Nc"
        Nc_temp=Nc[mask_tot]
        P3=np.copy(Nc_temp)
        P3min,P3max,P3step=expCCmin, expCCmax, expCCstep                    
    else:
        P3name="Fw"
        P3=np.copy(Fw_temp)
        P3min,P3max,P3step=Fwmin, Fwmax, Fwstep
        
    return P3, P3min, P3max, P3step, P3name

# ========== Compute NMOMENTS ===================
def compute_nmoments(micro,t):
    if((t =='ii') or ((micro =="LIMA" or micro =="LIMT" or micro =="LIMASG" or micro =="LIMAAG") and (t =='rr'))):
        NMOMENTS=2             
    else:
        NMOMENTS=1
    return NMOMENTS


# ========== Compute single type mask ===================
"""
Compute mask for each type and apply to el, Tc, M, Fw
 * in: M, el, Tc, Fw, mask_precip_dist, expMmin, micro, t
 * out: mask_tot, M_temp,el_temp, Tc_temp, Fw_temp, P3,
 P3min,P3max,P3step
 [mask_tot, M_temp, el_temp, Tc_temp, Fw_temp, P3, P3min, P3max, P3step]=\
     ope_lib.singletype_mask(M, el, Tc, Fw,mask_precip_dist, expMmin,micro)
"""
def singletype_mask(Mt, el, Tc, Fw, mask_precip_dist, expMmin,micro,t):

    # mask_M : selection of precip points only
    mask_M=(Mt>10**expMmin)
    mask_tot=(mask_precip_dist & mask_M)
    el_temp=el[mask_tot]
    Tc_temp=Tc[mask_tot]
    Fw_temp=Fw[mask_tot]
    M_temp=Mt[mask_tot]
    
    # if((t =='ii') or ((micro =="LIMA" or micro =="LIMT" or micro =="LIMASG" or micro =="LIMAAG") and (t =='rr'))):
    #     NMOMENTS=2
    #     if (t =='rr'):
    #         CC_temp=CC[mask_tot]
    #     else: # t='ii'
    #         CC_temp=CCI[mask_tot]
    #     P3=np.copy(CC_temp)
    #     P3min,P3max,P3step=expCCmin, expCCmax, expCCstep                    
    # else:
    #     NMOMENTS=1
    #     P3=np.copy(Fw_temp)
    #     P3min,P3max,P3step=Fwmin[t], Fwmax[t], Fwstep[t]
        
    return mask_tot, M_temp, el_temp, Tc_temp, Fw_temp #, P3, P3min, P3max, P3step

# ============== Compute radar geometry ==========
"""
Compute radar geometry variables in model grid

input: * X, Y, Z: model coordinate vectors + altitude
       * X0, Y0: radar location in model grid 
       * Z0: radar altitude
       * distmax_rad: max range (km) from radar location (no need to compute simulated data after this range)

output: * distmat_mod: 3D matrix in model geometry with distance to radar
        * mask_distmax

"""
def compute_radargeo(X,Y,Z,X0,Y0,Z0,distmax_rad,RT,elevmax):

    XX,YY=np.meshgrid(X-X0,Y-Y0)
        
    radardist=(XX**2+YY**2)**0.5 #2D array with radar distance
                                               # for each model point
    radardist3D=np.ones(Z.shape)
    for level in range(Z.shape[0]):
        radardist3D[level,:,:]=radardist
    
    mask_distmax=(radardist3D<distmax_rad)
    
    # radar elevation
    tanel = Z/radardist3D-3.*radardist3D/(8.*RT)
    el = np.arctan(tanel)*180./math.pi
    el[el<0] = 0.
    el[el>elevmax] = elevmax
    el[:,:,:]=0.
    
    return mask_distmax,el

# ============== Compute precipitation + maxdistance mask ==========
def mask_precip(mask_distmax,M,expMmin,micro):
    Mtot=np.copy(M['rr'])
    if(micro =="LIMAAG" or micro =="ICE4"):
        Mtot=M['rr']+M['gg']+M['ss']+M['ii']+M['hh']	
    else :											
        Mtot=M['rr']+M['gg']+M['ss']+M['ii']							
    mask_precip=(Mtot>10**expMmin)  
    
    # Precip + distance to radar mask
    mask_precip_dist=(mask_precip & mask_distmax)
    
    return mask_precip_dist    


    
# ============== Compute water fraction for wet species ==========
"""
Add wet types in hydrometeor contents (Mwg, Mwh)
and compute water fraction for wet species

* input: M, mixed phase option (Fwpos, Tpos, Fwposg)
 *output: Fw 3D, M with addition of wet hydrometeor types wg, wh 
"""
def compute_mixedphase(M, Nc, MixedPhase,expMmin,micro):
    
    # Bright band (or mixed phase) mask 
   # Mtot=np.copy(M['rr'])
    if(micro =="LIMAAG" or micro =="ICE4"):
        #Mtot=M['rr']+M['gg']+M['ss']+M['ii']+M['hh']	
        maskBB=((M["rr"] > 10**expMmin) & ((M["gg"]> 10**expMmin) | (M["hh"]> 10**expMmin)))
    else :											
        #Mtot=M['rr']+M['gg']+M['ss']+M['ii']							
        maskBB=((M["rr"] > 10**expMmin) & (M["gg"]> 10**expMmin))
    
    print("Calculation of Fw for wet graupel")
    Fw = np.zeros(np.shape(M["rr"]))                          
 
    if(micro =="LIMAAG" or micro =="ICE4"):
         Fw[maskBB] = (M["rr"]/(M["rr"]+M["gg"]+M["hh"]))[maskBB]
         M["wh"] = np.copy(M["hh"])
         Nc["wh"] = np.copy(Nc["hh"])
    else :
         Fw[maskBB] = (M["rr"]/(M["rr"]+M["gg"]))[maskBB]
    
    M["wg"] = np.copy(M["gg"])
    Nc["wg"] = np.copy(Nc["gg"])
    
#    # MixedPhase=="Tpos" => Graupel is transferred to melting graupel if T>=0   # => idem for hail           
#    if (MixedPhase=="Tpos"):  
#        M["wg"][Tc < 0] = 0
#        M["gg"][Tc >= 0] = 0
#        M["wh"][Tc < 0] = 0
#        M["hh"][Tc >= 0] = 0

    # MixedPhase=="Fwpos" => Graupel is transferred to melting graupel if Fw>=0     
    if (MixedPhase=="Fwpos"):
        M["wg"][Fw == 0] = 0
        M["wg"][maskBB] = M["gg"][maskBB]+M["rr"][maskBB] # If M["rr"] > 10**expMmin) & (M["gg"]> 10**expMmin)
                                   # the rainwater is added to the wet graupel content         
        M["rr"][maskBB] = 0        # and removed from the rain content  
        M["gg"][maskBB] = 0
       
        if(micro =="LIMAAG" or micro =="ICE4"):
            M["wh"] = M["hh"] + ( (M["rr"]*M["hh"])/(M["hh"]+M["gg"]) )	# d'après Wolfensberger, 2018

    if (MixedPhase=="Fwposg"):
        M["wg"][Fw == 0] = 0
        M["wg"][maskBB] = M["gg"][maskBB]
        M["gg"][maskBB] = 0
       
        if(micro =="LIMAAG" or micro =="ICE4"):
            M["wh"][Fw == 0] = 0
            M["wh"][maskBB] = M["hh"][maskBB]
            M["hh"][maskBB] = 0 # addition Clotilde (in the mixed phase, there is no dry hail)
            
    return M, Nc, Fw


#============== Interpolation linéaire de y1, y2 en x1,x2 ==============
def lin_interpol(x1,x2,y1,y2,x):
    if (x1==x2):
        res=0.5*(y1+y2)
    else:
        res=(y1*(x2-x)+y2*(x-x1))/(x2-x1)
    return res
    

#=============== Linear to dBZ conversion =====   
"""
Linear to dBZ conversion
""" 
def Z2dBZ(Z):
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Z[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp
#=============================


def latlon2XY(latrad,lonrad,latitude,longitude):
    """
    Fonction de conversion lat lon vers indices x y dans la grille modèle
    dx-y = resolution
    
    """
    dx = (longitude[2]-longitude[1])/2.
    dy = (latitude[2]-latitude[1])/2.
    
    lat0 = np.where((latrad-dy <latitude) & (latitude< latrad+dy))
    lon0 = np.where((lonrad-dx <longitude) & (longitude< lonrad+dx))
    jlat=lat0[0][0]
    ilon = lon0[0][0]
    
    if np.size(lat0)>1:
        if abs(latitude[lat0[0][0]]-latrad)<abs(latitude[lat0[0][1]]-latrad):
            jlat=lat0[0][0]
        else:
            jlat=lat0[0][1]
            
    if np.size(lon0)>1:
        if abs(longitude[lon0[0][0]]-lonrad)<abs(longitude[lon0[0][1]]-lonrad):
            ilon=lon0[0][0]
        else:
            ilon=lon0[0][1]
            
    return jlat, ilon


def compute_grid_alt(var3D,ztop,orography,level):
    altitude = np.zeros_like(var3D[0,:,:,:])
    nx = len(var3D[0,0,0,:])
    ny = len(var3D[0,0,:,0])
    for ii in range(0,ny):
        for jj in range(0,nx):
            if ZS[ii, jj] == 0:
                altitude[:, ii,jj] = level[:]
            else:
                for k in range(len(level)):
                    altitude[k, ii,jj] = orography[ii,jj] + level[k] * ((ztop - orography[ii,jj]) / ztop)
    return altitude


#======================================
#Opening net_cdf file
#======================================
def open_nc(file_name,variable,dim):
    from netCDF4 import Dataset
    
    fic=Dataset(file_name,'r')
    var=[]
    
    if dim<2 or dim>4:
        print('Error dim should be 3 or 4')
        var=99999 ; print('var=',var)
    if dim==2: 
        var=fic.variables[variable][:,:] # for lat and lon
        print(var.shape)
    if dim==3: 
        var=fic.variables[variable][:,1:-1,1:-1]
        print(var.shape)
    if dim==4 : 
        var=fic.variables[variable][:,:,1:-1,1:-1]  
        print(var.shape)
        
    return var

#==============================================================================
# Fonction get indice
#==============================================================================
def getlevel(zr,vec,eps):
    zz0 = np.where((zr-eps <vec)&(vec< zr+eps))
    zz=zz0[0][0]
    if np.size(zz0)>1:
        if abs(vec[zz0[0][0]]-zr)<abs(vec[zz0[0][1]]-zr):
            zz=zz0[0][0]
        else:
            zz=zz0[0][1]
    return zz


#============== Fonction Cropping lat lon ==============

def crop_latlon(file,latmin,latmax,lonmin,lonmax):
                 
    LAT = open_nc(file,'latitude',2) ;  #print(np.min(LAT[:,0]),np.max(LAT[:,0]))
    LON = open_nc(file,'longitude',2) ; #print(np.min(LON[0,:]),np.max(LON[0,:]))

    dlat = abs(LAT[2,0]-LAT[1,0]) ; dlon = abs(LON[0,2]-LON[0,1]) ; #print(dlat,dlon)

    ii_latmin=getlevel(latmin,LAT[:,0],dlat) ; #print(ii_latmin,LAT[ii_latmin,0])
    ii_latmax=getlevel(latmax,LAT[:,0],dlat) ; #print(ii_latmax,LAT[ii_latmax,0])

    ii_lonmax= getlevel(lonmax,LON[0,:],dlon) ; #print(ii_lonmax,LON[0,ii_lonmax])
    ii_lonmin = getlevel(lonmin,LON[0,:],dlon) ; #print(ii_lonmin,LON[0,ii_lonmin])
    
    return ii_latmin,ii_latmax,ii_lonmin,ii_lonmax





##"""======================================================================
#                                      #Test fonctions
# #======================================================================"""
#
#if __name__ == '__main__':
#
#    LAMmin,LAMstep,LAMmax=53.2,0.1,53.2
#    ELEVmin,ELEVmax,ELEVstep=0.0,4.0,12.0
#    Tcmin,Tcstep,Tcmax=-20.0,1.0,40.0
#    Fwmin,Fwstep,Fwmax=0.0,0.1,0.0
#
#    nk=2
#    nj=2
#    ni=2
#
#    LAMm=0.0532
#    el=1.5
#    ELEVrad=np.full(shape=(nk,nj,ni),fill_value=el/180*math.pi)
#    Tc=np.full(shape=(nk,nj,ni),fill_value=20.3)
#    Fw=np.full(shape=(nk,nj,ni),fill_value=0)
#    M=np.full(shape=(nk,nj,ni),fill_value=10.1**(-5))


