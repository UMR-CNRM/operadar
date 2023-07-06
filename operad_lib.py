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


import operad_conf as cf

# ===========================================================================



# ========== Define P3 ===================
"""
Define P3 : CC (2 moments) or Fw (1 moment)
"""
def defineP3(t,NMOMENTS,CC,CCI,mask_tot,Fw_temp,expCCmin, expCCmax, expCCstep,Fwmin, Fwmax, Fwstep):
    if(NMOMENTS==2):
        if (t =='rr'):
            CC_temp=CC[mask_tot]
        else: # t='ii'
            CC_temp=CCI[mask_tot]
        P3=np.copy(CC_temp)
        P3min,P3max,P3step=expCCmin, expCCmax, expCCstep                    
    elif (NMOMENTS==1):
        P3=np.copy(Fw_temp)
        P3min,P3max,P3step=Fwmin, Fwmax, Fwstep
        
    return P3, P3min, P3max, P3step

# ========== Compute NMOMENTS ===================
"""
Compute nmoments
"""
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
    
    
    # if((t =='ii') or ((cf.micro =="LIMA" or cf.micro =="LIMT" or cf.micro =="LIMASG" or cf.micro =="LIMAAG") and (t =='rr'))):
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

#============== Compute radar geometry ==========
"""
Compute radar geometry variables in model grid

input: * X, Y, Z: model coordinate vectors 
       * X0, Y0: radar location in model grid 
       * distmax_rad: max range (km) from radar location (no need to compute simulated data after this range)

output: * distmat_mod: 3D matrix in model geometry with distance to radar
        * mask_distmax

"""
def compute_radargeo(X,Y,Z,X0,Y0,distmax_rad,RT,elevmax):

    XX,YY=np.meshgrid(X-X0,Y-Y0)
        
    # distmat_mod : ground distance to the radar location for each model point
    distmat_mod=np.array([(XX**2+YY**2)**0.5]*len(Z))
    Z=np.array([np.ones(np.shape(XX))*Z[k] for k in range(len(Z))]) 
    mask_distmax=(distmat_mod<cf.distmax_rad)
    
    # radar elevation
    tanel = Z/distmat_mod-3.*distmat_mod/(8.*RT)
    el = np.arctan(tanel)*180./math.pi
    el[el<0] = 0.
    el[el>elevmax] = elevmax
    
    return mask_distmax,el
# =================================================

#==================================================
"""
Compute precipitation + maxdistance mask
"""
def mask_precip(mask_distmax,M,expMmin):
    Mtot=np.copy(M['rr'])
    if(cf.micro =="LIMAAG" or cf.micro =="ICE4"):
        Mtot=M['rr']+M['gg']+M['ss']+M['ii']+M['hh']	
    else :											
        Mtot=M['rr']+M['gg']+M['ss']+M['ii']							
    mask_precip=(Mtot>10**expMmin)  
    
    # Precip + distance to radar mask
    mask_precip_dist=(mask_precip & mask_distmax)
    
    return mask_precip_dist
#===================================================    



    
#=====================================================================
"""
Add wet types in hydrometeor contents (Mwg, Mwh)
and compute water fraction for wet species

* input: M, mixed phase option (Fwpos, Tpos, Fwposg)
 *output: Fw 3D, M with addition of wet hydrometeor types wg, wh 
"""
def compute_mixedphase(M,MixedPhase,expMmin):
    
    # Bright band (or mixed phase) mask 
   # Mtot=np.copy(M['rr'])
    if(cf.micro =="LIMAAG" or cf.micro =="ICE4"):
        #Mtot=M['rr']+M['gg']+M['ss']+M['ii']+M['hh']	
        maskBB=((M["rr"] > 10**expMmin) & ((M["gg"]> 10**expMmin) | (M["hh"]> 10**expMmin)))
    else :											
        #Mtot=M['rr']+M['gg']+M['ss']+M['ii']							
        maskBB=((M["rr"] > 10**expMmin) & (M["gg"]> 10**expMmin))
    
    print("Calculation of Fw for wet graupel")
    Fw = np.zeros(np.shape(M["rr"]))                          
 
    if(cf.micro =="LIMAAG" or cf.micro =="ICE4"):
         Fw[maskBB] = (M["rr"]/(M["rr"]+M["gg"]+M["hh"]))[maskBB]
         M["wh"] = np.copy(M["hh"])
    else :
         Fw[maskBB] = (M["rr"]/(M["rr"]+M["gg"]))[maskBB]
    
    M["wg"] = np.copy(M["gg"])
    
    # cf.MixedPhase=="Tpos" => Graupel is transferred to melting graupel if T>=0   # => idem for hail           
    if (cf.MixedPhase=="Tpos"):  
        M["wg"][Tc < 0] = 0
        M["gg"][Tc >= 0] = 0
        M["wh"][Tc < 0] = 0
        M["hh"][Tc >= 0] = 0

    # cf.MixedPhase=="Fwpos" => Graupel is transferred to melting graupel if Fw>=0     
    if (cf.MixedPhase=="Fwpos"):
        M["wg"][Fw == 0] = 0
        M["wg"][maskBB] = M["gg"][maskBB]+M["rr"][maskBB] # If M["rr"] > 10**expMmin) & (M["gg"]> 10**expMmin)
                                   # the rainwater is added to the wet graupel content         
        M["rr"][maskBB] = 0        # and removed from the rain content  
        M["gg"][maskBB] = 0
       
        if(cf.micro =="LIMAAG" or cf.micro =="ICE4"):
            M["wh"] = M["hh"] + ( (M["rr"]*M["hh"])/(M["hh"]+M["gg"]) )	# d'après Wolfensberger, 2018

    if (cf.MixedPhase=="Fwposg"):
        M["wg"][Fw == 0] = 0
        M["wg"][maskBB] = M["gg"][maskBB]
        M["gg"][maskBB] = 0
       
        if(cf.micro =="LIMAAG" or cf.micro =="ICE4"):
            M["wh"][Fw == 0] = 0
            M["wh"][maskBB] = M["hh"][maskBB]
            M["hh"][maskBB] = 0 # addition Clotilde (in the mixed phase, there is no dry hail)
            
    return M, Fw
#==========================================================================


#============== Fonction lin_interpol ==============
"""
Interpolation linéaire de y1, y2 en x1,x2
"""
def lin_interpol(x1,x2,y1,y2,x):
    if (x1==x2):
        res=0.5*(y1+y2)
    else:
        res=(y1*(x2-x)+y2*(x-x1))/(x2-x1)
    return res
#=============================
    

    
def Z2dBZ(Z):
    Ztemp = np.copy(Z)
    Ztemp[Z > 0.] = 10.*np.log10(Z[Z > 0.])
    Ztemp[Z == 0.] = -999.
    Ztemp[Z < 0.] = -9999.
    return Ztemp



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


