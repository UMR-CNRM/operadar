#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018

@author: augrosc & lebastardt & montangonm & sinhorin & davidcl

Dual Polarization Radar Forward operator :
calculation of dpol radar variables for each model grid point
- INPUT : 
    * a model file (MesoNH netcdf or AROME fa): modelfile => contains Z (altitude), Tc (temperature), hydrometeor contents M and concentration CC
    * the scattering coefficients tables (output files of Tmatrix) recorded for
    a range of T, M, CC (if 2 moments) 
- OUTPUT : 
    * ascii or .npz file with modele points i,j,k and dpol var Zh, Zdr, Kdp, Rhohv 
    
The modifications are recorded via git:
    git status (to see the modified and not recorded files)
    git log (for the historic of all modifications)
    git gui (visual tool for git)
"""

#===================================================================================
# C. Augros 12/05/2020
# - selection of Tmat version using repTmat
# - new variables: CCI, mu for LIMT, Dv
# - new functions to calculate Dm and SlopeParameter for LIMT :
# MeanMassDiameterLIMT SlopeparameterLIMT in plot_MesoNH_bib.py
# - 2 options for mixed phase: Tpos or Fwpos


# C. Augros 20/05/2020
# - simplification/reorganization of the code
# - for each hydrometeor type, dpol variables are calculated only for M > Mmin
# - mask for all hydromet (Mtot > Mmin) and for radar distance 
# - removal of the subdomains
# - bug corrections in the mixed phase with a new option: cf.MixedPhase="Fwpos" or "Tpos"
# => higher Zh and Zdr for T<0 with Fwpos !!!
# # - new option for the output file: native python compressed format (.npz):
# faster !!!!!
# - calculation of dpol variables for each hydromet type if option singletype=True
# => variables are saved in separate files

# C. Augros 06/2020
# - 3 options available for the mixed phase:
# * Tpos : wet graupel only for T >= 0
# * Fwpos : wet graupel for Fw > 0 and Mwg=Mg+Mr if Fw >0
# * Fwposg: wet graupel for Fw > 0 and Mwg=Mg if Fw >0
# - new LIMToption="" or "cstmu" => the model variables are taken from LIMT simulation
                # but a constant mu is applied in the PSD  for the dpol variables calculation 

# C. Augros 7/04/2023
# New version with addition of functions (in operad_lib) + configuration file (operad_conf.py)
# Objective = 
# * include in this unique code the MesoNH/AROME options with ICE3/ICE4 or LIMA microphysics
# * improve the lisibility and efficiency
# 
# Next step = simplify also plot_Hcut... ==> gather everything on git
                
# C. Augros 19/04/2023
# Code ok for Arome (ICE3/ICE4?) and MesoNH
# output = netcdf file (dataset with xarray)                
# implemented on github  https://github.com/augros/operadar               
#=====================================================================================



import sys
import os
import math
import numpy as np


# operad modules
import operad_lib as ope_lib
import read_arome as aro
import read_mesonh as meso
import read_tmatrix as read_tmat
import save_dpolvar as save
from pathlib import Path


#============= Parameters to configure =========================

configfile="operad_conf_AROME_ICE4.py"
#configfile="operad_conf_MesoNH_ICE3idpx.py"

os.system("cp "+configfile+" operad_conf.py")


import operad_conf as cf


#============= Programm ====================================

# ----- Test existence of pathfick => if not: creation of directory
if (os.path.exists(cf.pathfick)): 
    print ('pathfick exists : '+cf.pathfick)            
else:
    try:
        os.system("mkdir "+cf.pathfick)
    except:    
        print ('error in creation of '+cf.pathfick)
        sys.exit()

# ----------------------


liste_var_pol = ["Zhh", "Zdr", "Kdp","Rhohv"]
liste_var_calc=["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]



# -----------------------------------------------------------------------------
# ---------------------- Reading Tmatrix tables -------------------------------
#------------------------------------------------------------------------------
print("Reading Tmatrix tables")

if (cf.micro=="LIMT" and cf.LIMToption=="cstmu"):
    [LAMmin,LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax, 
 Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, 
 Fwmax,expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, 
 Tc_t, ELEV_t, Fw_t, M_t, S11carre_t, S22carre_t, ReS22S11_t,
 ImS22S11_t, ReS22fmS11f_t, ImS22ft_t, ImS11ft_t]=read_tmat.Read_TmatrixClotilde(cf.pathTmat,cf.band,"LIMA",cf.table_ind)

else:        
    [LAMmin,LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax, 
     Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, 
     Fwmax,expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, 
     Tc_t, ELEV_t, Fw_t, M_t, S11carre_t, S22carre_t, ReS22S11_t,
     ImS22S11_t, ReS22fmS11f_t, ImS22ft_t, ImS11ft_t]=read_tmat.Read_TmatrixClotilde(cf.pathTmat,cf.band,cf.micro,cf.table_ind)
LAM=LAMmin["rr"]/1000.
print("End reading Tmatrix tables")


# ----------------------------------------------------------------------------
# ----------------------- Loop over timesteps -----------------------------
# ----------------------------------------------------------------------------
for time in cf.timelist: 
    print("-------",time,"-------")
    fick = cf.pathfick+"k_"+cf.model+"_"+ cf.band+'_'+str(int(cf.distmax_rad/1000.))+"_ech"+time+"_2"
    
    if Path(fick+".nc").exists():
        print("operad netcdf file for",time,"already exists")
        continue   
    
    # =========== Reading model variables =============
    print("Reading model variables")    
    
    # Return 3D model variables + coordinates
    if (cf.model=="MesoNH"):
        [M, Tc, CC, CCI, X, Y, Z]=meso.read_mesonh(cf.micro,time)
    
    elif (cf.model=="Arome"):
        [M, Tc, CC, CCI, lon, lat, Z]=aro.read_arome(cf.micro,time)
    else:
        print("cf.model="+cf.model+" => needs to be either Arome or MesoNH")
      
    # =========== Compute radar geometry variables =====
    print("Compute radar geometry: elevation (el) and distance mask (mask_distmax)")
    # TODO: add latlon2XY function in ope_lib
    # else: 
    #     X0,Y0=ope_lib.latlon2XY(cf.latrad,cf.lonrad)
    # ----------------
    if (cf.model=="Arome"):
        el=np.zeros(Tc.shape)
        mask_distmax=(el >= 0.)
    
    elif (cf.model=="MesoNH"):
        if (cf.radarloc=="center"):
            X0=np.nanmean(X)
            Y0=np.nanmean(Y)        
        [mask_distmax,el]=ope_lib.compute_radargeo(X,Y,Z,X0,Y0,cf.distmax_rad,cf.RT,ELEVmax["rr"])

    # ============= Precip mask ========================= 
    mask_precip_dist=ope_lib.mask_precip(mask_distmax,M,expMmin) 


    # ============= Mixed phase parametrization ========    
    [M,Fw]=ope_lib.compute_mixedphase(M,cf.MixedPhase,expMmin)
  
    
    #  ============ Initialization of the dictionnary Vm_k 
    # contains all 3D dpol variables (all hydrometeor included)
    Vm_k={}
    for var in liste_var_calc:
        Vm_k[var]=np.zeros(Tc.shape)  
        
    # ==================== Loop over hydromet types ===========================
    print("Loop over hydrometeor types in Tmatrix tables:",cf.list_types_tot)
    for t in cf.list_types_tot:
        if (cf.singletype):
            Vm_t={}
            for var in liste_var_calc:
                Vm_t[var]=np.zeros(Tc.shape) 
         
        #============== Compute NMOMENTS ========================
        NMOMENTS=ope_lib.compute_nmoments(cf.micro,t)
        
        
        # ========== Compute single type mask ===================
        [mask_tot, M_temp, el_temp, Tc_temp, Fw_temp]=\
             ope_lib.singletype_mask(M[t], el, Tc, Fw,mask_precip_dist, expMmin,cf.micro,t)
        
        # ========= Define P3 : CC (2 moments) or Fw (1 moment) =
        [P3, P3min, P3max, P3step]=ope_lib.defineP3(t,NMOMENTS,CC,CCI,mask_tot,Fw_temp, \
          expCCmin, expCCmax, expCCstep,Fwmin[t], Fwmax[t], Fwstep[t])
        
        # ========= Extract scattering coefficients for singletype
        [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11]=read_tmat.get_scatcoef(\
          S11carre_t[t],S22carre_t[t],ReS22fmS11f_t[t],ReS22S11_t[t],ImS22S11_t[t],\
          LAMmin[t], LAMmax[t], LAMstep[t], ELEVmin[t], ELEVmax[t], ELEVstep[t],\
          Tcmin[t], Tcmax[t], Tcstep[t], P3min, P3max, P3step, expMmin,expMstep,expMmax,\
          NMOMENTS, el_temp,Tc_temp,P3, M_temp)
       
            
        # ========== Single type dpol var computation ===================       
        Vm_temp={}
        for var in liste_var_calc:
            Vm_temp[var]=Vm_k[var][mask_tot]
        Vm_temp["Zhhlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S22carre
        Vm_temp["Zvvlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S11carre
        Vm_temp["Kdp"] = 180.*1e3/math.pi*LAM*ReS22fmS11f
        Vm_temp["S11S22"] = ReS22S11**2+ImS22S11**2
        Vm_temp["S11S11"] = np.copy(S11carre)
        Vm_temp["S22S22"] = np.copy(S22carre)
        
        # =========== Addition of scattering coef for all hydromet ========
        for var in liste_var_calc:
            Vm_k[var][mask_tot]+=Vm_temp[var]

            if (cf.singletype):
                Vm_t[var][mask_tot]=Vm_temp[var]
                Vm_t[var][~mask_tot] = np.NaN 

        
        del S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
        del el_temp, Tc_temp, Fw_temp, M_temp, Vm_temp, mask_tot

        # ===========  Dpol variables for single hydrometeor types =========== 
        if (cf.singletype):
            
            Vm_t["Zhh"] = np.copy(Vm_t["Zhhlin"])
            Vm_t["Zhh"][Vm_t["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_t["Zhhlin"][Vm_t["Zhhlin"]>0])
            Vm_t["Zdr"] = np.copy(Vm_t["Zhhlin"])
            Vm_t["Zdr"][(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                    (Vm_t["Zhhlin"]/Vm_t["Zvvlin"])[(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)])
            Vm_t["Rhohv"] = np.sqrt(np.divide(Vm_t["S11S22"], Vm_t["S11S11"]*Vm_t["S22S22"]))
            
            # ========== Writing dpol var for a single hydrometeor type t           
            fick = cf.pathfick+"k_"+cf.model+"_"+cf.band+'_'+str(int(cf.distmax_rad/1000.))+"_ech"+time+"_"+t
               
            if (cf.model=="Arome"):
                save.save_dpolvar_arome(liste_var_pol, Vm_t, Tc, Z,lat,lon,fick)
            elif (cf.model=="MesoNH"):
                save.save_dpolvar_mesonh(liste_var_pol, Vm_t, Tc, Z, X, Y,fick)
            else:
                print("model="+cf.model," => the save dpolvar option is available for Arome or MesoNH only")
            
            del Vm_t

    # ===== end loop over hydromet types =====================================    
    
    for var in liste_var_calc:
        Vm_k[var][~mask_precip_dist] = np.NaN 
   
    
    # ============ dpol var calculation
    Vm_k["Zhh"] = np.copy(Vm_k["Zhhlin"])
    Vm_k["Zhh"][Vm_k["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_k["Zhhlin"][Vm_k["Zhhlin"]>0])
    Vm_k["Zdr"] = np.copy(Vm_k["Zhhlin"])
    Vm_k["Zdr"][(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
            (Vm_k["Zhhlin"]/Vm_k["Zvvlin"])[(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)])
    Vm_k["Rhohv"] = np.sqrt(np.divide(Vm_k["S11S22"], Vm_k["S11S11"]*Vm_k["S22S22"]))
    

    
    # ============= Save dpol var for all hydromet in txt or npz file
    
    if (cf.model=="Arome"):
        save.save_dpolvar_arome(liste_var_pol, Vm_k, Tc, Z,lat,lon,fick)
    elif (cf.model=="MesoNH"):
        save.save_dpolvar_mesonh(liste_var_pol, Vm_k, Tc, Z, X, Y,fick)
    else:
        print("model="+cf.model," => the save dpolvar option is available for Arome or MesoNH only")

          
    del Vm_k

# end loop over time ============================================================
