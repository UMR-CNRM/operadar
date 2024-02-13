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
    * netcdf or .npz file with modele points i,j,k and dpol var Zh, Zdr, Kdp, Rhohv 
    
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
# - bug corrections in the mixed phase with a new option: MixedPhase="Fwpos" or "Tpos"
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
# implemented on github, see last changes on : https://github.com/UMR-CNRM/operadar/               
#=====================================================================================


# External modules
import sys
import math
import time as tm
import numpy as np
import pandas as pd
import datetime as dt

from pathlib import Path
sys.path.insert(0, "./lib")
from vortex_experiments import set_vortex_experiments

# 0perad modules
import operad_lib as ope_lib
import read_tmatrix as read_tmat
import save_dpolvar as save


# ===== Configuration files ===== #
#import common_settings as settings
import operad_conf as cf


# ===== Begining of the program ===== #
begining_program_timer = tm.time()

model = sys.argv[1]
micro = sys.argv[3]

if (model=="MesoNH"):
    import read_mesonh as meso
elif (model=="Arome"):
    import read_arome as aro


# ===== Reading dates in study cases csv file 
df = pd.read_csv(cf.csvPath, delimiter=";")
time_columns = ["start_time", "end_time"]
df[time_columns] = df[time_columns].apply(lambda x: pd.to_datetime(x, format="%Y%m%d%H%M"))

if sys.argv[2] == "all" :
    studyCases = df
else :
    date  = dt.datetime.strptime(sys.argv[2], "%Y%m%d").date() # date à laisser si plusieurs cas pour une même date , plutot pas mettre la date + l'heure de début ?
    studyCases = df.loc[df['start_time'].dt.date == date]

# ===== Loop over study cases
for _,row in studyCases.iterrows():
    run  = str(row.run_model).zfill(2)
    deb = row.start_time
    fin = row.end_time
    radar_ids = "-".join([str(x) for x in row.radar_id_list.strip().split(',')])
    radar_band = str(row.radar_band)
    
    # Zoom
    lat_min = row.latmin ; lat_max = row.latmax
    lon_min = row.lonmin ; lon_max = row.lonmax

    # Time list
    datetimelist=[]
    ech = deb
    while ech <= fin :
        datetimelist += [ech]
        ech += cf.step


    # ===== Testing existence of the output directory (or creates it) ===== #
    output_dir = Path(cf.outPath)
    if not output_dir.exists():
        try:
            output_dir.mkdir(exist_ok=True, parents=True)
            print ('Creating output directories :',output_dir)
        except:    
            print ('Error in creation of',output_dir) ; sys.exit()
    else:
        print ('Output directories exist :',output_dir)            


    liste_var_pol = ["Zhh", "Zdr", "Kdp","Rhohv"]
    liste_var_calc=["Zhhlin","Zvvlin","S11S22","S11S11","S22S22","Kdp","Rhohv"]


    # ===== Loop over timesteps ===== #
    read_tmatrix = True
    extract_once = True
    for datetime in datetimelist: 
        day=datetime.strftime('%Y%m%d')
        time = datetime.strftime('%H%M')
        outFile = cf.outPath + f"/dpolvar_{model}_{micro}_{radar_band}_{day}{time}"

        # ----- Testing existence of the output file ----- #
        if Path(outFile + ".nc").exists():
            print("netcdf file for",time,"already exists")
            continue   
        print("-------",datetime,"-------")
        
        # ----- Reading Tmatrix tables ----- #
        if read_tmatrix :
            print("Reading Tmatrix tables")
            deb_timer = tm.time()
   
            [LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax, 
            Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, Fwmax,
            expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, 
            Tc_t, ELEV_t, Fw_t, M_t, S11carre_t, S22carre_t,
            ReS22S11_t, ImS22S11_t, ReS22fmS11f_t, ImS22ft_t, ImS11ft_t] = read_tmat.Read_TmatrixClotilde(cf.pathTmat,
                                                                                                          radar_band,
                                                                                                          micro,
                                                                                                          cf.list_types_tot)
            LAM = LAMmin["rr"]/1000.
            read_tmatrix = False
            print("End reading Tmatrix tables in",round(tm.time()- deb_timer,2),"seconds")
            
        # ----- Reading model variables ----- #
        # return 3D model variables + coordinates
        print("Reading model variables") ; deb_timer = tm.time()  
        
        
        if (model=="MesoNH"):
            pathmodel = cf.commonPath
            datetime_run=dt.datetime.strptime(run,'%Y%m%d%H%M')
            model_ech=((datetime - datetime_run).total_seconds())/60.0/15.0
            ech=(str(int(model_ech))).zfill(3)
            #ech="027"
            print(ech)
            modelfile=pathmodel+cf.commonFilename+ech+".nc"
            pathfick = f"{cf.outPath}/k{cf.MixedPhase}"
            [M, Tc, CC, CCI, lat,lon, X, Y, Z, time] = meso.read_mesonh(modelfile = modelfile,
                                                            microphysics = micro,
                                                            hydrometeors_list = cf.htypes_model)
        elif (model=="Arome") :
            model_hour = (datetime - dt.timedelta(hours=int(run))).strftime('%H:%M')
            # Vortex experiment name
            expeOLIVE = set_vortex_experiments(run,micro)

            # Paths
            pathmodel = cf.commonPath + f"{expeOLIVE}/{deb.strftime('%Y%m%dT')}{run}00P/forecast/"
            pathfick = f"{cf.outPath}/{deb.strftime('%Y%m%d')}/{run}Z_{micro}_k{cf.MixedPhase}/{radar_ids}"
            modelfile=pathmodel+cf.commonFilename+model_hour+".fa"
            if extract_once :
                [M, Tc, CC, CCI, Z, lon, lat] = aro.read_arome(modelfile = modelfile,
                                                            microphysics = micro,
                                                            extract_once = extract_once,
                                                            lonmin = lon_min, lonmax = lon_max,
                                                            latmin = lat_min, latmax = lat_max,
                                                            hydrometeors_list = cf.htypes_model,
                                                            )
            else :
                [M, Tc, CC, CCI, Z] = aro.read_arome(modelfile = modelfile,
                                                    microphysics = micro,
                                                    extract_once = extract_once,
                                                    lonmin = lon_min, lonmax = lon_max,
                                                    latmin = lat_min, latmax = lat_max,
                                                    hydrometeors_list = cf.htypes_model,
                                                    )
        else:
            print("model = "+model+" => needs to be either Arome or MesoNH")
        extract_once = False    
        print("  --> Done in",round(tm.time()- deb_timer,2),"seconds")
        
        # ----- Compute radar geometry variables ----- #
        print("Compute radar geometry: elevation (el) and distance mask (mask_distmax)") ; deb_timer = tm.time()
        # TODO: add latlon2XY function in ope_lib
        # else: 
        #     X0,Y0=ope_lib.latlon2XY(cf.latrad,cf.lonrad)
        # ----------------
        np.seterr(invalid='ignore') # silence warning of invalid division 0 by 0 (result in a nan)

        if (model=="Arome"):
            el = np.zeros(Tc.shape)
            mask_distmax = (el >= 0.)
        
        elif (model=="MesoNH"):
            if (cf.radarloc=="center"):
                X0 = np.nanmean(X)
                Y0 = np.nanmean(Y)        
            [mask_distmax, el] = ope_lib.compute_radargeo(X, Y, Z, X0, Y0,
                                                        cf.distmax_rad,
                                                        cf.RT,ELEVmax["rr"],
                                                        )

        mask_precip_dist = ope_lib.mask_precip(mask_distmax, M, expMmin, micro) # precip mask
        [M, Fw] = ope_lib.compute_mixedphase(M, cf.MixedPhase, expMmin, micro) # mixed phase parametrization
        print("  --> Done in",round(tm.time()- deb_timer,2),"seconds")
        
        # Initialization of dict(Vm_k) --> contains all 3D dpol variables (all hydrometeor included)
        Vm_k = {var:np.zeros(Tc.shape) for var in liste_var_calc}
            
        # ----- Loop over hydromet types ----- #
        print("Loop over hydrometeor types in Tmatrix tables:",cf.list_types_tot) ; deb_timer = tm.time()
        for hydromet in cf.list_types_tot:
            if (cf.singletype):
                Vm_t = {var:np.zeros(Tc.shape) for var in liste_var_calc}
            
            # Compute NMOMENTS
            NMOMENTS = ope_lib.compute_nmoments(micro,hydromet)
            
            # Compute single type mask
            [mask_tot, M_temp,
            el_temp, Tc_temp, Fw_temp] = ope_lib.singletype_mask(M[hydromet], el, Tc,
                                                                Fw,mask_precip_dist,
                                                                expMmin, micro, hydromet
                                                                )
            
            # Define P3 : CC (2 moments) or Fw (1 moment)
            [P3, P3min, P3max, P3step] = ope_lib.defineP3(hydromet, NMOMENTS, CC, CCI,
                                                        mask_tot, Fw_temp,
                                                        expCCmin,expCCmax, expCCstep,
                                                        Fwmin[hydromet], Fwmax[hydromet], Fwstep[hydromet]
                                                        )
            
            # Extract scattering coefficients for singletype
            [S11carre,
            S22carre,
            ReS22fmS11f,
            ReS22S11,
            ImS22S11] = read_tmat.get_scatcoef(S11carre_t[hydromet],S22carre_t[hydromet],
                                                ReS22fmS11f_t[hydromet],ReS22S11_t[hydromet],ImS22S11_t[hydromet],
                                                LAMmin[hydromet], LAMmax[hydromet], LAMstep[hydromet],
                                                ELEVmin[hydromet], ELEVmax[hydromet], ELEVstep[hydromet],
                                                Tcmin[hydromet], Tcmax[hydromet], Tcstep[hydromet],
                                                P3min, P3max, P3step, expMmin,expMstep,expMmax,
                                                NMOMENTS, el_temp,Tc_temp,P3, M_temp,cf.n_interpol,shutdown_warnings=True,
                                                )
                
            # Single type dpol var computation       
            Vm_temp = {var:Vm_k[var][mask_tot] for var in liste_var_calc}
            Vm_temp["Zhhlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S22carre
            Vm_temp["Zvvlin"]= 1e18*LAM**4./(math.pi**5.*0.93)*4.*math.pi*S11carre
            Vm_temp["Kdp"] = 180.*1e3/math.pi*LAM*ReS22fmS11f
            Vm_temp["S11S22"] = ReS22S11**2+ImS22S11**2
            Vm_temp["S11S11"] = np.copy(S11carre)
            Vm_temp["S22S22"] = np.copy(S22carre)
            
            # Addition of scattering coef for all hydromet
            for var in liste_var_calc:
                Vm_k[var][mask_tot]+=Vm_temp[var]
                if (cf.singletype):
                    Vm_t[var][mask_tot]=Vm_temp[var]
                    Vm_t[var][~mask_tot] = np.NaN 
            
            del S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
            del el_temp, Tc_temp, Fw_temp, M_temp, Vm_temp, mask_tot

            #  Dpol variables for single hydrometeor types 
            if (cf.singletype):
                Vm_t["Zhh"] = np.copy(Vm_t["Zhhlin"])
                Vm_t["Zhh"][Vm_t["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_t["Zhhlin"][Vm_t["Zhhlin"]>0])
                Vm_t["Zdr"] = np.copy(Vm_t["Zhhlin"])
                Vm_t["Zdr"][(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                        (Vm_t["Zhhlin"]/Vm_t["Zvvlin"])[(Vm_t["Zhhlin"]>0) & (Vm_t["Zvvlin"]>0)])
                Vm_t["Rhohv"] = np.sqrt(np.divide(Vm_t["S11S22"], Vm_t["S11S11"]*Vm_t["S22S22"]))
                
                # Writing dpol var for a single hydrometeor type hydromet           
                fick = pathfick+"k_"+model+"_"+radar_band+'_'+str(int(cf.distmax_rad/1000.))+"_ech"+time+"_"+hydromet
                
                if (model=="Arome") :
                    save.save_dpolvar_arome(liste_var_pol, Vm_t, Tc, Z,lat,lon,fick,datetime,cf.save_npz,cf.save_netcdf)
                elif (model=="MesoNH") :
                    save.save_dpolvar_mesonh(liste_var_pol, Vm_t, Tc, Z, X, Y,fick,cf.save_npz,cf.save_netcdf)
                else:
                    print("model = "+model," => the save dpolvar option is available for Arome or MesoNH only")
                del Vm_t
    
        for var in liste_var_calc:
            Vm_k[var][~mask_precip_dist] = np.NaN 
        
        # ----- dpol var calculation
        Vm_k["Zhh"] = np.copy(Vm_k["Zhhlin"])
        Vm_k["Zhh"][Vm_k["Zhhlin"]>0] = ope_lib.Z2dBZ(Vm_k["Zhhlin"][Vm_k["Zhhlin"]>0])
        Vm_k["Zdr"] = np.copy(Vm_k["Zhhlin"])
        Vm_k["Zdr"][(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)] = ope_lib.Z2dBZ( \
                (Vm_k["Zhhlin"]/Vm_k["Zvvlin"])[(Vm_k["Zhhlin"]>0) & (Vm_k["Zvvlin"]>0)])
        Vm_k["Rhohv"] = np.sqrt(np.divide(Vm_k["S11S22"], Vm_k["S11S11"]*Vm_k["S22S22"]))
        print("  --> Done in",round(tm.time() - deb_timer,2),"seconds")
        
        # ----- Save dpol var for all hydromet in netcdf and/or npz file
        if (model=="Arome"):
            save.save_dpolvar_arome(liste_var_pol,M, CC, CCI, Vm_k, Tc, Z,lat,lon,outFile,datetime,cf.save_npz,cf.save_netcdf)
        elif (model=="MesoNH"):
            save.save_dpolvar_mesonh(liste_var_pol, Vm_k, Tc, Z, X, Y,lat,lon,time,outFile,cf.save_npz,cf.save_netcdf)
        else:
            print("model = "+model," => the save dpolvar option is available for Arome or MesoNH only")
    
        del Vm_k

end_program_timer = tm.time()
elapsed_time = end_program_timer - begining_program_timer
print("End of the program in",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")
