#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 09:55:15 2023

@author: augros

Contains routines to :
* read tmatrix tables: 
  - Read_TmatrixClotilde
  - Read_VarTmatrixClotilde
* extract the scattering coefficient 
  - get_scatcoef
  - CALC_KTMAT
  - INTERPOL
        
"""

import numpy as np
import pandas as pd
import math
import sys

from operad_conf import micro_scheme, LIMToption

#============== Read_Tmatrix2020
"""
Reads Clotilde's 2020 Tmatrix tables:
    direct reading of min/step: max in coefficient files)
Output: 
min/step/max parameters
+ Tc_t, ELEV_t, Fw_t, M_t
+ S11carre_t, S22carre_t, ReS22S11_t,
  ImS22S11_t, ReS22fmS11f_t, ImS22ft_t, ImS11ft_tS11carre_t  
"""

def Read_TmatrixClotilde(pathTmat,bande,hydromet_list:list):

    # Choice of the right table depending on the microphysics
    if micro_scheme[0:3] == "ICE" :
        micro_for_Tmatrix = "ICE3"
    elif micro_scheme[0:3] == "LIM" or (micro_scheme=="LIMT" and LIMToption=="cstmu") :
        micro_for_Tmatrix = "LIMA"
    else :
        print('_____________')
        print('/!\ ERROR /!\ :',micro_scheme,'is not a valid name for Tmatrix computation')
        sys.exit()

    # Dictionnaries initialization
    LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax = {}, {}, {}, {}, {}, {}
    Tcmin, Tcstep, Tcmax = {}, {}, {}
    Fwmin, Fwstep, Fwmax = {}, {}, {}
    # LAM_t = {}
    Tc_t, ELEV_t, Fw_t, M_t = {}, {}, {}, {}
    S11carre_t, S22carre_t = {}, {}
    ReS22S11_t, ImS22S11_t = {}, {}
    ReS22fmS11f_t, ImS22ft_t, ImS11ft_t,RRint_t = {}, {}, {},{}
        
#    if bande =='S':
#        LAMstr='106.2'
#    elif bande=='C':
#        LAMstr="053.2"
    
    for t in hydromet_list: 
        #nomfileCoefInt = pathTmat+'TmatCoefInt_'+micro_for_Tmatrix+'_'+bande+LAMstr+'_'+t+table_ind
        #nomfileCoefInt = pathTmat+'TmatCoefInt_'+micro_for_Tmatrix+'_'+bande+'_'+t
        nomfileCoefInt = pathTmat+'TmatCoefInt_'+micro_for_Tmatrix+'_'+bande+t
        
        print("\tReading min/step/max for",t)
        df = pd.read_csv(nomfileCoefInt, sep=";",nrows = 1)
        LAMmin[t]= np.copy(df["LAMmin"])[0]
        LAMstep[t]= np.copy(df["LAMmin"])[0]
        LAMmax[t]= np.copy(df["LAMmin"])[0]
        ELEVmin[t]= np.copy(df["ELEVmin"])[0]
        ELEVstep[t]= np.copy(df["ELEVstep"])[0]
        ELEVmax[t]= np.copy(df["ELEVmax"])[0]
        Tcmin[t]= np.copy(df["Tcmin"])[0]
        Tcstep[t]= np.copy(df["Tcstep"])[0]
        Tcmax[t]= np.copy(df["Tcmax"])[0]
        Fwmin[t]= np.copy(df["Fwmin"])[0]
        Fwstep[t]= np.copy(df["Fwstep"])[0]
        Fwmax[t]= np.copy(df["Fwmax"])[0]

        print("\tReading scattering coef for",t)
        df_scat = pd.read_csv(nomfileCoefInt, sep=";",skiprows = [0, 1],dtype=np.float32)
        Tc_t[t] = df_scat['Tc'].to_numpy()
        ELEV_t[t] = df_scat['ELEV'].to_numpy()
        Fw_t[t] = df_scat['P3'].to_numpy()
        M_t[t] = df_scat['M'].to_numpy()
        S11carre_t[t] = df_scat['S11carre'].to_numpy()
        S22carre_t[t] = df_scat['S22carre'].to_numpy()
        ReS22S11_t[t] = df_scat['ReS22S11'].to_numpy()
        ImS22S11_t[t] = df_scat['ImS22S11'].to_numpy()
        ReS22fmS11f_t[t] = df_scat['ReS22fmS11f'].to_numpy()
        ImS22ft_t[t] = df_scat['ImS22ft'].to_numpy()
        ImS11ft_t[t] = df_scat['ImS11ft'].to_numpy()
        RRint_t[t] = df_scat['RRint'].to_numpy()
        del df_scat
        
    # For M and CC: same min/step/max for all types
    expMmin= np.copy(df["expMmin"])[0]
    expMstep= np.copy(df["expMstep"])[0]
    expMmax= np.copy(df["expMmax"])[0]
    expCCmin= np.copy(df["expCCmin"])[0]
    expCCstep= np.copy(df["expCCstep"])[0]
    expCCmax= np.copy(df["expCCmax"])[0]
    del df
   
    return LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax,Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, \
    Fwmax,expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, Tc_t, ELEV_t, Fw_t, M_t, S11carre_t, \
    S22carre_t, ReS22S11_t, ImS22S11_t, ReS22fmS11f_t, ImS22ft_t, ImS11ft_t


def Read_VarTmatrixClotilde(pathTmat,bande,schema_micro,table_ind,t):
        
    # Dictionnaries initialization
    LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax = {}, {}, {}, {}, {}, {}
    Tcmin, Tcstep, Tcmax = {}, {}, {}
    Fwmin, Fwstep, Fwmax = {}, {}, {}

    Tc_t, ELEV_t, Fw_t, M_t = {}, {}, {}, {}
    Zhh, Zdr, Rhv, Kdp = {}, {}, {}, {}
    
    nomfileVarInt = pathTmat+'TmatVarInt_'+schema_micro+'_'+bande+'_'+t+table_ind  

    print("reading min/step/max in : ", nomfileVarInt)
    df = pd.read_csv(nomfileVarInt, sep=";",nrows = 1)
    LAMmin[t]= np.copy(df["LAMmin"])[0]
    LAMstep[t]= np.copy(df["LAMmin"])[0]
    LAMmax[t]= np.copy(df["LAMmin"])[0]
    ELEVmin[t]= np.copy(df["ELEVmin"])[0]
    ELEVstep[t]= np.copy(df["ELEVstep"])[0]
    ELEVmax[t]= np.copy(df["ELEVmax"])[0]
    Tcmin[t]= np.copy(df["Tcmin"])[0]
    Tcstep[t]= np.copy(df["Tcstep"])[0]
    Tcmax[t]= np.copy(df["Tcmax"])[0]
    Fwmin[t]= np.copy(df["Fwmin"])[0]
    Fwstep[t]= np.copy(df["Fwstep"])[0]
    Fwmax[t]= np.copy(df["Fwmax"])[0]

    # For M and CC: same min/step/max for all types
    expMmin= np.copy(df["expMmin"])[0]
    expMstep= np.copy(df["expMstep"])[0]
    expMmax= np.copy(df["expMmax"])[0]
    expCCmin= np.copy(df["expCCmin"])[0]
    expCCstep= np.copy(df["expCCstep"])[0]
    expCCmax= np.copy(df["expCCmax"])[0]
    
    del df        
    
    print("reading dpol variables in : ", nomfileVarInt)
    df = pd.read_csv(nomfileVarInt, sep=";",skiprows = [0, 1])
    #df = pandas.read_csv(nomfileCoefInt, sep=";",names=["Tc_t", "ELEV_t", "Fw_t", "M_t", "S11carre_t", "S22carre_t", "ReS22S11_t", "ImS22S11_t", "ReS22fmS11f_t", "ImS22ft_t", "ImS11ft_t", "RRint_t"])
    Tc_t[t] = np.copy(df['Tc'])
    ELEV_t[t] = np.copy(df['ELEV'])
    Fw_t[t] = np.copy(df['P3'])
    M_t[t] = np.copy(df['M'])
    Zhh[t] = np.copy(df['Zhhlg'])
    Zdr[t] = np.copy(df['Zdrlg'])
    Rhv[t] = np.copy(df['Rhv'])
    Kdp[t] = np.copy(df['KDP'])

    del df
    # End loop over hydromet types    
   
    return LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax, \
    Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, \
    Fwmax,expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, \
    Tc_t, ELEV_t, Fw_t, M_t, Zhh, Zdr, Kdp, Rhv
# =============================================================================
    

# ========= Extract scattering coefficients for singletype ====================
"""
 * input:
     - full table of scattering coef for type t (S11carre_t[t] ...)
     - parameters of each column of the table (min, max, step): ELEVmin...
     - value of table columns in 3d arrays (el_temp, Tc_temp, P3, M_temp)
 * ouput :
     - scattering coef for type t within mask
"""
def get_scatcoef(S11carre_tt,S22carre_tt,ReS22fmS11f_tt,ReS22S11_tt,ImS22S11_tt,\
                 LAMmint, LAMmaxt, LAMstept,\
                 ELEVmint, ELEVmaxt, ELEVstept,\
                 Tcmint, Tcmaxt, Tcstept, P3min, P3max, P3step,\
                 expMmin,expMstep,expMmax,\
                 P3name, el_temp,Tc_temp,P3, M_temp,n_interpol,shutdown_warnings = False):    
                                                               
    # Find position in the T-matrix table
    [kTmat, LAMred, ELEVred, Tcred, P3red, Mred] = CALC_KTMAT(el_temp,\
        Tc_temp,P3, M_temp, LAMmint, LAMmaxt,\
        LAMstept, ELEVmint, ELEVmaxt, ELEVstept, Tcmint, Tcmaxt, Tcstept, P3min, P3max, P3step,\
        expMmin,expMstep,expMmax,P3name,shutdown_warnings)
    
     # Store scat coef values for each min/max born in Matcoef     
    MatCoef = {}
    
    for ind in list(range((n_interpol))):
        MatCoef[0, ind] = S11carre_tt[kTmat[ind]]
        MatCoef[1, ind] = S22carre_tt[kTmat[ind]]
        MatCoef[2, ind] = ReS22fmS11f_tt[kTmat[ind]]
        MatCoef[3, ind] = ReS22S11_tt[kTmat[ind]]
        MatCoef[4, ind] = ImS22S11_tt[kTmat[ind]]
    
    # Interpol scat coef values
    [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11] = INTERPOL(LAMred, ELEVred, Tcred, P3red, Mred, MatCoef)   

    return S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11
# =============================================================================    



#============== Fonction CALC_KTMAT

#def  CALC_KTMAT(LAMm,ELEV,Tc,P3r,M,LAMmin,LAMmax,LAMstep,
#    ELEVmin,ELEVmax,ELEVstep,Tcmin,Tcmax,Tcstep,
#    P3min,P3max,P3step,expMmin,expMstep,expMmax,NMOMENTS):
def  CALC_KTMAT(ELEV,Tc,P3r,M,LAMmin,LAMmax,LAMstep,
    ELEVmin,ELEVmax,ELEVstep,Tcmin,Tcmax,Tcstep,
    P3min,P3max,P3step,expMmin,expMstep,expMmax,P3name,shutdown_warnings):
    
    """
       ELEV = Elevation en radians (tab)
       Tc = température en °C (tab)
       P3r = liquid water fraction (1-moment) or concentration (2-moment)
       M = contenu kg/m3 (tab)
       LAMmin,max,step en mm (val)
       ELEVmin,max,step en degrees (val)
       Tcmin,max,step en °C (val)
       P3min,max,step entre 0 et 1 (val)
       expMmin,expMstep,expMmax
       P3name ("Fw" liquid water fraction or "Nc" number concentration)

       return: kTmat, dict with 32 tables with the M shape/size : indexes of the scattering coef (S11..)
       corresponding to the upper and lower bounds of LAM, ELEV, Tc, P3, M for 
       => used for the interpolation of these coefficients    
    """

     
    kTmat={}
    
    for n in range(32):
        #kTmat[n]=np.empty((nk,nj,ni))
        kTmat[n]=np.copy(M)*float('nan') 
    
    #Conversion de LAM de m en mm
#    LAM=LAMm*10**3
    LAM=LAMmin # No need interpolation with LAM (fixed in the Tmatrix tables)
    
    #Hydromet content
    #expM=np.ones((Min.shape[0]))*(-100.)
    expM=np.copy(M)*0.0-100
    expM[M>0]=(np.log10(M))[M>0]
    
    # 1 moment: P3=Liquid water fraction
    if (P3name=="Fw"):
        P3=np.copy(P3r)
    # 2 moments: P3=number concentration   
    else :
        if (P3name=="Nc"): 
            P3=np.copy(P3r)
            P3[P3r>0]=np.log10(P3r[P3r>0])
            P3[P3r==0]=P3min
        else:
            print("Error P3name (3d parameter in Tmatrix table), the only options are Fw (liq water fraction) or Nc (number concentration)")
            sys.exit(1)
    
    # If LAM, ELEV, Tc, P3 or M are outside min and max ranges:
    # warning and the values are set to the min (if below min) or max (if over max)

    #if (abs(LAM-LAMmin) < LAMstep/10):
    if (LAM<LAMmin):
        if shutdown_warnings == False : print("Warning : LAM = ",LAM, " < LAMmin=",LAMmin)
        LAM=LAMmin
    #if (abs(LAM-LAMmax) < LAMstep/10):
    if (LAM > LAMmax):
        if shutdown_warnings == False : print("Warning : LAM = ",LAM, " > LAMmax=",LAMmax)
        LAM=LAMmax
    
    if (len(ELEV[ELEV<ELEVmin])>0):
        if shutdown_warnings == False : print("Warning: ELEV < ELEVmin: ",ELEV[ELEV<ELEVmin])
    ELEV[ELEV<ELEVmin]=ELEVmin
    
    if (len(ELEV[ELEV>ELEVmax])>0):
        if shutdown_warnings == False : print("Warning: ELEV > ELEVmax: ",ELEV[ELEV>ELEVmax])
    ELEV[ELEV>ELEVmax]=ELEVmax
    
    if (len(Tc[Tc<Tcmin])>0):
        if shutdown_warnings == False : print("Warning: Tc < Tcmin: ",Tc[Tc<Tcmin])
    Tc[Tc<Tcmin]=Tcmin
    
    if (len(Tc[Tc>Tcmax])>0):
        if shutdown_warnings == False : print("Warning: Tc > Tcmax: ",Tc[Tc>Tcmax])
    Tc[Tc>Tcmax]=Tcmax
    
    if (len(P3[P3<P3min])>0):
        if shutdown_warnings == False : print("Warning: P3 < P3min: ",P3[P3<P3min])
    P3[P3<P3min]=P3min
 
    if (len(P3[P3>P3max])>0):
        if shutdown_warnings == False : print("Warning: P3 > P3max: ",P3[P3>P3max])
    P3[P3>P3max]=P3max

    if (len(expM[expM<expMmin])>0):
        if shutdown_warnings == False : print("Warning: expM < expMmin: ",expM[expM<expMmin])

    if (len(expM[expM>expMmax])>0):
        if shutdown_warnings == False : print("Warning: expM > expMmax: ",expM[expM>expMmax])
    expM[expM>expMmax]=expMmax
    
    #    condok=np.where((LAM >=LAMmin) & (LAM<=LAMmax) & (ELEV>=ELEVmin) & (ELEV<=ELEVmax) 
    #                     & (Tc >=Tcmin) & (Tc<=Tcmax) & (P3 >=P3min) & (P3<=P3max)
    #                     & (expM >=expMmin) & (expM<=expMmax) & (dist<=distmax_mod),1.,float('NaN'))
    condok=np.where((LAM >=LAMmin) & (LAM<=LAMmax) & (ELEV>=ELEVmin) & (ELEV<=ELEVmax) 
                     & (Tc >=Tcmin) & (Tc<=Tcmax) & (P3 >=P3min) & (P3<=P3max)
                     & (expM >=expMmin) & (expM<=expMmax),1.,float('NaN'))
    
    
    if shutdown_warnings == False : 
        print("Tc isnan: ",Tc[np.isnan(condok)])
        print("expM isnan: ",expM[np.isnan(condok)])
        print("P3 isnan: ",P3[np.isnan(condok)])
    
    # Looking for the location in the table of the values given as input in the 
    # function
    # kP3, kTc, kexpM are numpy arrays with the shape/size of P3, Tc, expM
    # with condok, values outside the min/max range are set to NaN
    #------- LAM ------------------
    kLAM=condok*math.floor((LAM-LAMmin)/LAMstep)
    LAMinf=LAMmin+kLAM*LAMstep
    kLAMs=np.copy(kLAM)
    kLAMs[LAM!=LAMinf]+=1
    
    LAMsup=LAMmin+kLAMs*LAMstep
    nLAM=condok*(((LAMmax-LAMmin)/LAMstep)+1)
    
    
    #------- ELEV ------------------
    kELEV=condok*(np.floor((ELEV-ELEVmin)/ELEVstep))
    ELEVinf=ELEVmin+kELEV*ELEVstep
    kELEVs=np.copy(kELEV)
    kELEVs[ELEV!=ELEVinf]+=1
    
    ELEVsup=ELEVmin+kELEVs*ELEVstep
    nELEV=condok*(int((ELEVmax-ELEVmin)/ELEVstep)+1)
    
    #------- Tc ------------------
    kTc=condok*(np.floor((Tc-Tcmin)/Tcstep))
    Tcinf=Tcmin+kTc*Tcstep
    kTcs=np.copy(kTc)
    kTcs[Tc!=Tcinf]+=1
    
    Tcsup=Tcmin+(kTcs)*Tcstep
    nTc=condok*(int((Tcmax-Tcmin)/Tcstep)+1)
    
    
    #------- P3 ------------------
    kP3=condok*(np.floor((P3-P3min)/P3step))
    P3inf=P3min+kP3*P3step
    kP3s=np.copy(kP3)
    kP3s[P3!=P3inf]+=1
    
    P3sup=P3min+(kP3s)*P3step
    nP3=condok*(int((P3max-P3min)/P3step)+1)
    
    
    #------- M ------------------
    kexpM=condok*(np.floor((expM-expMmin)/expMstep))
    expMinf=expMmin+kexpM*expMstep
    kexpMs=np.where(expM==expMinf,kexpM,kexpM+1)
    kexpMs=np.copy(kexpM)
    kexpMs[expM!=expMinf]+=1
    
    expMsup=expMmin+kexpMs*expMstep
    nM=condok*(int((expMmax-expMmin)/expMstep)+1)
    Minf=10**expMinf
    Msup=10**expMsup
    
    
    #-- Calcul des variables reduites (comprises entre 0 et 1)
    # pour l'interpolation linaire
    LAMred=np.copy(LAMinf)*0.0
    ELEVred=np.copy(ELEV)*0.0
    Tcred=np.copy(Tc)*0.0
    P3red=np.copy(P3)*0.0
    Mred=np.copy(M)*0.0
    
    msk=LAMsup>LAMinf
    LAMred[msk]=(LAM-LAMinf[msk])/(LAMsup[msk]-LAMinf[msk])
    
    msk=ELEVsup>ELEVinf
    ELEVred[msk]=(ELEV[msk]-ELEVinf[msk])/(ELEVsup[msk]-ELEVinf[msk])
    
    msk=Tcsup>Tcinf
    Tcred[msk]=(Tc[msk]-Tcinf[msk])/(Tcsup[msk]-Tcinf[msk])
    
    msk=P3sup>P3inf
    P3red[msk]=(P3[msk]-P3inf[msk])/(P3sup[msk]-P3inf[msk])
    
    msk=Msup>Minf
    Mred[msk]=(M[msk]-Minf[msk])/(Msup[msk]-Minf[msk])
    
    
        
    
    # On regroupe autant que possible les produits et sommes dans le calcul de kTmat. Les calculs originaux 
    # sont en commentaires un peu plus bas.
    fact1=nP3*nM
    fact2=nELEV*fact1
    fact3=nTc*fact2
    fact41=kLAM*fact3
    fact42=kLAMs*fact3
    fact51=kTc*fact2
    fact52=kTcs*fact2
    fact61=kELEV*fact1
    fact62=kELEVs*fact1
    fact71=kP3*nM
    fact72=kP3s*nM
    
    som11=fact41+fact51
    som111=som11+fact61
    som112=som11+fact62
    som1111=som111+fact71
    som1112=som111+fact72
    som1121=som112+fact71
    som1122=som112+fact72
    som12=fact41+fact52
    som121=som12+fact61
    som122=som12+fact62
    som1211=som121+fact71
    som1212=som121+fact72
    som1221=som122+fact71
    som1222=som122+fact72
    som21=fact42+fact51
    som211=som21+fact61
    som212=som21+fact62
    som2111=som211+fact71
    som2112=som211+fact72
    som2121=som212+fact71
    som2122=som212+fact72
    som22=fact42+fact52
    som221=som22+fact61
    som222=som22+fact62
    som2211=som221+fact71
    som2212=som221+fact72
    som2221=som222+fact71
    som2222=som222+fact72
    
    kTmat[0]=som1111+kexpM+1
    kTmat[1]=som1111+kexpMs+1
    kTmat[2]=som1112+kexpM+1
    kTmat[3]=som1112+kexpMs+1
    kTmat[4]=som1121+kexpM+1
    kTmat[5]=som1121+kexpMs+1
    kTmat[6]=som1122+kexpM+1
    kTmat[7]=som1122+kexpMs+1
    kTmat[8]=som1211+kexpM+1
    kTmat[9]=som1211+kexpMs+1
    kTmat[10]=som1212+kexpM+1
    kTmat[11]=som1212+kexpMs+1
    kTmat[12]=som1221+kexpM+1
    kTmat[13]=som1221+kexpMs+1
    kTmat[14]=som1222+kexpM+1
    kTmat[15]=som1222+kexpMs+1
    kTmat[16]=som2111+kexpM+1
    kTmat[17]=som2111+kexpMs+1
    kTmat[18]=som2112+kexpM+1
    kTmat[19]=som2112+kexpMs+1
    kTmat[20]=som2121+kexpM+1
    kTmat[21]=som2121+kexpMs+1
    kTmat[22]=som2122+kexpM+1
    kTmat[23]=som2122+kexpMs+1
    kTmat[24]=som2211+kexpM+1
    kTmat[25]=som2211+kexpMs+1
    kTmat[26]=som2212+kexpM+1
    kTmat[27]=som2212+kexpMs+1
    kTmat[28]=som2221+kexpM+1
    kTmat[29]=som2221+kexpMs+1
    kTmat[30]=som2222+kexpM+1
    kTmat[31]=som2222+kexpMs+1
    
    
    for i in range(32):
        kTmat[i][np.isnan(kTmat[i])]=-999
        kTmat[i]=kTmat[i].astype(int)
 
    return kTmat, LAMred, ELEVred,Tcred,P3red,Mred

    del ELEV,Tc,expM,kLAM,LAMinf,kLAMs,LAMsup,nLAM,kELEV,ELEVinf,kELEVs,ELEVsup,nELEV
    del kTc,Tcinf,kTcs,Tcsup,nTc,kP3,P3inf,kP3s,P3sup,nP3,kexpM,expMinf,kexpMs,expMsup,nM,Minf,Msup

#===== Fin fonction CALC_KTMAT


#========= Fonction INTERPOL ========================
"""
LAMred,ELEVred,Tcred,Fwred,Mred en sortie de CALC_KTMAT
MatCoef !matrice 5*32 contenant les coef à interpoler sur les 32 bornes: S11carre, S22carre,
ReS22fmS11f, ReS22S11 et ImS22S11
return: S11carre,S22carre,ReS22fmS11f, ReS22S11 et ImS22S11 
"""
def  INTERPOL(LAMred,ELEVred,Tcred,Fwred,Mred,MatCoef):
    
    ncoef=5
    nval=Mred.shape[0]
    print("  INTERPOL: nval=",nval)
    VectCoef={}

    #--- Interpolation linéaire ---
    for indcoef in range(ncoef):
        VectCoef[indcoef]=np.empty((nval))
        VectCoef[indcoef]= \
           ( (1-LAMred)*((1-Tcred)*((1-ELEVred)*((1-Fwred)*((1-Mred)*MatCoef[indcoef,0]+
                Mred*MatCoef[indcoef,1])+
                    Fwred*((1-Mred)*MatCoef[indcoef,2]+
                        Mred*MatCoef[indcoef,3]))+
                            ELEVred*((1-Fwred)*((1-Mred)*MatCoef[indcoef,4]+
                                Mred*MatCoef[indcoef,5])+
                                    Fwred*((1-Mred)*MatCoef[indcoef,6]+
                                        Mred*MatCoef[indcoef,7])))+
                                        Tcred*((1-ELEVred)*((1-Fwred)*((1-Mred)*MatCoef[indcoef,8]+
            Mred*MatCoef[indcoef,9])+
                Fwred*((1-Mred)*MatCoef[indcoef,10]+
                    Mred*MatCoef[indcoef,11]))+
                        ELEVred*((1-Fwred)*((1-Mred)*MatCoef[indcoef,12]+
                            Mred*MatCoef[indcoef,13])+
                                Fwred*((1-Mred)*MatCoef[indcoef,14]+
                                    Mred*MatCoef[indcoef,15]))))+
                                    LAMred*((1-Tcred)*((1-ELEVred)*((1-Fwred)*((1-Mred)*MatCoef[indcoef,16]+
            Mred*MatCoef[indcoef,17])+
                Fwred*((1-Mred)*MatCoef[indcoef,18]+
                    Mred*MatCoef[indcoef,19]))+
                        ELEVred*((1-Fwred)*((1-Mred)*MatCoef[indcoef,20]+
                            Mred*MatCoef[indcoef,21])+
                                Fwred*((1-Mred)*MatCoef[indcoef,22]+
                                    Mred*MatCoef[indcoef,23])))+
                                    Tcred*((1-ELEVred)*((1-Fwred)*((1-Mred)*MatCoef[indcoef,24]+
            Mred*MatCoef[indcoef,25])+
                Fwred*((1-Mred)*MatCoef[indcoef,26]+
                    Mred*MatCoef[indcoef,27]))+
                        ELEVred*((1-Fwred)*((1-Mred)*MatCoef[indcoef,28]+
                            Mred*MatCoef[indcoef,29])+
                                Fwred*((1-Mred)*MatCoef[indcoef,30]+
                                    Mred*MatCoef[indcoef,31]))))) 
        #--- fin interpolation lineaire -----------------------
        #print indcoef,np.count_nonzero(~np.isnan(VectCoef[indcoef]))

    S11carre=np.copy(VectCoef[0])
    S22carre=np.copy(VectCoef[1])
    ReS22fmS11f=np.copy(VectCoef[2])
    ReS22S11=np.copy(VectCoef[3])
    ImS22S11=np.copy(VectCoef[4])
                     
    return S11carre,S22carre,ReS22fmS11f, ReS22S11, ImS22S11

    del VectCoef
# ======= Fin fonction INTERPOL
    
