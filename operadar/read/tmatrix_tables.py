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

import sys
import math
import time as tm
import numpy as np
import pandas as pd

from operadar.operadar_conf import LIMToption



def initialize_Tmatrix_dictionary() -> dict :
    """Initializes the Tmatrix dictionary to store all the variables and columns of the tables"""
    Tmatrix_params = {}
    list_of_parameters = ['LAMmin', 'LAMstep', 'LAMmax', 'ELEVmin', 'ELEVstep', 'ELEVmax',
                          'Tcmin', 'Tcstep', 'Tcmax','Fwmin', 'Fwstep', 'Fwmax', 'Tc_h',
                          'ELEV_h', 'Fw_h', 'M_h', 'S11carre_h', 'S22carre_h', 'ReS22S11_h',
                          'ImS22S11_h', 'ReS22fmS11f_h', 'ImS22ft_h', 'ImS11ft_h','RRint_h',
                          'expMmin', 'expMstep', 'expMmax', 'expCCmin', 'expCCstep', 'expCCmax',
                          ]
    
    for param in list_of_parameters :
        Tmatrix_params[param] = {}
        
    return Tmatrix_params



def get_scheme_to_fetch_table(microphysics:str) -> str :
    """Choice of the right table depending on the microphysics (will be modified later ?)"""

    if microphysics[0:3] == "ICE" :
        return "ICE3"
    elif microphysics[0:3] == "LIM" or (microphysics=="LIMT" and LIMToption=="cstmu") :
        return "LIMA"
    else :
        print('_____________')
        print('/!\ ERROR /!\ :',microphysics,'is not a valid name for Tmatrix computation')
        sys.exit()



def read_Tmatrix_Clotilde(band:str,
                          hydrometeors:list,
                          scheme:str,
                          pathTmat:str,
                          verbose:bool,
                          )-> dict:
    """Extract min/step/max in coefficient tables and other parameters from Clotilde's 2020 Tmatrix tables.
    
    Args:
        band (str): radar band.
        hydrometeors (list): list of hydrometeors for which Tmatrix tables must be read.
        scheme (str): microphysics scheme.
        pathTmat (str): Tmatrix directory path.

    Returns:
        Tmatrix_params (dict) : dictionary containing min/step/max values for multiple parameters
    """
    print("Reading Tmatrix tables")
    deb_timer = tm.time()
    micro_for_Tmatrix = get_scheme_to_fetch_table(microphysics=scheme)
    Tmatrix_params = initialize_Tmatrix_dictionary()
    
    for h in hydrometeors: 
        nomfileCoefInt = f'{pathTmat}TmatCoefInt_{micro_for_Tmatrix}_{band}{h}'
        
        if verbose : print("\tReading min/step/max for",h)
        df = pd.read_csv(nomfileCoefInt, sep=";",nrows = 1) 
        df['LAMstep'],df['LAMmax'] = df['LAMmin'],df['LAMmin'] # currently no loop over LAM so no LAMstep and LAMmax column (maybe later)
        values_to_get = ['LAMmin', 'LAMstep', 'LAMmax', 'ELEVmin', 'ELEVstep', 'ELEVmax',
                         'Tcmin', 'Tcstep', 'Tcmax','Fwmin', 'Fwstep', 'Fwmax',
                         'expMmin', 'expMstep', 'expMmax', 'expCCmin', 'expCCstep', 'expCCmax']
        for value in values_to_get :
            Tmatrix_params[value][h] = np.copy(df[value])[0]

        if verbose : print("\tReading scattering coefficients for",h)
        df_scat = pd.read_csv(nomfileCoefInt, sep=";",skiprows = [0, 1],dtype=np.float32)
        table_columns_to_read = ['Tc_h', 'ELEV_h', 'Fw_h', 'M_h', 'S11carre_h', 'S22carre_h', 'ReS22S11_h',
                                 'ImS22S11_h', 'ReS22fmS11f_h', 'ImS22ft_h', 'ImS11ft_h','RRint_h']
        for column in table_columns_to_read :
            if column == 'Fw_h':
                Tmatrix_params[column][h] = df_scat['P3'].to_numpy()
            else :
                Tmatrix_params[column][h] = df_scat[column[:-2]].to_numpy()
        del df_scat
    
    # WILL BE REMOVED LATER FOR COMPUTATION PURPOSE
    # for contents and number concentration: same min/step/max for all class of hydrometeors
    #same_value_for_all_hydrometeors = ['expMmin', 'expMstep', 'expMmax', 'expCCmin', 'expCCstep', 'expCCmax']
    #for value in same_value_for_all_hydrometeors :
    #    Tmatrix_params[value] = np.copy(df[value])[0]
    
    del df
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    return Tmatrix_params



def Read_VarTmatrixClotilde(pathTmat,band,schema_micro,table_ind,h):
        
    # Dictionnaries initialization
    LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax = {}, {}, {}, {}, {}, {}
    Tcmin, Tcstep, Tcmax = {}, {}, {}
    Fwmin, Fwstep, Fwmax = {}, {}, {}

    Tc_h, ELEV_h, Fw_h, M_h = {}, {}, {}, {}
    Zhh, Zdr, Rhv, Kdp = {}, {}, {}, {}
    
    nomfileVarInt = pathTmat+'TmatVarInt_'+schema_micro+'_'+band+'_'+h+table_ind  

    print("reading min/step/max in : ", nomfileVarInt)
    df = pd.read_csv(nomfileVarInt, sep=";",nrows = 1)
    LAMmin[h]= np.copy(df["LAMmin"])[0]
    LAMstep[h]= np.copy(df["LAMmin"])[0]
    LAMmax[h]= np.copy(df["LAMmin"])[0]
    ELEVmin[h]= np.copy(df["ELEVmin"])[0]
    ELEVstep[h]= np.copy(df["ELEVstep"])[0]
    ELEVmax[h]= np.copy(df["ELEVmax"])[0]
    Tcmin[h]= np.copy(df["Tcmin"])[0]
    Tcstep[h]= np.copy(df["Tcstep"])[0]
    Tcmax[h]= np.copy(df["Tcmax"])[0]
    Fwmin[h]= np.copy(df["Fwmin"])[0]
    Fwstep[h]= np.copy(df["Fwstep"])[0]
    Fwmax[h]= np.copy(df["Fwmax"])[0]

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
    Tc_h[h] = np.copy(df['Tc'])
    ELEV_h[h] = np.copy(df['ELEV'])
    Fw_h[h] = np.copy(df['P3'])
    M_h[h] = np.copy(df['M'])
    Zhh[h] = np.copy(df['Zhhlg'])
    Zdr[h] = np.copy(df['Zdrlg'])
    Rhv[h] = np.copy(df['Rhv'])
    Kdp[h] = np.copy(df['KDP'])

    del df
    # End loop over hydromet types    
   
    return LAMmin, LAMstep, LAMmax, ELEVmin, ELEVstep, ELEVmax, \
    Tcmin, Tcstep, Tcmax, Fwmin, Fwstep, \
    Fwmax,expMmin, expMstep, expMmax, expCCmin, expCCstep, expCCmax, \
    Tc_h, ELEV_h, Fw_h, M_h, Zhh, Zdr, Kdp, Rhv

    

def get_scatcoef(tmatDict:dict, hydrometeor:str, colName_tmatrix:str,
                 colMin:float, colStep:float, colMax:float,
                 el_temp:np.ndarray, Tc_temp:np.ndarray, colTmat:np.ndarray,
                 M_temp:np.ndarray, n_interpol:int, shutdown_warnings = True):    
    """Extract scattering coefficients for each class of hydrometeor
    
    Input:
    - full table of scattering coef for type t (S11carre_h[h] ...)
    - parameters of each column of the table (min, max, step): ELEVmin...
    - value of table columns in 3d arrays (el_temp, Tc_temp, colTmat, M_temp)
    
    Ouput :
    - scattering coef for type t within mask
    """
    S11carre_h = tmatDict['S11carre_h'][hydrometeor]
    S22carre_h = tmatDict['S22carre_h'][hydrometeor]
    ReS22fmS11f_h = tmatDict['ReS22fmS11f_h'][hydrometeor]
    ReS22S11_h = tmatDict['ReS22S11_h'][hydrometeor]
    ImS22S11_h = tmatDict['ImS22S11_h'][hydrometeor]
    LAMmin_h= tmatDict['LAMmin'][hydrometeor]
    LAMmax_h= tmatDict['LAMmax'][hydrometeor]
    LAMstep_h= tmatDict['LAMstep'][hydrometeor]
    ELEVmin_h= tmatDict['ELEVmin'][hydrometeor]
    ELEVmax_h= tmatDict['ELEVmax'][hydrometeor]
    ELEVstep_h= tmatDict['ELEVstep'][hydrometeor]
    Tcmin_h= tmatDict['Tcmin'][hydrometeor]
    Tcmax_h= tmatDict['Tcmax'][hydrometeor]
    Tcstep_h= tmatDict['Tcstep'][hydrometeor]
    expMmin = tmatDict['expMmin'][hydrometeor]
    expMstep = tmatDict['expMstep'][hydrometeor]
    expMmax = tmatDict['expMmax'][hydrometeor]
    
    
                                                            
    # Find position in the T-matrix table
    [kTmat, LAMred, ELEVred, Tcred, P3red, Mred] = CALC_KTMAT(el_temp, Tc_temp,colTmat, M_temp, LAMmin_h, LAMmax_h,
                                                              LAMstep_h, ELEVmin_h, ELEVmax_h, ELEVstep_h,
                                                              Tcmin_h, Tcmax_h, Tcstep_h, colMin, colMax, colStep,
                                                              expMmin, expMstep, expMmax, colName_tmatrix, shutdown_warnings)
    
     # Store scat coef values for each min/max born in Matcoef     
    MatCoef = {}
    
    for ind in list(range((n_interpol))):
        MatCoef[0, ind] = S11carre_h[kTmat[ind]]
        MatCoef[1, ind] = S22carre_h[kTmat[ind]]
        MatCoef[2, ind] = ReS22fmS11f_h[kTmat[ind]]
        MatCoef[3, ind] = ReS22S11_h[kTmat[ind]]
        MatCoef[4, ind] = ImS22S11_h[kTmat[ind]]
    
    # Interpol scat coef values
    [S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11] = INTERPOL(LAMred, ELEVred, Tcred, P3red, Mred, MatCoef)   

    return S11carre, S22carre, ReS22fmS11f, ReS22S11, ImS22S11



def  CALC_KTMAT(ELEV:float, Tc:np.ndarray, P3r, M:np.ndarray, LAMmin:float, LAMmax:float, LAMstep:float,
    ELEVmin:float, ELEVmax:float, ELEVstep:float, Tcmin:float, Tcmax:float, Tcstep:float,
    P3min, P3max, P3step, expMmin:float, expMstep:float, expMmax:float, P3name, shutdown_warnings) :
    """Return the indexes of the scattering coef (S11..) corresponding to the upper and lower bounds of
    LAM, ELEV, Tc, P3, M for => used for the interpolation of these coefficients.
       
    Args:
        ELEV (float) : elevation in radians
        Tc (np.ndarray) : temeprature in Celsius degree
        P3r : liquid water fraction (1-moment) or concentration (2-moment)
        M (np.ndarray) : contents in kg/m3
        LAMmin,max,step (float) : wavelength in mm
        ELEVmin,max,step (float) : degrees
        Tcmin,max,step (float) : Celsius degree
        P3min,max,step (float) : values between 0 and 1
        expMmin,step,max (float) : exponential of the concentration
        P3name : "Fw" liquid water fraction or "Nc" number concentration

    Returns:
        kTmat (dict) : dict with 32 tables with the M shape/size
    
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



def  INTERPOL(LAMred,ELEVred,Tcred,Fwred,Mred,MatCoef):
    """Multidimensional interpolation with outputs from CALC_KTMAT function

    Args:
        LAMred (_type_): wavelength
        ELEVred (_type_): elevation
        Tcred (_type_): temperature
        Fwred (_type_): liquid water fraction
        Mred (_type_): content
        MatCoef (_type_): 5*32 matrix with S11carre, S22carre, ReS22fmS11f, ReS22S11, and ImS22S11 coefficients to interpolate over the 32 bounds. 
    
    Returns:
        S11carre : description ?
        S22carre : description ?
        ReS22fmS11f : description ?
        ReS22S11 : description ?
        ImS22S11 : description ?
    """
    ncoef=5
    nval=Mred.shape[0]
    print("\t  INTERPOL func over nval=",nval)
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