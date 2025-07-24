#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augros & davidcl        
"""

import sys
import time as tm
import numpy as np
import pandas as pd



def initialize_table_dictionary(dpol2add:list,method:str='Tmatrix') -> dict :
    """Initializes the dictionary to store all the necessary parameters and columns of the tables"""
    empty_table_dict = {}
    list_of_parameters = ['LAM',
                          'ELEVmin', 'ELEVstep', 'ELEVmax',
                          'Tcmin', 'Tcstep', 'Tcmax',
                          'Fwmin', 'Fwstep', 'Fwmax', 
                          'expMmin', 'expMstep', 'expMmax',
                          'expCCmin', 'expCCstep', 'expCCmax',
                          ]
    columns_to_retrieve = retrieve_needed_columns(dpol2add=dpol2add,
                                                  scattering_method=method,
                                                  )
    for key in list_of_parameters+columns_to_retrieve :
        empty_table_dict[key] = {}
        
    return empty_table_dict, list_of_parameters, columns_to_retrieve



def retrieve_needed_columns(dpol2add:list,scattering_method:str='Tmatrix')->list:
    """Depending on the variables the user wants to compute and the chosen scattering method,
    creating a list of the table's column names to extract."""
    
    table_columnNames = ['Tc', 'ELEV', 'M', 'Fw', 'N']
    
    if scattering_method == 'Tmatrix' or scattering_method == 'both' :
        if 'Zh' in dpol2add :
            table_columnNames += ['sighh']
        if 'Zdr' in dpol2add :
            table_columnNames += ['sighh','sigvv']
        if 'Kdp' in dpol2add :
            table_columnNames += ['kdp']
        if 'Rhohv' in dpol2add :
            table_columnNames += ['REdeltaco','IMdeltaco','sighh','sigvv']
        if 'Ah' in dpol2add :
            table_columnNames += ['Ah']
        if 'Av' in dpol2add :
            table_columnNames += ['Av']
            
    if scattering_method == 'Rayleigh' or scattering_method == 'both' :
        if 'Zh' in dpol2add :
            table_columnNames += ['sighhR']
        if 'Zdr' in dpol2add :
            table_columnNames += ['sighhR'',sigvvR']
        if 'Kdp' in dpol2add :
            table_columnNames += ['kdpR']
        if 'Rhohv' in dpol2add :
            table_columnNames += ['REdeltacoR','IMdeltacoR','sighhR','sigvvR']
        if 'Ah' in dpol2add :
            table_columnNames += ['AhR']
        if 'Av' in dpol2add :
            table_columnNames += ['AvR']    
    
    return list(set(table_columnNames)) # return unique list of strings



def cloud_water_species(hydrometeor:str, cloud_water_over:str) -> str :
    """Handle the different name of the cloud water lookup tables (there is
    one for cloud water over sea and one for cloud water over land)"""
    if cloud_water_over == 'land' :
        hydrometeor = 'cl'
    elif cloud_water_over == 'sea' :
        hydrometeor = 'cs'
    else :
        print('_____________')
        print('/!\\ ERROR /!\\ :',cloud_water_over,'is not a valid option for cloud water')
        print('                can only be cloud water over "sea" or "land"')
        sys.exit()
    return hydrometeor



def read_and_extract_tables_content(band:str,
                                    hydrometeors:list,
                                    moments:dict,
                                    cloud_water_over:str,
                                    scheme:str,
                                    dpol2add:list,
                                    path_table:str,
                                    verbose:bool,
                                    )-> dict:
    """Extract min/step/max parameters and necessary columns in the table for later
    computation of the dual-pol variables.
    
    Args:
        band (str): radar band.
        hydrometeors (list): list of hydrometeors for which tables must be read.
        scheme (str): microphysics scheme.
        path_table (str): table's directory path.

    Returns:
        table_dict (dict) : dictionary containing min/step/max values for multiple parameters
    """
    print("Reading tables for",band,"band")
    deb_timer = tm.time()
    micro_for_table = scheme[0:4]
    table_dict, parameters_to_retrieve, columns_to_retrieve = initialize_table_dictionary(dpol2add=dpol2add)
    
    for h in hydrometeors:
        if h == 'cc':
            hfile = cloud_water_species(hydrometeor=h, cloud_water_over=cloud_water_over)
        else :
            hfile = h
        nomfileCoefInt = f'{path_table}TmatCoefInt_{micro_for_table}_{moments[h]}M_{band}{hfile}'
        
        if verbose : print("\tReading min/step/max for",h)
        df_params = pd.read_csv(nomfileCoefInt, sep=";",nrows = 1)
        for value in parameters_to_retrieve :
            table_dict[value][h] = np.copy(df_params[value])[0]
        del df_params
        
        if verbose : print("\tRetrieving necessary columns in the table for",h)
        df_columns = pd.read_csv(nomfileCoefInt, sep=";",skiprows = [0, 1])
        for columnName in columns_to_retrieve :
            if columnName == 'Fw' or columnName == 'N' :
                if moments[h]==1 :
                    table_dict['Fw'][h] = df_columns['P3'].to_numpy()
                elif moments[h]==2 :
                    table_dict['N'][h] = df_columns['P3'].to_numpy()
            else :
                table_dict[columnName][h] = df_columns[columnName].to_numpy()
        del df_columns
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")
    return table_dict

    

def perform_nD_interpolation(tableDict:dict, which_columns:list, hydrometeor:str, colName:str,
                 colMin:float, colStep:float, colMax:float,
                 el_temp:np.ndarray, Tc_temp:np.ndarray, colTable:np.ndarray,
                 M_temp:np.ndarray) -> dict :
    """Construct 3D fields of scattering coefficients, based on the tables. The interpolation is
    currently performed with 4 fields (elevation, temperature, P3 (=Fw or Nc), and content).

    Args:
        tableDict (dict): dictionnary containing the necessary columns of the table
        which_columns (list): list of the columns name on which to perform the interpolation
        hydrometeor (str): hydrometeor type
        colName (str): P3 name (either Fw or Nc)
        colMin (float): P3 min (either Fwmin or Ncmin)
        colStep (float): P3 step (either Fwstep or Ncstep)
        colMax (float): P3 min (either Fwmax or Ncmax)
        el_temp (np.ndarray): elevation 3D field
        Tc_temp (np.ndarray): temperature 3D field
        colTable (np.ndarray): P3 (=Fw or Nc) 3D field
        M_temp (np.ndarray): content 3D field

    Returns:
        scatCoefsDict (dict): dictionnary containing fields of interpolated scattering coefficients over the grid
    """

    columns_2use_for_interpolation = 4 # elevation, temperature, P3 (=Fw or Nc), and content 
                                      
    # Find position in the table
    [kTmat, ELEVred, Tcred, P3red, Mred] = CALC_KTMAT(ELEV = el_temp,
                                                      Tc = Tc_temp,
                                                      P3r = colTable,
                                                      M = M_temp,
                                                      ELEVmin = tableDict['ELEVmin'][hydrometeor],
                                                      ELEVmax = tableDict['ELEVmax'][hydrometeor],
                                                      ELEVstep = tableDict['ELEVstep'][hydrometeor],
                                                      Tcmin = tableDict['Tcmin'][hydrometeor],
                                                      Tcmax = tableDict['Tcmax'][hydrometeor],
                                                      Tcstep = tableDict['Tcstep'][hydrometeor],
                                                      expMmin = tableDict['expMmin'][hydrometeor],
                                                      expMmax = tableDict['expMmax'][hydrometeor],
                                                      expMstep = tableDict['expMstep'][hydrometeor],
                                                      P3min = colMin,
                                                      P3max = colMax,
                                                      P3step = colStep,
                                                      P3name = colName,
                                                      ncol_interpolation = columns_2use_for_interpolation,
                                                      shutdown_warnings = True,
                                                      )
    
     # Store scat coef values for each min/max born in Matcoef     
    MatCoef = {}
    idx_key_pair = []
    
    for idx,key in enumerate(which_columns) :
        idx_key_pair += [(idx,key)] 
        for ind in range((2**columns_2use_for_interpolation)):
            MatCoef[idx, ind] = tableDict[key][hydrometeor][kTmat[ind]]
    
    # Interpol scat coef values
    scatCoefsDict = INTERPOL(ELEVred, Tcred, P3red, Mred, MatCoef,which_columns,idx_key_pair)   

    return scatCoefsDict 



def  CALC_KTMAT(ELEV:float, Tc:np.ndarray, P3r, M:np.ndarray,
                ELEVmin:float, ELEVmax:float, ELEVstep:float,
                Tcmin:float, Tcmax:float, Tcstep:float,
                P3min:float, P3max:float, P3step:float,
                expMmin:float, expMstep:float, expMmax:float,
                P3name:str, ncol_interpolation:int, shutdown_warnings:bool) :
    """Return the indexes of the scattering coef corresponding to the upper and lower bounds
     of ELEV, Tc, P3, M for => used for the interpolation of these coefficients.
       
    Args:
        ELEV (float) : elevation in radians
        Tc (np.ndarray) : temperature in Celsius degree
        P3r : liquid water fraction (1-moment) or concentration (2-moment)
        M (np.ndarray) : contents in kg/m3
        ELEVmin,max,step (float) : degrees
        Tcmin,max,step (float) : Celsius degree
        P3min,max,step (float) : values between 0 and 1
        expMmin,step,max (float) : exponential of the concentration
        P3name : "Fw" liquid water fraction or "Nc" number concentration

    Returns:
        kTmat (dict) : dict with 2**ncol_interpolation arrays with the M shape/size
    """

    kTmat={}
    
    for n in range(2**ncol_interpolation):
        kTmat[n]=np.copy(M)*float('nan') 
    
    #Hydromet content
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
            print("Error P3name (3d parameter in the table), the only options are Fw (liq water fraction) or Nc (number concentration)")
            sys.exit(1)
    
    # If ELEV, Tc, P3 or M are outside min and max ranges:
    # warning and the values are set to the min (if below min) or max (if over max)    
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
    
    condok=np.where((ELEV>=ELEVmin) & (ELEV<=ELEVmax) 
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
    ELEVred=np.copy(ELEV)*0.0
    Tcred=np.copy(Tc)*0.0
    P3red=np.copy(P3)*0.0
    Mred=np.copy(M)*0.0
    
    msk=ELEVsup>ELEVinf
    ELEVred[msk]=(ELEV[msk]-ELEVinf[msk])/(ELEVsup[msk]-ELEVinf[msk])
    
    msk=Tcsup>Tcinf
    Tcred[msk]=(Tc[msk]-Tcinf[msk])/(Tcsup[msk]-Tcinf[msk])
    
    msk=P3sup>P3inf
    P3red[msk]=(P3[msk]-P3inf[msk])/(P3sup[msk]-P3inf[msk])
    
    msk=Msup>Minf
    Mred[msk]=(M[msk]-Minf[msk])/(Msup[msk]-Minf[msk])
        
    # On regroupe autant que possible les produits et sommes dans le calcul de kTmat.
    fact1=nP3*nM
    fact2=nELEV*fact1
    fact51=kTc*fact2
    fact52=kTcs*fact2
    fact61=kELEV*fact1
    fact62=kELEVs*fact1
    fact71=kP3*nM
    fact72=kP3s*nM
    
    som11=fact51
    som111=som11+fact61
    som112=som11+fact62
    som1111=som111+fact71
    som1112=som111+fact72
    som1121=som112+fact71
    som1122=som112+fact72
    som12=fact52
    som121=som12+fact61
    som122=som12+fact62
    som1211=som121+fact71
    som1212=som121+fact72
    som1221=som122+fact71
    som1222=som122+fact72
    
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
    
    
    for i in range(2**ncol_interpolation):
        kTmat[i][np.isnan(kTmat[i])]=-999
        kTmat[i]=kTmat[i].astype(int)
 
    return kTmat, ELEVred,Tcred,P3red,Mred



def  INTERPOL(ELEVred:np.ndarray, Tcred:np.ndarray, Fwred:np.ndarray, Mred:np.ndarray,
              MatCoef:dict, which_columns:list, idx_key_pair:list )-> dict :
    """Multidimensional interpolation with outputs from CALC_KTMAT function

    Args:
        ELEVred (np.ndarray): elevation
        Tcred (np.ndarray): temperature
        Fwred (np.ndarray): liquid water fraction
        Mred (np.ndarray): content
        MatCoef (np.ndarray): ncoefx2**4 matrix with max ncoef=7 table's coefficient columns to interpolate over the 2**4 bounds. 
    
    Returns:
        interpolated_fields (dict): dict containing the interpolated values
    """
    ncoef=len(which_columns)
    nval=Mred.shape[0]
    print("\t  INTERPOL func over nval=",nval)
    VectCoef={}

    #--- Linear interpolation ---
    for indcoef in range(ncoef):
        VectCoef[indcoef]=np.empty((nval))
        VectCoef[indcoef]= \
           ((1-Tcred)*(
                (1-ELEVred)*(
                    (1-Fwred)*(
                        (1-Mred)*MatCoef[indcoef,0] + Mred*MatCoef[indcoef,1] )+
                    Fwred*(
                        (1-Mred)*MatCoef[indcoef,2] + Mred*MatCoef[indcoef,3] ) 
                )+ ELEVred*(
                    (1-Fwred)*(
                        (1-Mred)*MatCoef[indcoef,4] + Mred*MatCoef[indcoef,5] )+
                    Fwred*(
                        (1-Mred)*MatCoef[indcoef,6] + Mred*MatCoef[indcoef,7] ) 
                )
            )+ Tcred*(
                (1-ELEVred)*(
                    (1-Fwred)*(
                        (1-Mred)*MatCoef[indcoef,8] + Mred*MatCoef[indcoef,9] )+
                    Fwred*(
                        (1-Mred)*MatCoef[indcoef,10] + Mred*MatCoef[indcoef,11] )
                )+ ELEVred*(
                    (1-Fwred)*(
                        (1-Mred)*MatCoef[indcoef,12] + Mred*MatCoef[indcoef,13] )+
                    Fwred*(
                        (1-Mred)*MatCoef[indcoef,14] + Mred*MatCoef[indcoef,15] )
                )))
        #print indcoef,np.count_nonzero(~np.isnan(VectCoef[indcoef]))
    
    interpolated_fields = {key:np.copy(VectCoef[idx]) for (idx,key) in idx_key_pair}
              
    del VectCoef
    return interpolated_fields