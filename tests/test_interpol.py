from operadar.radar.dualpol_variables import hypercube_interpolation
from operadar.read.lookup_tables import read_and_extract_tables_content, retrieve_needed_columns


import numpy as np
import sys


band='C'
hydrometeor='cc'
micro='ICE3'
moment=1
dpol2add=['Zh']
path_table ='/home/davidcl/Programmation/operadar/tables_generator/tables/default/'
tableDict=read_and_extract_tables_content(band=band,
                                          hydrometeors=[hydrometeor],
                                          moments={hydrometeor:moment},
                                          scheme=micro,
                                          dpol2add=dpol2add,
                                          path_table=path_table,
                                          test_interpolation=True,
                                          )

colName='Fw'
#print(tableDict)

def perform_nD_interpolation(tableDict,hydrometeor,colName,elev,T,P3,content,dpol2add):
    [kTmat, ELEVred, Tcred, P3red, Mred] = CALC_KTMAT(ELEV = elev,
                                                      Tc = T,
                                                      P3r = P3,
                                                      M = content,
                                                      ELEVmin = tableDict['ELEVmin'][hydrometeor],
                                                      ELEVmax = tableDict['ELEVmax'][hydrometeor],
                                                      ELEVstep = tableDict['ELEVstep'][hydrometeor],
                                                      Tcmin = tableDict['Tcmin'][hydrometeor],
                                                      Tcmax = tableDict['Tcmax'][hydrometeor],
                                                      Tcstep = tableDict['Tcstep'][hydrometeor],
                                                      expMmin = tableDict['expMmin'][hydrometeor],
                                                      expMmax = tableDict['expMmax'][hydrometeor],
                                                      expMstep = tableDict['expMstep'][hydrometeor],
                                                      P3min = tableDict[f'{colName}min'][hydrometeor],
                                                      P3max = tableDict[f'{colName}max'][hydrometeor],
                                                      P3step = tableDict[f'{colName}step'][hydrometeor],
                                                      P3name = colName,
                                                      ncol_interpolation = 4,
                                                      shutdown_warnings = False,
                                                      )
    # Store scat coef values for each min/max born in Matcoef
    MatCoef = {}
    idx_key_pair = []
    scatCoef_columns = retrieve_needed_columns(dpol2add=dpol2add)
    for idx,key in enumerate(scatCoef_columns) :
        idx_key_pair += [(idx,key)] 
        for ind in range((2**4)):
            MatCoef[idx, ind] = tableDict[key][hydrometeor][kTmat[ind]]
            print('index',kTmat[ind],'---> borne T=',tableDict['Tc'][hydrometeor][kTmat[ind]],
                    ' /  borne M=',tableDict['M'][hydrometeor][kTmat[ind]],
                    ' /  borne ELEV=',tableDict['ELEV'][hydrometeor][kTmat[ind]],
                    ' /  borne P3=',tableDict[colName][hydrometeor][kTmat[ind]],
                    ' /  sighh =',tableDict['sighh'][hydrometeor][kTmat[ind]]) 
    # Interpol scat coef values
    scatCoefsDict = INTERPOL(ELEVred, Tcred, P3red, Mred, MatCoef,scatCoef_columns,idx_key_pair)   

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
    #print("\t  INTERPOL func over nval=",nval)
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



# Fake data                       
elev=np.array([0])
T=np.array([-10.5])
P3=np.array([0.])
content=np.array([0.1*1e-3])
print('======================================')
print(f'Point(s) to interpolate for {hydrometeor} (T,M,ELEV,P3):',T,content,elev,P3,)
print('======================================','\n')

print('--------------- KTMAT + INTERPOL ---------------')
scatCoefsDict2=perform_nD_interpolation(tableDict,hydrometeor,colName,elev,T,P3,content,dpol2add)
print('Final value (old way):',scatCoefsDict2,'\n')

print('--------------- NEW HYPERCUBE INTERPOL ---------------')
scatCoefsDict=hypercube_interpolation(tableDict,hydrometeor,colName,elev,T,P3,content,dpol2add,test_mode=True)
print('Final value (hypercube):',scatCoefsDict,'\n')