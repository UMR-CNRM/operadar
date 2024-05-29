import os
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd

#fileDir   = '/home/davidcl/Programmation/data/final_data/'
simuList = ["Arome_oper","MesoNH_ICE3","MesoNH_LIMA","MesoNH_LIMAAG"]
#simuList=[]
simuCloe = ["Arome_ICE3","Arome_LIMA"]

computeObs=True

dataDir = '/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/'
#dataDir = '/home/augros/DONNEES/MESONH/CORSE/'
statDir = dataDir

fileDirDict   = {'Arome_oper' : dataDir+'AROME_oper/dpolvar/',
             'MesoNH_ICE3': dataDir+'CT1KM/dpolvar/',
                 'MesoNH_LIMA': dataDir+'LIREF/dpolvar',
                 'MesoNH_LIMAAG': dataDir+'LIMAH/dpolvar/',
                 'obs' : dataDir+'OBS/',
                 'Arome_ICE3':dataDir+"AROME_ICE3/",
                 'Arome_LIMA':dataDir+"AROME_LIMA/",
             }



# === SETTINGS === #
inside = 'Zh' #'ZDR_columns' #'cell_core' 'cell_envelop'
threshold = {'Zh':0.,
             'ZDR_columns':0.,
             'cell_core':0.,
             'cell_envelop':0.,
        }
altiList = [int(x) for x in np.arange(0,15e3,500)]
listVarPol = ['Zh','Zdr','Kdp','Rhohv'] 
#listVarObs = ['zh','zdr','kdp','rhohv'] ] # le nom des variables diffère dans les netcdf de Cloe
listVarHydro = [] #['rr','ii','ss','gg','wg','vv','cc']


# === FUNCTIONS === #
def compute_stats(dataType:str,statDict,var,alti,valueList): 
    statDict['varname'] += [var]
    statDict['altitude'] += [alti]
    statDict['dataType'] += [dataType]
    try :
        statDict['Q5'] += [np.percentile(valueList,5)]
        statDict['Q25'] += [np.percentile(valueList,25)]
        statDict['Q50'] += [np.percentile(valueList,50)]
        statDict['Q75'] += [np.percentile(valueList,75)]
        statDict['Q95'] += [np.percentile(valueList,95)]
    except :
        statDict['Q5'] += [np.nan]
        statDict['Q25'] += [np.nan]
        statDict['Q50'] += [np.nan]
        statDict['Q75'] += [np.nan]
        statDict['Q95'] += [np.nan]
    return statDict
     
def extract_values(fileDir,dataType:str,tmp_array,var,alti):
    for filename in sorted(os.listdir(fileDir)):
        #print(filename)
        f = os.path.join(fileDir, filename)
        #if dataType in f :
        if '.nc' in f:    
            #print('  ',filename)
            ds = xr.open_dataset(f)
            #ds
            if np.isin(var,listVarHydro):
                arr = ds.M.sel(hydrometeor=var,Alt=alti).where(ds[inside] > threshold[inside]).values
            else :
                #arr = ds[var].sel(alti-250.<=ds["Alt"]<alti+250.).where(ds[inside] > threshold[inside]).values
                mask = (ds["Alt"]>=alti-250.) & (ds["Alt"]< alti+250.) & (ds[inside] > threshold[inside])
                arr = ds[var].where(mask).values
            tmp_array += arr[~np.isnan(arr)].tolist()
            ds.close() ; del ds
    return tmp_array

def extract_values_cloe(fileDir,dataType:str,tmp_array,var,alti):
    for filename in sorted(os.listdir(fileDir)):
        f = os.path.join(fileDir, filename)
        #if f'{dataType}.nc' in f :
        if '.nc' in f:    
            print('  ',filename)
            ds = xr.open_dataset(f)
            ds=ds.rename({'zh':'Zh','zdr':'Zdr','kdp':'Kdp','rhohv':'Rhohv'})
            if np.isin(var,listVarHydro):
                arr = ds.M.sel(hydrometeor=var,z=alti).where(ds[inside] > 0).values
            else :
                arr = ds[var].sel(z=alti).where(ds[inside] > 0).values
            tmp_array += arr[~np.isnan(arr)].tolist()
            ds.close() ; del ds
    return tmp_array


#==== Observations statistics ==========     
if (computeObs):
    
    listVarTot = listVarPol
    print('Extract and compute stats for: observations')
    
    
    pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_obs.pkl'
    statDict = {'altitude' : [],
                'varname' : [],
                'dataType' : [],
                'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
            }
    for alti in altiList :
        for var in listVarPol:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values_cloe(fileDirDict['obs'],'obs',tmp_arr,var,alti)
            statDict = compute_stats('obs',statDict,var,alti,tmp_arr)
            del tmp_arr 
    
    df = pd.DataFrame(data=statDict)
    df.to_pickle(pickle_name)        


# ==== Simulation statistics ==========              
listVarTot = listVarPol+listVarHydro

# === AROME oper and MesoNH simulations
for simu in simuList:
    pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_{simu}.pkl'
    print('Extract and compute stats for: ',simu)
    
    statDict = {'altitude' : [],
            'varname' : [],
            'dataType' : [],
            'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
        }
    
    for alti in altiList :
        for var in listVarTot:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values(fileDirDict[simu],simu,tmp_arr,var,alti)
            statDict = compute_stats(simu,statDict,var,alti,tmp_arr)
            del tmp_arr               

    df = pd.DataFrame(data=statDict)
    df.to_pickle(pickle_name)
    
# === AROME Cloe's simulations
for simu in simuCloe:
    pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_{simu}.pkl'
    print('Extract and compute stats for: ',simu)
    
    statDict = {'altitude' : [],
            'varname' : [],
            'dataType' : [],
            'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
        }
    
    for alti in altiList :
        for var in listVarTot:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values_cloe(fileDirDict[simu],simu,tmp_arr,var,alti)
            statDict = compute_stats(simu,statDict,var,alti,tmp_arr)
            del tmp_arr               

    df = pd.DataFrame(data=statDict)
    df.to_pickle(pickle_name)

