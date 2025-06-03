import os
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

fileDir   = '/home/davidcl/Programmation/data/final_data/'

statDict = {'altitude' : [],
            'varname' : [],
            'dataType' : [],
            'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
        }

# === SETTINGS === #
inside = 'ZDR_columns' #'cell_core' 'cell_envelop'
altiList = [int(x) for x in np.arange(0,15e3,500)]
listVarPol = ['zh','zdr','kdp','rhohv'] 
listVarHydro = ['rr','ii','ss','gg','wg','vv','cc']
pickle_name = f'{fileDir}df_for_stats_{inside}.pkl'



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
        f = os.path.join(fileDir, filename)
        if f'{dataType}.nc' in f :
            print('  ',filename)
            ds = xr.open_dataset(f)
            if np.isin(var,listVarHydro):
                arr = ds.M.sel(hydrometeor=var,z=alti).where(ds[inside] > 0).values
            else :
                arr = ds[var].sel(z=alti).where(ds[inside] > 0).values
            tmp_array += arr[~np.isnan(arr)].tolist()
            ds.close() ; del ds
    return tmp_array

if not Path(pickle_name).exists():
     
    listVarTot = listVarPol
    for alti in altiList :
        for var in listVarTot:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values(fileDir,'obs',tmp_arr,var,alti)
            statDict = compute_stats('obs',statDict,var,alti,tmp_arr)
            del tmp_arr
              
    listVarTot = listVarPol+listVarHydro
    for alti in altiList :
        for var in listVarTot:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values(fileDir,'ICE3',tmp_arr,var,alti)
            statDict = compute_stats('ICE3',statDict,var,alti,tmp_arr)
            del tmp_arr
            
    for alti in altiList :
        for var in listVarTot:
            print(alti,var)
            tmp_arr = []
            tmp_arr = extract_values(fileDir,'LIMASG',tmp_arr,var,alti)
            statDict = compute_stats('LIMASG',statDict,var,alti,tmp_arr)
            del tmp_arr

    df = pd.DataFrame(data=statDict)
    df.to_pickle(f'{fileDir}stats_{inside}.pkl')

else :
    print('Dataframe for statistics computation already exists')