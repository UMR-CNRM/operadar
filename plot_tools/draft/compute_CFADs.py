import os
import sys
#from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd

sys.path.insert(0, "./lib")
import cfad_lib as cfad
import level2alt as lev2alt


##fileDir   = '/home/davidcl/Programmation/data/final_data/'
#simuList = ["MesoNH_ICE3"] #"Arome_oper","MesoNH_ICE3","MesoNH_LIMA","MesoNH_LIMAAG"]
##simuList=["MesoNH_LIMA"] #"Arome_oper"]
##simuCloe = ["Arome_ICE3","Arome_LIMA"]
#simuCloe =[]

#simu=sys.argv[1]
simu="Arome_ICE3"

dataDir = '/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/'
#dataDir = '/home/augros/DONNEES/MESONH/CORSE/'
statDir = dataDir

fileDirDict   = {'Arome_oper' : dataDir+'AROME/dpolvar/',
             'MesoNH_ICE3': dataDir+'CT1KM/dpolvar/',
                 'MesoNH_LIMA': dataDir+'LIREF/dpolvar',
                 'MesoNH_LIMAAG': dataDir+'LIMAAG/dpolvar/',
                 'obs' : dataDir+'OBS/',
                 'Arome_ICE3':dataDir+"AROME_ICE3/",
                 'Arome_LIMA':dataDir+"AROME_LIMA/",
             }



# === SETTINGS === #
inside = 'zh' #'ZDR_columns' #'cell_core' 'cell_envelop'
threshold = {'zh':40.,
             'ZDR_columns':0.,
             'cell_core':0.,
             'cell_envelop':0.,
        }

listVarPol = ['zh','zdr','kdp','rhohv'] # le nom des variables diffère dans les netcdf de Cloe
listVarHydro = ['rr','ii','ss','gg','wg','vv','cc']
if (simu=="MesoNH_LIMAAG"):
    listVarHydro=['rr','ii','ss','gg','wg','hh','wh','vv','cc']


resolV,alti_min,alti_max=500,0,15e3
lvl_intervals=2
altiList = [int(x) for x in np.arange(alti_min,alti_max,resolV)]

#==== Observations or simu Cloe statistics ==========     
if (simu=="obs" or simu=="Arome_LIMA" or simu=="Arome_ICE3"):
    listVarTot = listVarPol
    print('Extract and compute stats for: ',simu)
    
    
    pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_{simu}.pkl'
    statDict = {'altitude' : [],
                'varname' : [],
                'dataType' : [],
                'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
            }
    for alti in altiList :
        for var in listVarPol:
            print(alti,var)
            tmp_arr = []
            tmp_arr = cfad.extract_values_cloe(fileDirDict[simu],simu,tmp_arr,var,inside,threshold[inside],alti)
            statDict = cfad.compute_stats(simu,statDict,var,alti,tmp_arr)
            del tmp_arr 
    
    df = pd.DataFrame(data=statDict)
    df.to_pickle(pickle_name)        


# ==== Simulation statistics ==========              
listVarTot = listVarPol+listVarHydro

# === AROME oper and MesoNH simulations
#for simu in simuList:
fileDir=fileDirDict[simu]
pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_{simu}.pkl'
print('Extract and compute stats for: ',simu)

statDict = {'altitude' : [],
        'varname' : [],
        'dataType' : [],
        'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
    }

# Interpolation of dpol variables in regular altitude levels
for filename in sorted(os.listdir(fileDir)):
    print(filename)
    f = os.path.join(fileDir, filename)
    if '0000.nc' in f: 
        fout=f[:-3] + "interp.nc"
        if os.path.isfile(fout):
            print("file ",fout," exists, no interpolation needed")
        else:   
            ds = xr.open_dataset(f)
            # MesoNH: select one out of 2 levels from 0 to 85 (Alt level 85 ~ 17 km)
            if ("MesoNH" in filename):
                sub_ds = ds.sel(level=slice(0,85,lvl_intervals))
            # Arome: 7 to 89 (Alt level 7 ~ 16 km)
            if ("Arome" in filename):
                sub_ds = ds.sel(level=slice(7,89,lvl_intervals))
            ds_interp=lev2alt.interpolate_dataset(sub_ds,resolV,alti_min,alti_max)
            ds_interp.to_netcdf(fout)

# Extraction of values to build the CFAD
for alti in altiList :
    for var in listVarTot:
        print(alti,var)
        tmp_arr = []
        tmp_arr = cfad.extract_values_cloe(fileDirDict[simu],simu,tmp_arr,var,inside,threshold[inside],alti)
        statDict = cfad.compute_stats(simu,statDict,var,alti,tmp_arr)
        del tmp_arr               

df = pd.DataFrame(data=statDict)
df.to_pickle(pickle_name)
    
## === AROME Cloe's simulations
#for simu in simuCloe:
#    pickle_name = f'{statDir}df_for_stats_{inside}sup{str(int(threshold[inside]))}_{simu}.pkl'
#    print('Extract and compute stats for: ',simu)
#    
#    statDict = {'altitude' : [],
#            'varname' : [],
#            'dataType' : [],
#            'Q5':[], 'Q25':[], 'Q50':[] , 'Q75':[] , 'Q95':[],
#        }
#    
#    for alti in altiList :
#        for var in listVarTot:
#            print(alti,var)
#            tmp_arr = []
#            tmp_arr = cfad.extract_values_cloe(fileDirDict[simu],simu,tmp_arr,var,inside,threshold[inside],alti)
#            statDict = cfad.compute_stats(simu,statDict,var,alti,tmp_arr)
#            del tmp_arr               
#
#    df = pd.DataFrame(data=statDict)
#    df.to_pickle(pickle_name)

