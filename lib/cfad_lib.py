import os
import numpy as np
import xarray as xr

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
     
def extract_values_cloe(fileDir,dataType:str,tmp_array,var,inside,threshold,alti):
    listVarHydro = ['rr','ii','ss','gg','wg','vv','cc']
    for filename in sorted(os.listdir(fileDir)):
        f = os.path.join(fileDir, filename)
        # For AROME and MesoNH dpol var after interpolation
        if ('interp.nc' in f):    
            print('  ',filename)
            ds = xr.open_dataset(f)
            if np.isin(var,listVarHydro):
                arr = ds.M.sel(hydrometeor=var,z=alti).where(np.nanmax(ds[inside],0) > threshold).values
            else :
                arr = ds[var].sel(z=alti).where(np.nanmax(ds[inside],0) > threshold).values
            tmp_array += arr[~np.isnan(arr)].tolist()
            ds.close() ; del ds    
    
        # For observation file: selection of times with minute=0 in obs file
        if ('obs_hours.nc' in f):   
            print('  ',filename)
            ds = xr.open_dataset(f)
            ds_hours=ds.sel(time=ds.time.dt.minute == 0)
            if np.isin(var,listVarHydro):
                arr = ds_hours.M.sel(hydrometeor=var,z=alti).where(ds_hours["zh_max_z"] > threshold).values
            else :
                arr = ds_hours[var].sel(z=alti).where(ds_hours["zh_max_z"] > threshold).values  
            tmp_array += arr[~np.isnan(arr)].tolist()
            ds.close() ; del ds
        # end obs_hours.nc    
    return tmp_array

