import os
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl ; mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as datetime
from netCDF4 import Dataset

import epygram

sys.path.insert(0, "./lib")


# === Declarations ====
    
interpole=False #True
plot_nointerpole=False    
 
band="C"
day="20220818"
lat_min,lat_max,lon_min,lon_max=41.,45.,4.,11.

simuList = ["Arome_ICE3","MesoNH_ICE3","MesoNH_LIMA","MesoNH_LIMAAG"]
simuList=["obs"]
#simuCloe = ["Arome_ICE3","Arome_LIMA","obs"]
#simuCloe =["obs"]
simuDiag = [] #"MesoNH_Rayleigh"]

timeList=["0300","0400","0500","0600","0700","0800"]
#timeList=["0500","0600","0700"]

dataDir = '/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/'
imgDir = dataDir+"IMG/"
#dataDir = '/home/augros/DONNEES/MESONH/CORSE/'
statDir = dataDir

fileDirDict   = {'Arome_oper' : dataDir+'AROME_oper/dpolvar/',
             'MesoNH_ICE3': dataDir+'CT1KM/dpolvar/',
                 'MesoNH_LIMA': dataDir+'LIREF/dpolvar/',
                 'MesoNH_LIMAAG': dataDir+'LIMAAG/dpolvar/',
                 'obs' : dataDir+'OBS/',
                 'Arome_LIMA':dataDir+"AROME_LIMA/",
                 'MesoNH_Rayleigh':"/scratch/work/straussc/CORSE/005_run1/REF/",
                 'Arome_ICE3':dataDir+"AROME_ICE3/",
             }

varDict = {
    'zh': {'name': 'Reflectivity', 'min': 40, 'max': 64, 'step':4,'unit': 'dBZ'},
    'zdr': {'name': 'Differential Reflectivity', 'min': 1, 'max': 6.5, 'step': 0.5, 'unit': 'dB'},
    'kdp': {'name': 'Specific Differential Phase', 'min': 0, 'max': 7, 'step':0.5, 'unit': '°/km'},
    'rhohv': {'name': 'Copolar cross correlation coefficient', 'min': 0.9, 'max': 1, 'step':0.01, 'unit': '/'}
}




# ============================================================================
#               Start programm
# ============================================================================


# ======= Colormap ======
# Colormap epygram radar (comme script M. Mandement)
epygram.init_env()
epygram.util.load_cmap('radar')

cmap=plt.get_cmap('radar');cmap.set_under('white') ; cmap.set_over('deeppink')
#cmap=plt.get_cmap('gist_ncar');cmap.set_under('white') ; cmap.set_over('deeppink')


# === SETTINGS === #

altiList=[1000,2000,3000,4000,5000,6000]
#altiList=["max"]

listVarPol = ['zh','zdr','kdp'] #'zdr','kdp','rhohv'] 
listVarPol = ['rhohv'] 
listVarHydro = [] #['rr','ii','ss','gg','wg','vv','cc']

listVarTot = listVarPol+listVarHydro

             
# ======= Cloe nc files after interpolation over altitude z ===
for simu in simuList:
    print('Plot 2D maps for: ',simu)
    micro=simu.split("_")[-1]
    model=simu.split("_")[0]
    if (simu=="Arome_oper"):
        micro="ICE3"
    
    fileDir=fileDirDict[simu]
    if (simu=="obs" or simu=="Arome_LIMA" or simu=="Arome_ICE3"):
        filename=day+"_AJAC-COLL-NIME_"+micro+".nc"
        print(fileDir,filename)
        f = os.path.join(fileDir, filename)
        dsobs = xr.open_dataset(f)
    #ds=ds.rename({'zh':'Zh','zdr':'Zdr','kdp':'Kdp','rhohv':'Rhohv'})

    # loop variable
    for var in listVarTot:                
        vmin,vmax,step=varDict[var]['min'],varDict[var]['max'],varDict[var]['step']
        bounds = np.arange(vmin, vmax + step, step)  
        
        # loop over altitudes
        for alti in altiList:
        
            plt.figure(figsize=(10,8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            
            # loop over time
            day_dt = datetime.datetime.strptime(day, "%Y%m%d")
    
            for time in timeList:
                time_dt = datetime.datetime.strptime(time, "%H%M")
                datetime_combined = datetime.datetime.combine(day_dt.date(), time_dt.time())
                time_str = datetime_combined.strftime("%Y-%m-%dT%H:%M:%S")
                
                if (simu=="obs" or simu=="Arome_LIMA" or simu=="Arome_ICE3"):  
                    if (alti=="max"):
                        data=dsobs[var].sel(time=time_str).max(dim='z')
                        datazh=dsobs["zh"].sel(time=time_str).max(dim='z')
                    else:
                        data=dsobs[var].sel(z=alti).sel(time=time_str)
                        datazh=dsobs["zh"].sel(z=alti).sel(time=time_str)
                else:
                    filename="dpolvar_"+model+"_"+micro+"_"+band+"_"+day+time+"00interp.nc"
                    print(fileDir+filename)
                    if os.path.isfile(fileDir+filename):
                        print('  ',filename,' ok')
                        f = os.path.join(fileDir, filename)
                        ds = xr.open_dataset(f) 

                        if (alti=="max"):
                            data=ds[var].max(dim='z')
                            datazh=ds["zh"].max(dim='z')
                        else:
                            data=ds[var].sel(z=alti)
                            datazh=ds["zh"].sel(z=alti)
                        #end if altimax
                   # else:
		   # 	print('missing file: ',filename)
                    #end if file available
                #end if obs                

                if (time==timeList[0]):
                    datamax=np.copy(data)
                else:
                    datamax=np.fmax(data,datamax)
     
                # reflectivity contours
                contour_plot = ax.contour(datazh.lon, datazh.lat, datazh,levels=[40.],colors='black', linewidths=1)
 
            #end loop over time
            
            # Plot zh, zdr or kdp field with contourf
            contour_plot = ax.contourf(data.lon, data.lat, datamax,levels=bounds,cmap=cmap, extend='both')
            
            # Add colorbar with label
            cbar = plt.colorbar(contour_plot,orientation='vertical', ticks=bounds, shrink=0.7)
            cbar.set_label(var+" ("+varDict[var]['unit']+")")
            
            # Définir les limites du domaine
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            
            # Add coastlines and gridlines
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
            ax.add_feature(cfeature.STATES.with_scale('10m'), linestyle='-', edgecolor='gray')
            
            # Afficher les labels des latitudes et des longitudes uniquement à gauche et en bas
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            
            # Add title and labels
            varTitle=varDict[var]['name'] +' at '+str(alti)+' m'
            filepng=simu+'_'+day+timeList[0]+"to"+time+f'_{var}_lev{alti}'
            
            if (alti=="max"):
                varTitle=varDict[var]['name']+" max"
                filepng=simu+'_'+day+time+f'_{var}max'
            
            title=varTitle+" - "+simu+" - "+timeList[0]+" to "+ time+" UTC"    
            plt.title(title)
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
        
            
            plt.savefig(imgDir+'/'+filepng+'.png',bbox_inches='tight',pad_inches=0,dpi=100)
            plt.close()
        # end loop over altitudes
    
    # end loop over variables            
    
    if (simu=="obs" or simu=="Arome_LIMA" or simu=="Arome_ICE3"):
        dsobs.close();del dsobs
    else:    
        ds.close() ; del ds
        
## ======= AROME oper and MesoNH operadar nc files (no interpolation) ===
#if (plot_nointerpole):
#    for simu in simuList:
#        print('Plot 2D maps for: ',simu)
#        fileDir=fileDirDict[simu]
#        dataType=simu
#        if (simu=='Arome_oper'):
#            dataType='Arome_ICE3'
#         
#        for time in timeList:
#            
#            filename="dpolvar_"+dataType+"_"+band+"_"+day+time+"00.nc"
#            
#            if not os.path.isfile(fileDir+filename):
#                old_name = filename[:-5] + ".nc"
#                os.rename(fileDir+old_name, fileDir+filename)
#                print(f"Le fichier a été renommé en : {filename}")
#            
#            if os.path.isfile(fileDir+filename):
#                print('  ',filename)
#                f = os.path.join(fileDir, filename)
#                ds = xr.open_dataset(f)  
#                
#                if (interpole):
#                    ds_interp=lev2alt.interpolate_dataset(ds,lvl_intervals=1,resolV=500,alti_min=0,alti_max=15e3)
#                    ds_interp.to_netcdf(fileDir+dataType+"_"+band+"_"+day+time+"00interp.nc")
#                
#                for var in listVarTot:                
#                    vmin,vmax,step=varDict[var]['min'],varDict[var]['max'],varDict[var]['step']
#                    bounds = np.arange(vmin, vmax + step, step)           
#                    
#                    # Levels
#                    for lev in levelList :
#                        print("Plot ",var," at level ",lev)            
#                        
#                        data = ds[var].sel(level=lev)                    
#                        plot_map(var,data.lon, data.lat, data, bounds,cmap,lev,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
#                        
#                    # Max over all levels
#                    print("Plot ",var," max")    
#                    datamax = ds[var].max(dim='level')
#                    plot_map(var,datamax.lon, datamax.lat, datamax, bounds,cmap,-1,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
#                
#                if (interpole):
#                    ds_interp.close() ; del ds_interp
#                ds.close() ; del ds    
                
