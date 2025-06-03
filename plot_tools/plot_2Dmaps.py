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


# ==== Functions ======
#plot_map(var,var_lev.lon, var_lev.lat, var_lev, bounds,cmap,lev, fileDir, filename):

def plot_map(var,lon,lat,data,bounds,cmap,lev,fileDir,filename,imgDir,simu,day, time,lon_min,lon_max,lat_min,lat_max):
    """
    Tracer une carte avec les données 2D.

    Args:
        var (str): Le nom de la variable radar.
        lon (numpy.array): Les valeurs de longitude.
        lat (numpy.array): Les valeurs de latitude.
        data (numpy.array): Les données 2D.
        bounds (numpy.array): Les limites des intervalles de couleur.
        cmap (str ou Colormap): La colormap à utiliser.
        lev (int): Le niveau de la variable radar. Si level=-1, on affiche ma valeur max sur tous les niveaux.
        fileDir (str): Le répertoire où sauvegarder l'image.
        filename (str): Le nom du fichier pour sauvegarder l'image.
        lon_min (float): Valeur minimale de longitude du domaine.
        lon_max (float): Valeur maximale de longitude du domaine.
        lat_min (float): Valeur minimale de latitude du domaine.
        lat_max (float): Valeur maximale de latitude du domaine.
    """

    plt.figure(figsize=(10,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    contour_plot = ax.contourf(lon, lat, data,levels=bounds,cmap=cmap, extend='both')
    
    # Add colorbar with label
    cbar = plt.colorbar(contour_plot,orientation='vertical', ticks=bounds, shrink=0.7)
    cbar.set_label(var+" ("+varDict[var]['unit']+")")
    
    # Définir les limites du domaine
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Add coastlines and gridlines
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax.add_feature(cfeature.STATES.with_scale('10m'), linestyle='-', edgecolor='black')
    
    # Afficher les labels des latitudes et des longitudes uniquement à gauche et en bas
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
    # Add title and labels
    varTitle=varDict[var]['name'] +' at level '+str(lev)
    filepng=simu+'_'+day+time+f'_{var}_lev{lev}'
    
    if (lev==-1):
        varTitle=varDict[var]['name']+" max"
        filepng=simu+'_'+day+time+f'_{var}max'
    
    title=varTitle+" - "+simu+" - "+time+" UTC"    
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    
    plt.savefig(imgDir+'/'+filepng+'.png',bbox_inches='tight',pad_inches=0,dpi=100)
    plt.close()
# ===========================================================================


# === Declarations ====
    
interpole=False #True
plot_nointerpole=True    
 
band="C"
day="20220818"
lat_min,lat_max,lon_min,lon_max=41.,45.,4.,11.

simuList = ["Arome_oper","MesoNH_ICE3","MesoNH_LIMA","MesoNH_LIMAAG"]
simuList=["Arome_oper"]
#simuCloe = ["Arome_ICE3","Arome_LIMA","obs"]
simuCloe =["obs"]
simuDiag = [] #"MesoNH_Rayleigh"]

#timeList=["0100","0200","0300","0400","0500","0600","0700","0800","0900","1000"]
timeList=["0800"]

dataDir = '/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/'
imgDir = dataDir+"IMG/"
#dataDir = '/home/augros/DONNEES/MESONH/CORSE/'
statDir = dataDir

fileDirDict   = {'Arome_oper' : dataDir+'AROME/dpolvar/',
             'MesoNH_ICE3': dataDir+'CT1KM/dpolvar/',
                 'MesoNH_LIMA': dataDir+'LIREF/dpolvar/',
                 'MesoNH_LIMAAG': dataDir+'LIMAAG/dpolvar/',
                 'obs' : dataDir+'OBS/',
                 'Arome_ICE3':dataDir+"AROME_ICE3/",
                 'Arome_LIMA':dataDir+"AROME_LIMA/",
                 'MesoNH_Rayleigh':"/scratch/work/straussc/CORSE/005_run1/REF/",
             }

varDict = {
    'zh': {'name': 'Reflectivity', 'min': 8, 'max': 64, 'step':4,'unit': 'dBZ'},
    'zdr': {'name': 'Differential Reflectivity', 'min': 0, 'max': 8, 'step': 1, 'unit': 'dB'},
    'kdp': {'name': 'Specific Differential Phase', 'min': 0, 'max': 10, 'step':1, 'unit': '°/km'}
}


# options for MesoNH diag file (from operadar config files: conf_MesoNH_ICE3_CORSEbe.py et CORSE_MesoNG.csv)
step_seconds = 900.
commonFilename = "CT1KM.1.SEG01."
run="20220818000000"
deb="20220818060000"
fin="20220818060000"
start = datetime.datetime.strptime(deb, "%Y%m%d%H%M%S")
end = datetime.datetime.strptime(fin, "%Y%m%d%H%M%S")

datetimelist = []
current_date = start
while current_date <= end:
    datetimelist.append(current_date)
    current_date += datetime.timedelta(seconds=step_seconds)

print([date.strftime("%d-%m-%Y %H:%M:%S") for date in datetimelist])



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
altiList = [1000] #1000,2000] #2000] #[int(x) for x in np.arange(0,15e3,500)]
levelList=[] #20] #[89] #,87,85,80,75,70,60,50,40,30,20,10,0]
listVarPol = ['zh','zdr','kdp'] #'zdr','kdp','rhohv'] 
listVarPol = ['zh'] 
listVarHydro = [] #['rr','ii','ss','gg','wg','vv','cc']

listVarTot = listVarPol+listVarHydro

             
# ======= Cloe nc files after interpolation over altitude z ===
for simu in simuCloe:
    print('Plot 2D maps for: ',simu)
    micro=simu.split("_")[-1]
    
    fileDir=fileDirDict[simu]
    filename=day+"_AJAC-COLL-NIME_"+micro+".nc"
    f = os.path.join(fileDir, filename)
    ds = xr.open_dataset(f)
    #ds=ds.rename({'zh':'Zh','zdr':'Zdr','kdp':'Kdp','rhohv':'Rhohv'})
        
    # loop time
    day_dt = datetime.datetime.strptime(day, "%Y%m%d")
    for time in timeList:
        time_dt = datetime.datetime.strptime(time, "%H%M")
        datetime_combined = datetime.datetime.combine(day_dt.date(), time_dt.time())
        time_str = datetime_combined.strftime("%Y-%m-%dT%H:%M:%S")
    
        # loop variable
        for var in listVarTot:                
            vmin,vmax,step=varDict[var]['min'],varDict[var]['max'],varDict[var]['step']
            bounds = np.arange(vmin, vmax + step, step)  
            # loop alti
            for alti in altiList :
                data=ds[var].sel(z=alti).sel(time=time_str)
                plot_map(var,data.lon, data.lat, data, bounds,cmap,alti,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
                
#            # Max over all levels
            print("Plot ",var," max")    
            datamax = ds[var].sel(time=time_str).max(dim='z')
            plot_map(var,datamax.lon, datamax.lat, datamax, bounds,cmap,-1,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
            
    ds.close() ; del ds
        
# ======= AROME oper and MesoNH operadar nc files (no interpolation) ===
if (plot_nointerpole):
    for simu in simuList:
        print('Plot 2D maps for: ',simu)
        fileDir=fileDirDict[simu]
        dataType=simu
        if (simu=='Arome_oper'):
            dataType='Arome_ICE3'
         
        for time in timeList:
            
            filename="dpolvar_"+dataType+"_"+band+"_"+day+time+"00.nc"
            
            if not os.path.isfile(fileDir+filename):
                old_name = filename[:-5] + ".nc"
                os.rename(fileDir+old_name, fileDir+filename)
                print(f"Le fichier a été renommé en : {filename}")
            
            if os.path.isfile(fileDir+filename):
                print('  ',filename)
                f = os.path.join(fileDir, filename)
                ds = xr.open_dataset(f)  
                
                if (interpole):
                    ds_interp=lev2alt.interpolate_dataset(ds,lvl_intervals=1,resolV=500,alti_min=0,alti_max=15e3)
                    ds_interp.to_netcdf(fileDir+dataType+"_"+band+"_"+day+time+"00interp.nc")
                
                for var in listVarTot:                
                    vmin,vmax,step=varDict[var]['min'],varDict[var]['max'],varDict[var]['step']
                    bounds = np.arange(vmin, vmax + step, step)           
                    
                    # Levels
                    for lev in levelList :
                        print("Plot ",var," at level ",lev)            
                        
                        data = ds[var].sel(level=lev)                    
                        plot_map(var,data.lon, data.lat, data, bounds,cmap,lev,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
                        
                    # Max over all levels
                    print("Plot ",var," max")    
                    datamax = ds[var].max(dim='level')
                    plot_map(var,datamax.lon, datamax.lat, datamax, bounds,cmap,-1,fileDir,filename,imgDir,simu,day,time,lon_min,lon_max,lat_min,lat_max)
                
                if (interpole):
                    ds_interp.close() ; del ds_interp
                ds.close() ; del ds    
                
# ======= MesoNH DIAG files (RARE with Rayleigh scattering) ============
for simu in simuDiag:
    print('Plot 2D maps for: ',simu)
    fileDir=fileDirDict[simu]  
    datetime_run=datetime.datetime.strptime(run,'%Y%m%d%H%M%S')
    
    for date in datetimelist: 
        print(date.strftime("%Y%m%d %H%M%S"))
        heure=date.strftime("%H%M")
    
        model_ech=((date - datetime_run).total_seconds())/step_seconds
        ech=(str(int(model_ech))).zfill(3)
        print(ech)
        modelfile=fileDir+commonFilename+ech+"DIA.nc"
        ncfile = Dataset(modelfile,'r')
        #print(ncfile.variables.keys())

        time=ncfile.variables['time'][:]
        lat = ncfile.variables['latitude'][:,0]
        lon = ncfile.variables['longitude'][0,:]
        rare = ncfile.variables['RARE'][0,:,:,:]
        rare[rare>998.]=np.nan
        
        # Max Zh over all levels
        print("Plot Zh Rayleigh RARE max")    
        var="Zh"
        vmin,vmax,step=varDict[var]['min'],varDict[var]['max'],varDict[var]['step']
        bounds = np.arange(vmin, vmax + step, step)
        datamax = np.nanmax(rare,axis=0)
        plot_map(var,lon, lat, datamax, bounds,cmap,-1,fileDir,modelfile,imgDir,simu,day,heure,lon_min,lon_max,lat_min,lat_max)

