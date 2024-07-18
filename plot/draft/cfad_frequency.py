import sys
sys.path.append('/home/davidcl/Programmation/utilPythonFunc/')
from get_data_paths import final_data, out_repo_path, data_repo_path

import os
from pathlib import Path
import dask
import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

normalize_cbar = True
variable='zh'
xlabel = 'Reflectivity (dBZ)'#'Copolar correlation coefficient' 'Differential reflectivity (dB)' 'Specific differential phase (°/km)'
xbins = np.arange(25,80,2) # rhohv(0.6, 1.01, 0.01) zh zdr(0,9,0.25) kdp(0,9,0.25)
xlim = (24.5,80)
height = np.arange(0,15500,500)

fileDir = '/home/davidcl/Programmation/data/final_data/'

def obtain_count_hist(dataType:str) :
    dataList = []
    dsk = []

    if dataType == 'obs' :
        for filename in sorted(os.listdir(fileDir)):
            f = os.path.join(fileDir, filename)
            if f[-6:] == f'obs.nc' :
                ds = xr.open_dataset(f)
                dataList += [ds[variable].where((ds.cell_envelop !=0) & (ds.zh >= 25))]
                ds.close() ; del ds
            else :
                pass
    elif dataType == 'ICE3' :
        for filename in sorted(os.listdir(fileDir)):
            f = os.path.join(fileDir, filename)
            if f[-7:] == 'ICE3.nc' :
                ds = xr.open_dataset(f)
                dataList += [ds[variable].where((ds.cell_envelop !=0) & (ds.zh >= 25))]
                ds.close() ; del ds
            else :
                pass
    elif dataType == 'LIMASG' :
        for filename in sorted(os.listdir(fileDir)):
            f = os.path.join(fileDir, filename)
            if f[-9:] == 'LIMASG.nc' :
                ds = xr.open_dataset(f)
                dataList += [ds[variable].where((ds.cell_envelop !=0) & (ds.zh >= 25))]
                ds.close() ; del ds
            else :
                pass

    for j in range(len(height)):
        a = dask.delayed(np.histogram)(dataList[0][:, j], bins=xbins)
        count = a[0]
        for data in dataList[1:] :
            b = dask.delayed(np.histogram)(data[:, j], bins=xbins)
            count += b[0]
        dsk.append(count)

    hist = np.array(dask.compute(*dsk))
    coords = {'x': xbins[:-1], 'z': height}
    dims = ['z', 'x']
    attrs = {'long_name': f'CFAD for {variable}', 'units': '1'}
    data_array = xr.DataArray(data=hist, dims=dims, coords=coords, attrs=attrs)
    data_array['x'].attrs = {'long_name': 'X bins for CFAD', 'units': '1'}
    return data_array

if not Path(f'{data_repo_path}/tmp/cfad_obs_{variable}.nc').exists():
    print(f'Creating {variable} obs dataset')
    cfad_obs = obtain_count_hist('obs')
    cfad_obs.to_netcdf(f'{data_repo_path}/tmp/cfad_obs_{variable}.nc')
else:
    cfad_obs = xr.open_dataarray(f'{data_repo_path}/tmp/cfad_obs_{variable}.nc')

if not Path(f'{data_repo_path}/tmp/cfad_ice3_{variable}.nc').exists():
    print(f'Creating {variable} ICE3 dataset')
    cfad_ice3 = obtain_count_hist('ICE3')
    cfad_ice3.to_netcdf(f'{data_repo_path}/tmp/cfad_ice3_{variable}.nc')
else:
    cfad_ice3 = xr.open_dataarray(f'{data_repo_path}/tmp/cfad_ice3_{variable}.nc')
    
if not Path(f'{data_repo_path}/tmp/cfad_lima_{variable}.nc').exists():
    print(f'Creating {variable} LIMASG dataset')
    cfad_lima = obtain_count_hist('LIMASG')
    cfad_lima.to_netcdf(f'{data_repo_path}/tmp/cfad_lima_{variable}.nc')
else:
    cfad_lima = xr.open_dataarray(f'{data_repo_path}/tmp/cfad_lima_{variable}.nc')


print('Plotting')
fig, axs = plt.subplots(1,3,figsize=(15,6),sharey=True, layout='constrained')
if normalize_cbar :
    # OBS
    axs[0].pcolor(cfad_obs['x'], cfad_obs['z'], cfad_obs,cmap='jet',norm=colors.LogNorm())
    axs[0].set_title('Radar observations',fontsize=14)
    axs[0].set_ylabel("Height (m)",fontsize=12)
    # ICE3
    axs[1].pcolor(cfad_ice3['x'], cfad_ice3['z'], cfad_ice3,cmap='jet',norm=colors.LogNorm())
    axs[1].set_title('AROME + ICE3',fontsize=14)
    # LIMAsg
    cs = axs[2].pcolor(cfad_lima['x'], cfad_lima['z'], cfad_lima,cmap='jet',norm=colors.LogNorm()) 
    axs[2].set_title('AROME + LIMA',fontsize=14)
else :
    # OBS
    axs[0].pcolor(cfad_obs['x'], cfad_obs['z'], cfad_obs,cmap='jet')
    axs[0].set_title('Radar observations',fontsize=14)
    axs[0].set_ylabel("Height (m)",fontsize=12)
    # ICE3
    axs[1].pcolor(cfad_ice3['x'], cfad_ice3['z'], cfad_ice3,cmap='jet')
    axs[1].set_title('AROME + ICE3',fontsize=14)
    # LIMAsg
    cs = axs[2].pcolor(cfad_lima['x'], cfad_lima['z'], cfad_lima,cmap='jet') 
    axs[2].set_title('AROME + LIMA',fontsize=14)

cb = fig.colorbar(cs, label='Count')
fig.supxlabel(xlabel)
for ax in axs :
    ax.set_xlim(xlim)

saveFigPath = '/home/davidcl/Programmation/output/statistiques/'
if normalize_cbar :
    fig.savefig(saveFigPath+f'CFAD_{variable}_all_studyCases_full_cell_normalized_cbar.png', bbox_inches='tight',dpi=100)
else :
    fig.savefig(saveFigPath+f'CFAD_{variable}_all_studyCases_full_cell.png', bbox_inches='tight',dpi=100)
