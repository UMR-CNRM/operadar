import os
from pathlib import Path
import dask
import warnings
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

logNorm_cbar = False
compare = True
variable='zh'
xlabel = 'Reflectivity (dBZ)'#'Copolar correlation coefficient' 'Differential reflectivity (dB)' 'Specific differential phase (°/km)'
xbins = np.arange(25,80,2) # rhohv(0.6, 1.01, 0.01) zh zdr(0,9,0.25) kdp(0,9,0.25)
xlim = (24,80)
height = np.arange(0,15500,500)

data_repo_path = '/home/davidcl/Programmation/data/final_data/' #'/cnrm/precip/users/davidcl/final_data/'
saveFigPath = '/home/davidcl/Programmation/output/statistiques/'

def obtain_count_hist(dataType:str) :
    dataList = []
    dsk = []

    if dataType == 'obs' :
        for filename in sorted(os.listdir(data_repo_path)):
            f = os.path.join(data_repo_path, filename)
            if f[-6:] == f'obs.nc' :
                ds = xr.open_dataset(f)
                dataList += [ds[variable].where((ds.cell_envelop !=0) & (ds.zh >= 25))]
                ds.close() ; del ds
            else :
                pass
    elif dataType == 'ICE3' :
        for filename in sorted(os.listdir(data_repo_path)):
            f = os.path.join(data_repo_path, filename)
            if f[-7:] == 'ICE3.nc' :
                ds = xr.open_dataset(f)
                dataList += [ds[variable].where((ds.cell_envelop !=0) & (ds.zh >= 25))]
                ds.close() ; del ds
            else :
                pass
    elif dataType == 'LIMASG' :
        for filename in sorted(os.listdir(data_repo_path)):
            f = os.path.join(data_repo_path, filename)
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

def create_plot(da,ax,logNorm_cbar,title,vmax) :
    if logNorm_cbar : 
        plot = ax.pcolor(da['x'], da['z'], da,cmap='jet',norm=colors.LogNorm(vmin=1,vmax=vmax))
    else :
        plot = ax.pcolor(da['x'], da['z'], da,cmap='jet',vmax=vmax)
    ax.set_title(title,fontsize=14)
    ax.set_ylabel("Height (m)",fontsize=12)
    return plot

def plot_diff(diff,ax,title,vmin,vmax,logNorm_cbar):
    vmax_abs = np.max([np.abs(vmin),np.abs(vmax)])
    if logNorm_cbar : 
        plot = ax.pcolor(diff['x'], diff['z'], diff,cmap='bwr',norm=colors.SymLogNorm(vmin=-vmax_abs,vmax=vmax_abs,linthresh=10,base=10))#vmax=vmax_abs TwoSlopeNorm(vcenter = 0,vmin=vmin,vmax=vmax)
    else :
        plot = ax.pcolor(diff['x'], diff['z'], diff,cmap='seismic',vmax=vmax_abs)
    ax.set_title(title,fontsize=14)
    ax.set_ylabel("Height (m)",fontsize=12)
    return plot

def read_or_compute_cfad_dataset(dataType):
    if not Path(f'{data_repo_path}tmp/cfad_{dataType}_{variable}.nc').exists():
        print(f'Creating {variable} {dataType} dataset')
        cfad_ds = obtain_count_hist(dataType)
        cfad_ds.to_netcdf(f'{data_repo_path}tmp/cfad_{dataType}_{variable}.nc')
    else:
        cfad_ds = xr.open_dataarray(f'{data_repo_path}tmp/cfad_{dataType}_{variable}.nc')
    return cfad_ds

cfad_obs = read_or_compute_cfad_dataset('obs')
cfad_ice3 = read_or_compute_cfad_dataset('ICE3')
cfad_lima = read_or_compute_cfad_dataset('LIMASG')

print('Plotting')
vmax = np.max([cfad_obs.max().values,cfad_ice3.max().values,cfad_lima.max().values])

if compare :
    diff_ice3 = cfad_ice3 - cfad_obs
    diff_lima = cfad_lima - cfad_obs
    vmax_diff = np.max([diff_ice3.max().values,diff_lima.max().values])
    vmin_diff = np.min([diff_ice3.min().values,diff_lima.min().values])
    fig, axs = plt.subplots(2,3,figsize=(15,12),sharey=True, layout='constrained')
    create_plot(cfad_obs,axs[0,0],logNorm_cbar,'Radar observations',vmax)    # OBS
    create_plot(cfad_ice3,axs[0,1],logNorm_cbar,'AROME + ICE3',vmax)         # ICE3
    plot_diff(diff_ice3,axs[1,1],'ICE3 - OBS',vmin_diff,vmax_diff,logNorm_cbar)
    cs = create_plot(cfad_lima,axs[0,2],logNorm_cbar,'AROME + LIMA',vmax)    # LIMAsg
    diff = plot_diff(diff_lima,axs[1,2],'LIMASG - OBS',vmin_diff,vmax_diff,logNorm_cbar)
    if logNorm_cbar : fig.colorbar(cs,ax=axs[0,:], label='Count (LogNorm)',pad=0.01)
    else : fig.colorbar(cs,ax=axs[0,:], label='Count',pad=0.01)
    fig.colorbar(diff,ax=axs[1,:], label='Difference',pad=0.01)
    axs[1,0].set_axis_off()
else :
    fig, axs = plt.subplots(1,3,figsize=(15,6),sharey=True, layout='constrained')
    create_plot(cfad_obs,axs[0],logNorm_cbar,'Radar observations',vmax)    # OBS
    create_plot(cfad_ice3,axs[1],logNorm_cbar,'AROME + ICE3',vmax)         # ICE3
    cs = create_plot(cfad_lima,axs[2],logNorm_cbar,'AROME + LIMA',vmax)    # LIMAsg
    if logNorm_cbar : fig.colorbar(cs,ax=axs, label='Count (LogNorm)',pad=0.01)
    else : fig.colorbar(cs,ax=axs, label='Count',pad=0.01)
    
fig.supxlabel(xlabel)
for ax in axs.ravel() :
    ax.set_xlim(xlim)

figname = f'CFAD_{variable}_all_studyCases'
if logNorm_cbar :
    figname += '_logNorm_cbar'
if compare :
    figname += '_with_diff_plot'
      
fig.savefig(saveFigPath+figname+'.png', bbox_inches='tight',dpi=100)
