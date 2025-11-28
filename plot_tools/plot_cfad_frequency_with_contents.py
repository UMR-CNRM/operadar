#!/usr/bin/env python
# coding: utf-8

import sys
import dask  
import argparse
import numpy as np
import xarray as xr
import pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator



# Plot settings
bins = {'Zh' : np.arange(-15, 70, 5),
        'Zdr': np.arange(-1, 4, 0.2),
        'Kdp': np.arange(-1, 4, 0.2),
       }
xticks_label = {'Zh' : '$Z_{H} (dBZ)$',
                'Zdr':'$Z_{DR} (dB)$',
                'Kdp':'$K_{DP} (°/km)$',
               }
xticks_step = {'Zh' : 10,
               'Zdr': 1,
               'Kdp': 1,
              }
minor_ticks = {'Zh' : MultipleLocator(2.5),
               'Zdr': MultipleLocator(0.2),
               'Kdp': MultipleLocator(0.2),
              }
percent_max = {'Zh' : 25,
               'Zdr': 35,
               'Kdp': 50,
              }
level = np.arange(0,90,1)


# ============== FUNCTIONS ============== #
def cfad_frequency_with_quantiles(data,xbins,zbins):
    # Frequency distribution per zbin
    task_list=[]
    for z in zbins:
        hist_z = dask.delayed(np.histogram)(data.sel(level=z), bins=xbins)
        task_list.append(hist_z[0])
    hist = np.array(dask.compute(*task_list))
    # Quantiles
    qDict = dict(altitude=zbins)
    qDict['Q25'] = np.nanpercentile(data,q=25,axis=(1,2))[:len(zbins)]
    qDict['Q50'] = np.nanpercentile(data,q=50,axis=(1,2))[:len(zbins)]
    qDict['Q75'] = np.nanpercentile(data,q=75,axis=(1,2))[:len(zbins)]
    qDict['Q95'] = np.nanpercentile(data,q=95,axis=(1,2))[:len(zbins)]
    
    return hist, qDict


def plot_cfad(data,zbins,h,ax,color,add_extreme=False,add_Q25_Q75_interval=True,alpha=0.2):
    # Q50
    q50 = np.nanpercentile(data,q=50,axis=(1,2))[:len(zbins)]
    plot, = ax.plot(q50,zbins,label=f'{h} Q50',c=color)
    # Q25 - Q75
    if add_Q25_Q75_interval :
        q25 = np.nanpercentile(data,q=25,axis=(1,2))[:len(zbins)]
        q75 = np.nanpercentile(data,q=75,axis=(1,2))[:len(zbins)]
        plot = ax.fill_betweenx(zbins,q25,q75,alpha=alpha,label='[Q25-Q75]',color=color)
    # Q95
    if add_extreme :
        q95 = np.nanpercentile(data,q=95,axis=(1,2))[:len(zbins)]
        plot = ax.plot(q95,zbins,label='Q95',c=color,ls=':',lw=1.5)
    return plot



def main(var:str, micro:str, filepath:str, improved:str, percentage:bool):
    print(f'----- var:{var} improved:{improved} percentage:{percentage} --------------------------')
    
    if not improved:
        extension = ''
        subtitle = f'{micro} (REF)'
    else : 
        extension = '_improved'
        subtitle = f'{micro} (NEW)'

    # Data selection
    print('Reading data')
    ds = xr.open_dataset(filepath+'.nc')
    time = ds.time.values.astype('datetime64[us]').item().strftime("%Y%m%d%H%M")
    mask = (ds.Zh >= 5)
    data = ds[var]#.where(mask) # A CHANGER
    
    # Obtaining frequency CFAD and quantiles
    print('Filling CFAD bins')
    hist, qDict = cfad_frequency_with_quantiles(data, xbins=bins[var], zbins=level)
    
    # Plot figure
    print('Building figure...')
    fig, axs = plt.subplots(1,2,figsize=(12,9),sharey=True,layout='constrained')
    
    # Frequencies
    if percentage :
        print('   Plot in percentage (% per altitude level)')
        xdim, ydim = len(bins[var])-1, len(level)
        sum_per_alt_level = np.repeat(hist.sum(axis=1),xdim).reshape(ydim,xdim)
        proportion_per_alt_level = hist/sum_per_alt_level*100
        cs = axs[0].pcolor(bins[var][:-1], level, proportion_per_alt_level,
                            cmap='jet', vmin=0,vmax=percent_max[var],
                            ) ; cs.set_edgecolor('face')
        label = 'Distribution by altitude level (%)'
    else :
        print('   Plot in absolute frequency')
        cs = axs[0].pcolor(bins[var][:-1], level, hist,
                            cmap='jet', norm = LogNorm(vmin=1,vmax=1e6),
                            ) ; cs.set_edgecolor('face')
        label = 'Absolute frequency'
    # Percentiles
    print('   Overlaying the 25, 50 and 75th percentiles')
    axs[0].plot(qDict['Q25'], level, 'k', label='Q25', ls='-.')
    axs[0].plot(qDict['Q50'], level, 'k', label='Q50')
    axs[0].plot(qDict['Q75'], level, 'k', label='Q75', ls='--')
    
    # X-axis
    axs[0].set_title(subtitle,fontsize=14)
    axs[0].xaxis.set_minor_locator(minor_ticks[var])
    axs[0].set_xticks(np.arange(bins[var][0],bins[var][-1]+xticks_step[var],xticks_step[var]),
                        labels = np.arange(bins[var][0],bins[var][-1]+xticks_step[var],xticks_step[var]).astype(np.int32),fontsize=11)
    axs[0].set_xlim(bins[var][0],bins[var][-1])
    axs[0].set_xlabel(f'{xticks_label[var]}',fontsize=12)
    axs[0].tick_params(axis='both', which='major', labelsize=11)

    cbar = fig.colorbar(cs, ax=axs[0], pad=0.02, extend='max') ; cbar.set_label(label=label,size=11) ; cbar.ax.tick_params(labelsize=11) 
    axs[0].legend(fontsize=11)

    # Plot contents
    print('Adding the corresponding model hydrometeor contents...')
    listHydro = [('cc','tab:cyan'),('rr','tab:blue'),('gg','maroon'),('ss','deeppink'),('ii','red')]
    for (h,c) in listHydro :
        da = ds['Contents'].sel(hydrometeor=h).where(ds['Contents'].sel(hydrometeor=h)>0)#.where(mask) # A CHANGER
        plot_cfad(data=da, zbins=level, h=h, ax=axs[1], color=c )
        del da
    # Grid
    axs[1].grid(which='major',ls='--',lw=0.6)
    axs[1].grid(True,ls='--',lw=0.3,which='minor',axis='both')
    # X-axis
    axs[1].set_title(f'Contents where content[hydrometeor]>0',fontsize=14)
    axs[1].xaxis.set_minor_locator(MultipleLocator(0.05))
    axs[1].set_xticks(np.arange(0,1.2,0.2))
    axs[1].set_xlim(0,1)
    axs[1].set_xlabel(f'Contents (g/m$^3$)',fontsize=12)
    axs[1].tick_params(axis='both', which='major', labelsize=11)
    #  Legend
    axs[1].legend(fontsize=11)

    # Y-axis
    axs[0].yaxis.set_minor_locator(MultipleLocator(10))
    axs[0].set_ylim(89,9) ; axs[1].set_ylim(89,9)
    axs[0].set_ylabel('AROME level',fontsize=12)

    print('Saving...')
    if percentage :
        plt.savefig(f'{fig_output_path}cfad_{time}_level_percentage_{var}_{micro}{extension}.png',bbox_inches='tight',dpi=100)
    else :
        plt.savefig(f'{fig_output_path}cfad_{time}_level_absFreq_{var}_{micro}{extension}.png',bbox_inches='tight',dpi=100)
    

# ================= MAIN ================= #    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="------ CFADs with frequency plot tool ------")
    parser.add_argument("filepath", type=str,
                        help='Complete path file (without .nc extension)')
    parser.add_argument("var", type=str, default='Zh', 
                        help='Chose one polarimetric variable among : Zh, Zdr, Kdp')
    parser.add_argument("micro", type=str, default='ICE3',
                        help='Chose one microphysics among : ICE3 ICJW LIMASG LIMC')
    parser.add_argument("--improved", action="store_true",default=False,
                        help='Wether to plot with the improved OPERADAR configuration or not')
    parser.add_argument("--percentage", action="store_true",default=False,
                        help='Wether to plot or not the data and colorbar in %% per altitude level instead of absolute frequency.')
    
    fig_output_path='/home/davidcl/Programmation/output/chapter_4_manuscript/' # A CHANGER
    parser_args = parser.parse_args()
    
    main(filepath=parser_args.filepath,
         var=parser_args.var,
         micro=parser_args.micro,
         improved=parser_args.improved,
         percentage=parser_args.percentage,
         )
