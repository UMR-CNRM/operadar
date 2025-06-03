#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt


# # CFAD
# Faire tourner le script `stats_for_CFADs.py`
# 
# **La conversion en pickle est bien mieux que csv car elle permet le stockage des valeurs pol sous forme de array dans le dataframe**

# ======== Functions ========================
def plot_CFAD(ax,df,var,dataType,color='tab:blue',add_extreme=False,add_Q25_Q75_interval=True):
    subdf = df[(df.varname==var)&(df.dataType==dataType)]
    alti = subdf.altitude
    #plot, = ax.plot(subdf.Q50,alti,label=f'{var} Q50',c=color)
    plot, = ax.plot(subdf.Q50,alti,label=f'{dataType}',c=color,lw=2)
    if add_Q25_Q75_interval :
        plot = ax.fill_betweenx(alti,subdf.Q25,subdf.Q75,alpha=0.2,label='[Q25-Q75]',color=color)
        #ax.plot(subdf.Q25,alti,c=color,ls=':',lw=1.5) # commenter cette ligne pour ne pas l'afficher
        #plot = ax.plot(subdf.Q75,alti,c=color,ls=':',lw=1.5)
    if add_extreme :
        ax.plot(subdf.Q5,alti,c=color,ls='--',lw=1) # commenter cette ligne pour ne pas l'afficher
        plot = ax.plot(subdf.Q95,alti,c=color,ls=':',lw=1.5)
    return plot


def legend_hydrometeor(figure,h,l):
    legend_name = {'rr':'rain','ii':'ice','ss':'snow','gg':'graupel','wg':'wet graupel','cc':'cloud water'}
    labels = []
    for lgd_elmt in l :
        if lgd_elmt[0:2] in legend_name.keys():
            labels += [legend_name[lgd_elmt[0:2]]]
        else :
            labels += [lgd_elmt]
    figure.legend(handles=h,
               labels=labels,
               bbox_to_anchor=(0.5, 0.5),loc='center left',ncol=1,fontsize=12)
# ============================================    


#========== Main programm =======

# ==== Data path and simulations list
dataPath = "/home/cnrm_other/ge/mrmp/augros/WKD/CORSE/"

simuList = ["obs","MesoNH_ICE3","Arome_ICE3","Arome_LIMA"] #,"MesoNH_LIMAAG"] #,"Arome_ICE3","Arome_LIMA","MesoNH_ICE3","MesoNH_LIMA","MesoNH_LIMAAG"]
colorDict = {"obs":"tab:blue",
             "Arome_oper":"tab:red",
             "Arome_ICE3":"tab:orange",
             "Arome_LIMA":"tab:green",
             "MesoNH_ICE3":'tab:red',
             "MesoNH_LIMA":'tab:blue',
             "MesoNH_LIMAAG":'tab:green',
    }

# ==== Radar variables
#df = pd.read_pickle(f'{final_data}stats_cell_core.pkl')
#df = pd.read_pickle(f'{final_data}stats_cell_envelop.pkl')
#df = pd.read_pickle(f'{final_data}stats_ZDR_columns.pkl')



var = 'zh' ; title = '$Z_{H}$' ; unit = 'dBZ' ; minv=10 ; maxv=60;
var = 'zdr' ; title = '$Z_{dr}$' ; unit = 'dB' ; minv=0 ; maxv=3;
var = 'kdp' ; title = '$K_{dp}$' ; unit = 'deg/km' ; minv=0 ; maxv=3;
#var = 'rhohv' ; title = '$rho_{hv}$' ; unit = '/' ; minv=0.9 ; maxv=1;
# subdf = df[(df.varname==var)&(df.dataType==dataType)]

threshold="40"

fig, ax = plt.subplots(1,1,figsize=(8,5))

for simu in simuList:
    df = pd.read_pickle(f'{dataPath}df_for_stats_zhsup{threshold}_{simu}.pkl')
    plot_CFAD(ax,df,var,simu,colorDict[simu])
    del df


ax.legend(loc='upper right',ncol=1,prop={'size':13})

ax.grid(ls='--',lw=0.5)
ax.set_xlim(minv,maxv)
#ax.set_yticks(np.arange(0,16e3,1e3),labels=np.arange(0,16,1),fontsize=13)
#ax.set_xticks(np.arange(10,60,5),labels=np.arange(10,60,5),fontsize=13)
ax.set_ylim(0, 12e3)
ax.set_ylabel('Height (m)',fontsize=14)
ax.set_xlabel(f'{title} ({unit})',fontsize=14)
#ax.set_title(f'Distribution of {title} values over 9 study cases inside storms cores',fontsize=15,y=1.02)
plt.savefig(dataPath+f'cfad_zhsup{threshold}_{var}',bbox_inches='tight',dpi=100)


# ### Contenus en hydrométéores (uniquement pour les modèles)

# # In[7]:


altiList = [int(x) for x in np.arange(0,15e3,500)]
fig, axs = plt.subplots(1,3,figsize=(15,6),sharey=True)

for index,simu in enumerate(['MesoNH_ICE3','Arome_ICE3','Arome_LIMA']):
    df = pd.read_pickle(f'{dataPath}df_for_stats_zhsup{threshold}_{simu}.pkl')
    rr = plot_CFAD(axs[index],df,'rr',simu,'tab:blue',add_extreme=False,add_Q25_Q75_interval=True)
    ii = plot_CFAD(axs[index],df,'ii',simu,'tab:olive',add_extreme=False,add_Q25_Q75_interval=True)
    ss = plot_CFAD(axs[index],df,'ss',simu,'tab:pink',add_extreme=False,add_Q25_Q75_interval=True)
    gg = plot_CFAD(axs[index],df,'gg',simu,'tab:gray',add_extreme=False,add_Q25_Q75_interval=True)
    wg = plot_CFAD(axs[index],df,'wg',simu,'tab:purple',add_extreme=False,add_Q25_Q75_interval=True)
    cc = plot_CFAD(axs[index],df,'cc',simu,'tab:cyan',add_extreme=False,add_Q25_Q75_interval=True)
    axs[index].legend(handles=[rr,ii,ss,gg,wg,cc], labels=['rain','ice','snow','graupel','wet graupel','cloud water'])
    axs[index].grid(ls='--',lw=0.5)
    axs[index].set_xlim(0,3)
    #axs[index].set_xscale('log')
    #axs[index].set_yticks(np.arange(0,16e3,1e3),labels=np.arange(0,16,1))
    axs[index].set_ylim(0, 15e3)
    axs[0].set_ylabel('Height (km)',fontsize=11)
    axs[index].set_xlabel('Hydrometeor content ($g/m^{3}$)',fontsize=11)
    axs[index].set_title(simu)
    h, l = axs[index].get_legend_handles_labels()
    del df
#axs[index].legend(loc='upper right',ncol=1,prop={'size':13})
#legend_hydrometeor(fig,h,l)
fig.suptitle(f'Distribution of hydrometeors content inside cell cores (zh>{threshold} dBZ)',fontsize=13)
plt.savefig(dataPath+f'cfad_hydrometeor_contents_zhsup{threshold}',dpi=120, bbox_inches='tight')


# # In[8]:


# altiList = [int(x) for x in np.arange(0,15e3,500)]
# fig, axs = plt.subplots(1,2,figsize=(10,6),sharey=True,layout='constrained')

# for index,model in enumerate(['ICE3','LIMASG']):
#     rr = plot_CFAD(axs[index],df,'rr',model,'tab:green',add_extreme=False,add_Q25_Q75_interval=True)
#     ii = plot_CFAD(axs[index],df,'ii',model,'tab:blue',add_extreme=False,add_Q25_Q75_interval=False)
#     ss = plot_CFAD(axs[index],df,'ss',model,'tab:red',add_extreme=False,add_Q25_Q75_interval=False)
#     gg = plot_CFAD(axs[index],df,'gg',model,'tab:orange',add_extreme=False,add_Q25_Q75_interval=False)
#     wg = plot_CFAD(axs[index],df,'wg',model,'tab:purple',add_extreme=True,add_Q25_Q75_interval=True)
#     cc = plot_CFAD(axs[index],df,'cc',model,'tab:cyan',add_extreme=True,add_Q25_Q75_interval=False)
#     #axs[index].legend(handles=[rr,ii,ss,gg,wg,cc], labels=['rain','ice','snow','graupel','wet graupel','cloud water'])
#     axs[index].grid(ls='--',lw=0.5)
#     axs[index].set_xlim(0,2)
#     #axs[index].set_xscale('log')
#     axs[index].set_yticks(np.arange(0,16e3,1e3),labels=np.arange(0,16,1))
#     axs[index].set_ylim(0, 15e3)
#     axs[0].set_ylabel('Height (km)',fontsize=11)
#     axs[index].set_xlabel('Hydrometeor content ($g/m^{3}$)',fontsize=11)
#     axs[index].set_title(model)
#     h, l = axs[index].get_legend_handles_labels()
# legend_hydrometeor(fig,h,l)
# fig.suptitle(f'Distribution of hydrometeors content over 10 study cases inside storm cells',fontsize=13)
# #plt.savefig(f'{out_repo_path}statistiques/CFAD_hydrometeor_contents_all_studyCases_cell_envelop',dpi=150, bbox_inches='tight')


# # In[9]:


# altiList = [int(x) for x in np.arange(0,15e3,500)]
# fig, axs = plt.subplots(1,2,figsize=(10,6),sharey=True,layout='constrained')

# for index,model in enumerate(['ICE3','LIMASG']):
#     rr = plot_CFAD(axs[index],df,'rr',model,'tab:green',add_extreme=True,add_Q25_Q75_interval=True)
#     ii = plot_CFAD(axs[index],df,'ii',model,'tab:blue',add_extreme=False,add_Q25_Q75_interval=False)
#     ss = plot_CFAD(axs[index],df,'ss',model,'tab:red',add_extreme=False,add_Q25_Q75_interval=False)
#     gg = plot_CFAD(axs[index],df,'gg',model,'tab:orange',add_extreme=False,add_Q25_Q75_interval=False)
#     wg = plot_CFAD(axs[index],df,'wg',model,'tab:purple',add_extreme=True,add_Q25_Q75_interval=True)
#     cc = plot_CFAD(axs[index],df,'cc',model,'tab:cyan',add_extreme=True,add_Q25_Q75_interval=False)
#     #axs[index].legend(handles=[rr,ii,ss,gg,wg,cc], labels=['rain','ice','snow','graupel','wet graupel','cloud water'])
#     axs[index].grid(ls='--',lw=0.5)
#     axs[index].set_xlim(0,3.5)
#     #axs[index].set_xscale('log')
#     axs[index].set_yticks(np.arange(0,16e3,1e3),labels=np.arange(0,16,1))
#     axs[index].set_ylim(0, 15e3)
#     axs[0].set_ylabel('Height (km)',fontsize=11)
#     axs[index].set_xlabel('Hydrometeor content ($g/m^{3}$)',fontsize=11)
#     axs[index].set_title(model)
#     h, l = axs[index].get_legend_handles_labels()
# legend_hydrometeor(fig,h,l)
# fig.suptitle(f'Distribution of hydrometeors content over 10 study cases inside storms cores',fontsize=13)
# #plt.savefig(f'{out_repo_path}statistiques/CFAD_hydrometeor_contents_all_studyCases_cell_core.png',dpi=150, bbox_inches='tight')


# # In[ ]:





# # In[ ]:




