# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:23:04 2020

@author: borderiesm and augrosc 
Trace distributions D-Zhh et M-Zhh (ou Zdr, Kdp, Rhv)
pour bandes de frequence specifiee dans band_list
"""
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import math
import os


micro="ICE3" # ICE3_1M or ICJW_1M
moments={"ICE3":{"rr":"1M","ss":"1M","gg":"1M","wg":"1M","cl":"1M","ii":"1M"},
         "ICJW":{"rr":"1M","ss":"1M","gg":"1M","wg":"1M","cl":"1M","ii":"1M"},
         "LIMA":{"rr":"2M","ss":"1M","gg":"1M","wg":"1M","cl":"2M","ii":"2M"},
         "LIMC":{"rr":"2M","ss":"1M","gg":"1M","wg":"1M","cl":"2M","ii":"1M"}
         }

Path_tables = "../tables_generator/tables"
    
dir_fig="IMG/"

band,color_list,typeh_list ='C',['b','k','c','g','gray','gray'],['rr','cl','ss','ii','gg']#,'wg']
#band,color_list,typeh_list ='C',['b','c','gray'],['rr','hh','gg']

listplot=['M'] #,'M']
TmatOption={'C':'default','K':'default','W':'default'}
pltunit={'D':'mm','M':r'kg m$^{-3}$'}
pltX={'D':'Deq','M':'M'}
listvar=['Zh','Zdr','Kdp'] #,'Rhohv','Ah','Av']
unit={'Zh':'dBZ','Zdr':'dB','Kdp':u'\u00B0'+r' km$^{-1}$','Rhohv':'/','Ah':'dB','Av':'dB'}
vn={}


tick_size=16
pol_label=20
pol_legend=18
pol_title=26
pol_suptitle=22
lw=3


typeName = {'ii':'Pristine ice','ss':'Dry Snow','gg':'Dry Graupel','cl':'Cloud Water','rr':'Rain','wg':'Wet Graupel','hh':'Dry Hail','wh':'Wet Hail'}
typeName = {'ii':'ice','ss':'snow','gg':'graupel','cl':'cloud water','rr':'rain','wg':'wetgr','hh':'hail','wh':'wet hail'}
ymin_dict,ymax_dict={},{}
for var in ['Zh','Zdr','Rhohv','Kdp']:
    ymin_dict[var],ymax_dict[var]={},{}

ymin_dict["Zh"] = -30
ymax_dict["Zh"] = 70

ymin_dict["Zdr"] = -1
ymax_dict["Zdr"] = 6

ymin_dict["Kdp"] = -1
ymax_dict["Kdp"] = 6


for typeh in typeh_list:
    ymin_dict["Rhohv"][typeh]=0.9
    ymax_dict["Rhohv"][typeh]=1.0


dmax_dict={'ii':10,'ss':20,'gg':50,'cl':2,'rr':10,'wg':50,'hh':100,'wh':100}

Fwsel=0 # selected Fw for all species except graupel (= 0 because only graupe can be wet in ICE3/LIMA)
Fw_list,Fw_ls=[0.0,0.1,0.6,1.0],['-.',':','--','-']
ELEVsel=0 #0 for ground radars, 90 for RASTA or BASTA
Nii=800 #selected number concentration for primary ice
T_dict = {'ii':-30,'ss':-10,'gg':0,'cl':5,'rr':10,'wg':0,'hh':1,'wh':0}
Fw_list,Fw_ls=[0.1,0.6,0.9],['-.',':','--']
#T_list={'ii':[-30,-20,-10],'ss':[-20,-10,0],'gg':[-20,-10,0],'tt':[-10,0,10],'rr':[0,10,25],'wg':[-10,0],'hh':[-15,0,15],'wh':[-10,0,10]}



delim={'D':r"\s+",'M':';'}

legend_list = typeName

for plot in listplot:  
    # fig,axes = plt.subplots(nrows= 2,ncols=2,figsize=(15,15))
    fig,axes = plt.subplots(nrows= 1,ncols=3,figsize=(20,8))
    #fig,axes = plt.subplots(nrows= 1,ncols=1,figsize=(10,10))
    if (len(listvar))>1:
        ax = axes.flat
    else:
        ax=[axes,]
    
    for ivar,var in enumerate(listvar): 
        print("Plotting "+var+"-"+plot)
        for itypeh,typeh in enumerate(typeh_list):
            print("type : "+typeh)
            begtable={'D':'TmatCoefDiff_','M':'TmatCoefInt_'+micro+'_'+moments[micro][typeh]+'_'}
            
            print("band : "+band)
            table=Path_tables+"/"+TmatOption[band]+"/"+begtable[plot]+band+typeh

            # reading table options (from config file)
            df_param=pd.read_csv(table, sep=delim[plot],nrows=1,engine='python')
            SIGBETA=int(df_param['SIGBETA'][0])
            ARfunc=df_param['ARfunc'][0]
            ARcnst=df_param['ARcnst'][0]                
                
            # reading table variables
            df=pd.read_csv(table, sep=delim[plot],skiprows =2,engine='python')
            Tcol=df['Tc'].to_numpy()
            ELEV=df['ELEV'].to_numpy()
            if (plot=='M'):
                P3col=df['P3'].to_numpy()                           
            else:
                P3col=df['Fw'].to_numpy()

            vn['Zh'],vn['Zdr']=df['zhh'].to_numpy(),df['zdr'].to_numpy()
            vn['Rhohv'],vn['Kdp']=df['rhohv'].to_numpy(),df['kdp'].to_numpy()
            vn['Ah'],vn['Av']=df['Ah'].to_numpy(),df['Av'].to_numpy()
            x=df[pltX[plot]].to_numpy()

            
            if (typeh=='wg' or typeh=='wh'):
                for iFw,Fw in enumerate(Fw_list):
                    ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(P3col == Fw)) 
                    label=legend_list[typeh]+" Fw="+str(Fw)+" T="+str(T_dict[typeh])+u'\u00B0C'
                    ax[ivar].plot(x[ind],vn[var][ind],label = label,color=color_list[itypeh],ls=Fw_ls[iFw],linewidth=lw)#,s=0.5)
                if (ivar==2):
                    ax[ivar].legend(loc = 'best',fontsize=pol_legend)   
            else:                
                if (typeh=='ii' and plot=='M') :
                    id=np.argmin(np.abs(P3col-Nii)) # position of closest concentration to Nii
                    P3sel=P3col[id]
                    label=legend_list[typeh]+" Nii="+str(P3sel)+" T="+str(T_dict[typeh])+u'\u00B0C'
                    ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(P3col == P3sel))           
                    ax[ivar].plot(x[ind],vn[var][ind],label = label,color=color_list[itypeh],linewidth=lw)#,s=0.5)  
                else:
                    label=legend_list[typeh]+" T="+str(T_dict[typeh])+u'\u00B0C'
                    ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(P3col == Fwsel))           
                    ax[ivar].plot(x[ind],vn[var][ind],label = label,color=color_list[itypeh],linewidth=lw)#,s=0.5)               
                
            #end choice wg
        # end loop over typeh          

        ax[ivar].set_ylabel(var+"("+unit[var]+")",fontsize = pol_label)
        ax[ivar].set_xlabel(plot+"("+pltunit[plot]+")",fontsize = pol_label)
        ax[ivar].tick_params(axis='both', which='major', labelsize=tick_size)
        #ax[itypeh].set_xscale('log')
        ymin=ymin_dict[var]
        ymax=ymax_dict[var]
        ax[ivar].set_ylim(ymin,ymax)
        if plot=="D":
            ax[ivar].set_xlim(0,dmax_dict[typeh])
        if plot=="M":
            ax[ivar].set_xlim(1e-05,1e-02)
            ax[ivar].set_xscale('log')
        if (typeh != 'wg' and typeh!='wh'):
            ax[ivar].legend(loc = 'best',fontsize=pol_legend)
        ax[ivar].grid()
  
    # end loop over var        
    figtitle="Zh, Zdr, Kdp as a function of "+plot+" ("+pltunit[plot]+") "+"\n"+band+" band, "+"1-moment PSD (ICE3)"
    fig.suptitle(figtitle,fontsize=pol_suptitle)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    varpol=''.join(listvar)
    types=''.join(typeh_list)        
   
    micron=micro
    if (plot=='D'):
        micron=''
    nomfig=dir_fig +'DistTmat'+band+micro+"T"+str(T_dict['rr'])+varpol+types+'_'+plot+'varpol_alltypes'+TmatOption[band]+'.png'      
 
    fig.savefig(nomfig,dpi=200, bbox_inches='tight')
    plt.clf()
    plt.close('all')
    print("Figure saved in: "+nomfig)
    
#end loop over plot (D or M - var)

  

    
