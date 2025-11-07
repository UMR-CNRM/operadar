# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:23:04 2020

@author: borderiesm and augrosc 
Trace distributions D-Zhh et M-Zhh (ou Zdr, Kdp, Rhohv ou Ah)
pour bandes de frequence specifiee dans band_list
"""
import numpy as np
import matplotlib;
matplotlib.use('Agg')
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import pandas as pd
import math
import os


micro="LIMC" # ICE3_1M or ICJW_1M
moments={"ICE3":{"rr":"1M","ss":"1M","gg":"1M","wg":"1M","cl":"1M","cs":"1M","ii":"1M"},
         "ICJW":{"rr":"1M","ss":"1M","gg":"1M","wg":"1M","cl":"1M","cs":"1M","ii":"1M"},
         "LIMA":{"rr":"2M","ss":"1M","gg":"1M","wg":"1M","cl":"2M","cs":"2M","ii":"2M"},
         "LIMC":{"rr":"2M","ss":"1M","gg":"1M","wg":"1M","cl":"2M","cs":"2M","ii":"1M"}
         }
plotR=True #True 

Path_tables = "../tables_generator/tables"
    
dir_fig="IMG/"

#band_list,typeh_list =['S','C','X','K','W'],['rr','ss','gg','hh','cc','ii','wg','wh']
band_list,typeh_list =['C'],['rr','ss','gg','cl','ii','wg'] #'rr','ss','gg','cc','ii','wg'] #,'cc','ii','wg'] #,'ss','gg','ii'] #,'gg'] #,'ii'] #,'wg'] #,'gg','wg','cc','ii']

band_name={'S':'S','C':'C','X':'X','K':'Ka','W':'W'}
color={'S':'k','C':'blue','X':'green','K':'gray','W':'green' }#'lightgray'}
style={'S':'-','C':'-','X':'-','K':'-','W':'-'}

listplot=['M'] #,'M']
TmatOption={'C':'default','K':'default','W':'default'}
pltunit={'D':'mm','M':r'kg m$^{-3}$'}
pltX={'D':'Deq','M':'M'}
listvar=['Zh']
unit={'Zh':'dBZ','Zdr':'dB','Kdp':u'\u00B0'+r' km$^{-1}$','Rhohv':'/','Ah':'dB'+r' km$^{-1}$','Av':'dB'+r' km$^{-1}$'}
vn={}
vn_R={}

pol_label=20
pol_legend=20
pol_title=22
pol_suptitle=30
lw=3


typeName = {'ii':'Pristine ice','ss':'Dry Snow','gg':'Dry Graupel','cl':'Cloud Water','cs':'Cloud Water','rr':'Rain','wg':'Wet Graupel','hh':'Dry Hail','wh':'Wet Hail'}
ymin_dict,ymax_dict={},{}
for var in ['Zh','Zdr','Rhohv','Kdp','Ah','Av']:
    ymin_dict[var],ymax_dict[var]={},{}

#ymin_dict["Zh"] = {'ii':-50,'ss':-60,'gg':-40,'cl':-60,'cc':-60,'rr':-20,'wg':-10,'hh':20,'wh':20}
#ymax_dict["Zh"] = {'ii':30,'ss':30,'gg':80,'cl':20,'cc':20,'rr':70,'wg':80,'hh':100,'wh':100}
ymin_dict["Zh"] = {'ii':-20,'ss':-20,'gg':-20,'cl':-20,'cs':-20,'rr':-20,'wg':-20,'hh':-20,'wh':-20}
ymax_dict["Zh"] = {'ii':70,'ss':70,'gg':70,'cl':70,'cs':70,'rr':70,'wg':70,'hh':70,'wh':70}

ymin_dict["Zdr"] = {'ii':0,'ss':-2,'gg':-2,'cl':0,'cs':0,'rr':-4,'wg':-4,'hh':-2,'wh':-4}
ymax_dict["Zdr"] = {'ii':6,'ss':2,'gg':2,'cl':1,'cs':1,'rr':10,'wg':10,'hh':2,'wh':10}


ymin_dict["Kdp"] = {'ii':0,'ss':0,'gg':-0.2,'cl':0,'cs':0,'rr':-2,'wg':-4,'hh':-30,'wh':-30}
ymax_dict["Kdp"] = {'ii':0.2,'ss':0.2,'gg':0.2,'cl':1,'cs':1,'rr':5,'wg':4,'hh':20,'wh':20}


ymin_dict["Ah"] = {'ii':0,'ss':0,'gg':0,'cl':0,'rr':0,'wg':0,'hh':0,'wh':0}
ymax_dict["Ah"] = {'ii':10,'ss':10,'gg':10,'cl':10,'rr':10,'wg':10,'hh':10,'wh':10}

ymin_dict["Av"] = {'ii':0,'ss':0,'gg':0,'cl':0,'rr':0,'wg':0,'hh':0,'wh':0}
ymax_dict["Av"] = {'ii':10,'ss':10,'gg':10,'cl':10,'rr':10,'wg':10,'hh':10,'wh':10}

for typeh in typeh_list:
    ymin_dict["Rhohv"][typeh]=0.9
    ymax_dict["Rhohv"][typeh]=1.0


dmax_dict={'ii':10,'ss':20,'gg':50,'cl':2,'cs':2,'rr':10,'wg':50,'hh':100,'wh':100}

Fwsel=0 # selected Fw for all species except graupel (= 0 because only graupe can be wet in ICE3/LIMA)
Fw_list,Fw_ls=[0.0,0.1,0.6,1.0],['-.',':','--','-']
ELEVsel=0 #0 for ground radars, 90 for RASTA or BASTA
Nii=800 #selected number concentration for primary ice
expNmin,expNmax=1,8
nexpN=expNmax-expNmin+1
expN_list,color_list= np.arange(expNmin,expNmax+1), plt.cm.jet(np.linspace(0,1,nexpN))
T_dict = {'ii':-30,'ss':-10,'gg':0,'cl':5,'cs':5,'rr':10,'wg':10,'hh':1,'wh':10}

#T_list={'ii':[-30,-20,-10],'ss':[-20,-10,0],'gg':[-20,-10,0],'cc':[-10,0,10],'rr':[0,10,25],'wg':[-10,0],'hh':[-15,0,15],'wh':[-10,0,10]}


delim={'D':r"\s+",'M':';'}

for var in listvar:
    for plot in listplot:    
        print("Plotting "+var+"-"+plot)
        fig,axes = plt.subplots(nrows= 2,ncols=3,figsize=(24,12))
        ax = axes.flat
        for itypeh,typeh in enumerate(typeh_list):
            print("type : "+typeh)
            begtable={'D':'TmatCoefDiff_','M':'TmatCoefInt_'+micro+'_'+moments[micro][typeh]+'_'}
                        
            for iband,band in enumerate(band_list): 
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
                vn_R['Zh'],vn_R['Zdr']=df['zhhR'].to_numpy(),df['zdrR'].to_numpy()
                vn_R['Rhohv'],vn_R['Kdp']=df['rhohvR'].to_numpy(),df['kdpR'].to_numpy()
                vn_R['Ah'],vn_R['Av']=df['AhR'].to_numpy(),df['AvR'].to_numpy()
                x=df[pltX[plot]].to_numpy()
                
                if (typeh=='wg' or typeh=='wh'):
                    for iFw,Fw in enumerate(Fw_list):
                        ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(P3col == Fw)) 
                        label=band_name[band]+" Fw="+str(Fw)
                        ax[itypeh].plot(x[ind],vn[var][ind],label = label,color=color[band],ls=Fw_ls[iFw],linewidth=lw)#,s=0.5)
                    if (iband==0):
                        ax[itypeh].legend(loc = 'best',fontsize=pol_legend)   
                else:
                    if ((typeh=='ii' and plot=='M') or (plot =='M' and moments[micro][typeh]=="2M")) :
                        expNcol=np.log10(P3col) 
                        for iN,expN in enumerate(expN_list):
                            id=np.argmin(np.abs(expNcol-expN)) # position of closest concentration to expN
                            expNsel=expNcol[id]
                            ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(expNcol == expNsel))  
                            label=" expN="+str(expNsel)
                            ax[itypeh].plot(x[ind],vn[var][ind],label = label,color=color_list[iN],linewidth=lw)#,s=0.5)
                            ax[itypeh].legend(loc = 'best',fontsize=pol_legend) 
                    else:
                        ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVsel)*(P3col == Fwsel))           
                        ax[itypeh].plot(x[ind],vn[var][ind],label = band_name[band],color=color[band],linewidth=lw)#,s=0.5)
                        if (plotR):
                            ax[itypeh].plot(x[ind],vn_R[var][ind],label = band_name[band]+', Rayleigh',color=color[band],linewidth=lw,ls='--')
                    
                #end choice wg
            # end loop over bands
            title1= typeName[typeh]+" T="+str(int(T_dict[typeh]))+u'\u00B0C'+'\n'
            title2= r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0'
            title3=" AR="+str(ARcnst)+" elev="+str(ELEVsel)
            if (typeh=="rr"):
                title3=" AR="+ARfunc
            ax[itypeh].set_title(title1+title2+title3,fontsize=pol_title)
            ax[itypeh].set_ylabel(var+"("+unit[var]+")",fontsize = pol_label)
            ax[itypeh].set_xlabel(plot+"("+pltunit[plot]+")",fontsize = pol_label)
            ax[itypeh].tick_params(axis='x', labelsize=pol_label)
            ax[itypeh].tick_params(axis='y', labelsize=pol_label)
            #ax[itypeh].set_xscale('log')
            ymin=ymin_dict[var][typeh]
            ymax=ymax_dict[var][typeh]
            ax[itypeh].set_ylim(ymin,ymax)
            if plot=="D":
                ax[itypeh].set_xlim(0,dmax_dict[typeh])
            if plot=="M":
                ax[itypeh].set_xlim(1e-05,1e-02)
                ax[itypeh].set_xscale('log')
            if (typeh != 'wg' and typeh!='wh'):
                ax[itypeh].legend(loc = 'best',fontsize=pol_legend)
            ax[itypeh].grid()
        # end loop over typeh    
        
        figtitle=var+"("+unit[var]+") as a function of "+plot+"("+pltunit[plot]+") - ELEV="+str(ELEVsel)+u'\u00B0'+" - "+micro
        fig.suptitle(figtitle,fontsize=pol_suptitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        bandtit=''
        for band in band_list:
            bandtit+=band
        
        micron=micro
        if (plot=='D'):
            micron=''

        nomfig=dir_fig +'DistTmat'+bandtit+micron+'_'+plot+var+TmatOption[band]
        if (plotR):
            nomfig+="_R"
        
        if not os.path.exists(dir_fig):
            os.makedirs(dir_fig)    
        fig.savefig(nomfig+".png",dpi=200, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        print("Figure saved in: "+nomfig)
    #end loop over plot (D or M - var)
#end loop over var

# # Histogram zh
# fig, ax = plt.subplots( nrows=1, ncols=1 )  
# plt.hist(vn["Zh"])
# fig.savefig(nomfig+"zh_hist.png")

    
