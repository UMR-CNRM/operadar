# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:23:04 2020

@author: borderiesm and augrosc 
Trace distributions D-Zhh et M-Zhh (ou Zdr, Kdp, Rhohv)
pour bandes de frequence specifiee dans band_list
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import math
import os


micro="ICE3"
plotR=True #True 

Path_tables = "../tables_generator/tables"
    
dir_fig="IMG/"

#band_list,typeh_list =['S','C','X','K','W'],['rr','ss','gg','hh','cc','ii','wg','wh']
band_list,typeh_list =['C','K','W'],['rr','ss','gg','cc','ii','wg'] #,'cc','ii','wg'] #,'ss','gg','ii'] #,'gg'] #,'ii'] #,'wg'] #,'gg','wg','cc','ii']

band_name={'S':'S','C':'C','X':'X','K':'Ka','W':'W'}
color={'S':'k','C':'blue','X':'green','K':'gray','W':'green' }#'lightgray'}
style={'S':'-','C':'-','X':'-','K':'-','W':'-'}

listplot=['D','M'] #,'M']
TmatOption={'C':"elev90",'K':'default','W':'default'}
pltunit={'D':'mm','M':r'kg m$^{-3}$'}
pltX={'D':'Deq','M':'M'}
listvar=['Zh','Ah']
unit={'Zh':'dBZ','Zdr':'dB','Kdp':u'\u00B0'+r' km$^{-1}$','Rhohv':'/','Ah':'dB'+r' km$^{-1}$','Av':'dB'+r' km$^{-1}$'}
vn={}
vn_R={}

pol_label=20
pol_legend=20
pol_title=22
pol_suptitle=30
lw=3


typeName = {'ii':'Pristine ice','ss':'Dry Snow','gg':'Dry Graupel','cc':'Cloud Water','rr':'Rain','wg':'Wet Graupel','hh':'Dry Hail','wh':'Wet Hail'}
ymin_dict,ymax_dict={},{}
for var in ['Zh','Zdr','Rhohv','Kdp','Ah','Av']:
    ymin_dict[var],ymax_dict[var]={},{}

ymin_dict["Zh"] = {'ii':-50,'ss':-60,'gg':-40,'cc':-60,'rr':-20,'wg':-10,'hh':20,'wh':20}
ymax_dict["Zh"] = {'ii':30,'ss':30,'gg':80,'cc':20,'rr':70,'wg':80,'hh':100,'wh':100}

ymin_dict["Zdr"] = {'ii':0,'ss':-2,'gg':-2,'cc':0,'rr':-4,'wg':-4,'hh':-2,'wh':-4}
ymax_dict["Zdr"] = {'ii':6,'ss':2,'gg':2,'cc':1,'rr':10,'wg':10,'hh':2,'wh':10}

ymin_dict["Kdp"] = {'ii':0,'ss':0,'gg':-0.2,'cc':0,'rr':-2,'wg':-4,'hh':-30,'wh':-30}
ymax_dict["Kdp"] = {'ii':0.2,'ss':0.2,'gg':0.2,'cc':1,'rr':5,'wg':4,'hh':20,'wh':20}

ymin_dict["Ah"] = {'ii':0,'ss':0,'gg':0,'cc':0,'rr':0,'wg':0,'hh':0,'wh':0}
ymax_dict["Ah"] = {'ii':10,'ss':10,'gg':10,'cc':10,'rr':10,'wg':10,'hh':10,'wh':10}

ymin_dict["Av"] = {'ii':0,'ss':0,'gg':0,'cc':0,'rr':0,'wg':0,'hh':0,'wh':0}
ymax_dict["Av"] = {'ii':10,'ss':10,'gg':10,'cc':10,'rr':10,'wg':10,'hh':10,'wh':10}

for typeh in typeh_list:
    ymin_dict["Rhohv"][typeh]=0.9
    ymax_dict["Rhohv"][typeh]=1.0


dmax_dict={'ii':10,'ss':20,'gg':50,'cc':2,'rr':10,'wg':50,'hh':100,'wh':100}

Fwchoix=0
ELEVchoix=90 #0 pour radars sol, 90 pour rasta
Nii=800 #selected number concentration for primary ice
T_dict = {'ii':-30,'ss':-10,'gg':0,'cc':5,'rr':10,'wg':10,'hh':1,'wh':10}
Fw_list,Fw_ls=[0.0,0.1,0.6,1.0],['-.',':','--','-']
#T_list={'ii':[-30,-20,-10],'ss':[-20,-10,0],'gg':[-20,-10,0],'cc':[-10,0,10],'rr':[0,10,25],'wg':[-10,0],'hh':[-15,0,15],'wh':[-10,0,10]}

begtable={'D':'TmatCoefDiff_','M':'TmatCoefInt_'+micro+'_'}
delim={'D':r"\s+",'M':';'}

for var in listvar:
    for plot in listplot:    
        print("Plotting "+var+"-"+plot)
        fig,axes = plt.subplots(nrows= 2,ncols=3,figsize=(24,12))
        ax = axes.flat
        for itypeh,typeh in enumerate(typeh_list):
            print("type : "+typeh)          
            
            
            for iband,band in enumerate(band_list): 
                print("band : "+band)
                #table=Path_tables+"/"+typeh+"/"+begtable[plot]+band+typeh
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
                    Fwcol=df['P3'].to_numpy()                            
                else:
                    Fwcol=df['Fw'].to_numpy()
                vn['Zh'],vn['Zdr']=df['zhh'].to_numpy(),df['zdr'].to_numpy()
                vn['Rhohv'],vn['Kdp']=df['rhohv'].to_numpy(),df['kdp'].to_numpy()
                vn['Ah'],vn['Av']=df['Ah'].to_numpy(),df['Av'].to_numpy()
                vn_R['Zh'],vn_R['Zdr']=df['zhhR'].to_numpy(),df['zdrR'].to_numpy()
                vn_R['Rhohv'],vn_R['Kdp']=df['rhohvR'].to_numpy(),df['kdpR'].to_numpy()
                vn_R['Ah'],vn_R['Av']=df['AhR'].to_numpy(),df['AvR'].to_numpy()
                x=df[pltX[plot]].to_numpy()
                
                if (typeh=='wg' or typeh=='wh'):
                    for iFw,Fw in enumerate(Fw_list):
                        ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVchoix)*(Fwcol == Fw)) 
                        label=band_name[band]+" Fw="+str(Fw)
                        ax[itypeh].plot(x[ind],vn[var][ind],label = label,color=color[band],ls=Fw_ls[iFw],linewidth=lw)#,s=0.5)
                    if (iband==0):
                        ax[itypeh].legend(loc = 'best',fontsize=pol_legend)   
                else:
                    if (typeh=='ii' and plot=='M') :
                        id=np.argmin(np.abs(Fwcol-math.log10(Nii))) # position of closest concentration to Nii
                        P3choix=Fwcol[id]
                        ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVchoix)*(Fwcol == P3choix))           
                        ax[itypeh].plot(x[ind],vn[var][ind],label = band_name[band],color=color[band],linewidth=lw)#,s=0.5)  
                    else:
                        ind = np.where((Tcol == T_dict[typeh])*(ELEV == ELEVchoix)*(Fwcol == Fwchoix))           
                        ax[itypeh].plot(x[ind],vn[var][ind],label = band_name[band],color=color[band],linewidth=lw)#,s=0.5)
                        if (plotR):
                            ax[itypeh].plot(x[ind],vn_R[var][ind],label = band_name[band]+', Rayleigh',color=color[band],linewidth=lw,ls='--')
                    
                #end choice wg
            # end loop over bands
            title1= typeName[typeh]+" T="+str(int(T_dict[typeh]))+u'\u00B0C'+'\n'
            title2= r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0'
            title3=" AR="+str(ARcnst)+" elev="+str(ELEVchoix)
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
        
        figtitle=var+"("+unit[var]+") as a function of "+plot+"("+pltunit[plot]+") - ELEV="+str(ELEVchoix)+u'\u00B0'+" - "+micro
        fig.suptitle(figtitle,fontsize=pol_suptitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        bandtit=''
        for band in band_list:
            bandtit+=band
        
        micron=micro
        if (plot=='D'):
            micron=''

        nomfig=dir_fig +'DistTmat'+bandtit+micron+'_'+plot+var
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
  

    
