# -*- coding: utf-8 -*-
import os, sys
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl



# ====== BASIC DEFINITIONS ====== #
pol_label=10
pol_legend=10
pol_title=14
pol_suptitle=12
lw=1

pltX={'D':'Deq',
      'M':'M',
      }
pltunit={'D':'mm',
         'M':r'kg m$^{-3}$',
         }
dmax_dict={'ii':5,
           'ss':20,
           'gg':50,
           'tt':2,
           'rr':10,
           'wg':50,
           'hh':100,
           'wh':100,
           }
T_dict = {'ii':-30,
          'ss':-10,
          'gg':0,
          'tt':5,
          'rr':10,
          'wg':10,
          'hh':1,
          'wh':10,
          }
nameHydrometeors = {'ii':'Pristine ice',
                    'ss':'Dry Snow',
                    'gg':'Dry Graupel',
                    'tt':'Cloud Water',
                    'rr':'Rain',
                    'wg':'Wet Graupel',
                    'hh':'Dry Hail',
                    'wh':'Wet Hail',
                    }
unit={'Zh':'dBZ',
        'Zdr':'dB',
        'Kdp':u'\u00B0'+r' km$^{-1}$',
        'Rhohv':'/',
        }
base_cfg = {'ii':{'SIGBETA':0.0,
                  'ARfunc':'CNST',
                  'ARvalue':1.0,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'RY19dry',
                  },
            'ss':{'SIGBETA':0.0,
                  'ARfunc':'AUds',
                  'ARvalue':0.8,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'RY19dry',
                  },
            'gg':{'SIGBETA':0.0,
                  'ARfunc':'RYdg',
                  'ARvalue':0.8,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'RY19dry',
                  },
            'wg':{'SIGBETA':0.0,
                  'ARfunc':'RYwg',
                  'ARvalue':0.8,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'LBwetgr',
                  },
            'rr':{'SIGBETA':0.0,
                  'ARfunc':'BR02',
                  'ARvalue':1.0,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'Liebe91',
                  },
            'tt':{'SIGBETA':0.0,
                  'ARfunc':'CNST',
                  'ARvalue':1.0,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'Liebe91',
                  },
            'hh':{'SIGBETA':0.0,
                  'ARfunc':'RYdg',
                  'ARvalue':0.8,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'RY19dry',
                  },
            'wh':{'SIGBETA':0.0,
                  'ARfunc':'RYwg',
                  'ARvalue':0.8,
                  'DSTYfunc':'RHOX',
                  'Frim':1.0,
                  'DIEL':'LBwetgr',
                  },
            }
ls={'Tmatrix':'-',
    'Rayleigh':'--',
    }



# ====== FUNCTIONS ====== #
def get_Y_min_max(var,hydro,xAxis):
    ymin_dict,ymax_dict={},{}
    for dpol in ['Zh','Zdr','Kdp','Rhohv']: ymin_dict[dpol],ymax_dict[dpol]={'D':{},'M':{}},{'D':{},'M':{}}
    # ZH
    ymin_dict["Zh"]['M'] = {'ii':-50,'ss':-60,'gg':-40,'tt':-60,'rr':-20,'wg':-10,'hh': 20,'wh': 20}
    ymax_dict["Zh"]['M'] = {'ii': 30,'ss': 30,'gg': 80,'tt': 20,'rr': 70,'wg': 80,'hh':100,'wh':100}
    ymin_dict["Zh"]['D'] = {'ii':-50,'ss':-60,'gg':-40,'tt':-60,'rr':-20,'wg':-10,'hh': 20,'wh': 20}
    ymax_dict["Zh"]['D'] = {'ii': 30,'ss': 30,'gg': 80,'tt': 20,'rr': 70,'wg': 80,'hh':100,'wh':100}
    # ZDR
    ymin_dict["Zdr"]['M'] = {'ii':0,'ss':-2,'gg':-2,'tt':0,'rr':-4,'wg':-4,'hh':-2,'wh':-4}
    ymax_dict["Zdr"]['M'] = {'ii':6,'ss':2,'gg':2,'tt':1,'rr':10,'wg':10,'hh':2,'wh':10}
    ymin_dict["Zdr"]['D'] = {'ii':0,'ss':0,'gg':-2,'tt':0,'rr':-4,'wg':-4,'hh':-2,'wh':-4}
    ymax_dict["Zdr"]['D'] = {'ii':6,'ss':0.5,'gg':2,'tt':1,'rr':10,'wg':10,'hh':2,'wh':10}
    # KDP
    ymin_dict["Kdp"]['M'] = {'ii':0,'ss':0,'gg':-0.2,'tt':0,'rr':-2,'wg':-4,'hh':-30,'wh':-30}
    ymax_dict["Kdp"]['M'] = {'ii':0.2,'ss':0.2,'gg':0.2,'tt':1,'rr':5,'wg':4,'hh':20,'wh':20}
    ymin_dict["Kdp"]['D'] = {'ii':0,'ss':0,'gg':-0.2,'tt':0,'rr':-2,'wg':-4,'hh':-30,'wh':-30}
    ymax_dict["Kdp"]['D'] = {'ii':0.05,'ss':0.2,'gg':0.2,'tt':1,'rr':5,'wg':4,'hh':20,'wh':20}
    # RHOHV
    ymin_dict["Rhohv"]['M'][hydro]=0.9
    ymax_dict["Rhohv"]['M'][hydro]=1.0
    ymin_dict["Rhohv"]['D'][hydro]=0.9
    ymax_dict["Rhohv"]['D'][hydro]=1.0
    
    return ymin_dict[var][xAxis][hydro], ymax_dict[var][xAxis][hydro]



def get_table_name(axeX,band,hydro,micro='ICE3'):
    tableName={'D':f'TmatCoefDiff_{band}{hydro}',
               'M':f'TmatCoefInt_{micro}_{band}{hydro}',
               }
    delim={'D':r"\s+",
           'M':';',
           }
    return tableName[axeX], delim[axeX]



def read_table(table_path,delim,axeX,method):
    dpol={'Tmatrix':{},'Rayleigh':{}}
    
    # reading table options (from config file)
    df_param=pd.read_csv(table_path, sep=delim,nrows=1,engine='python')
    SIGBETA=int(df_param['SIGBETA'][0])
    ARfunc=df_param['ARfunc'][0]
    ARcnst=df_param['ARcnst'][0]
    DSTYfunc=df_param['DSTYfunc'][0]
    DIEL=df_param['DIEL'][0]
    Frim=df_param['Frim'][0]
    
    # reading table variables
    df=pd.read_csv(table_path, sep=delim,skiprows =2,engine='python')
    Tcol=df['Tc'].to_numpy()
    ELEV=df['ELEV'].to_numpy()
    if (axeX=='M'): Fwcol=df['P3'].to_numpy()                            
    else: Fwcol=df['Fw'].to_numpy()

    if method == 'Tmatrix':
        dpol[method]['Zh'],dpol[method]['Zdr']=df['zhh'].to_numpy(),df['zdr'].to_numpy()
        dpol[method]['Rhohv'],dpol[method]['Kdp']=df['rhohv'].to_numpy(),df['kdp'].to_numpy()
    elif method == 'Rayleigh' :
        dpol[method]['Zh'],dpol[method]['Zdr']=df['zhhR'].to_numpy(),df['zdrR'].to_numpy()
        dpol[method]['Rhohv'],dpol[method]['Kdp']=df['rhohvR'].to_numpy(),df['kdpR'].to_numpy()
    elif method == 'Both' :
        dpol['Tmatrix']['Zh'],dpol['Tmatrix']['Zdr']=df['zhh'].to_numpy(),df['zdr'].to_numpy()
        dpol['Tmatrix']['Rhohv'],dpol['Tmatrix']['Kdp']=df['rhohv'].to_numpy(),df['kdp'].to_numpy()
        dpol['Rayleigh']['Zh'],dpol['Rayleigh']['Zdr']=df['zhhR'].to_numpy(),df['zdrR'].to_numpy()
        dpol['Rayleigh']['Rhohv'],dpol['Rayleigh']['Kdp']=df['rhohvR'].to_numpy(),df['kdpR'].to_numpy()
        
    xaxis=df[pltX[axeX]].to_numpy()
    return SIGBETA, DIEL, ARfunc, ARcnst, DSTYfunc, Frim, Tcol, ELEV, Fwcol, dpol,xaxis



def set_legend(whichLegend,SIGBETA,ARfunc,ARcnst,DSTYfunc,Frim,DIELfunc,Fw):
    if type(whichLegend)==str :
        whichLegend = [whichLegend]
    label = ''
    for lgd_elmt in whichLegend :
        if lgd_elmt == 'canting_angle' :
            label += r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0'+' '
        elif lgd_elmt == 'axis_ratio_func' :
            label += f'ARfunc:{ARfunc} '
        elif lgd_elmt == 'axis_ratio' :
            label += f'AR={ARcnst} '
        elif lgd_elmt == 'liquid_water_fraction' :
            label += f'Fw={Fw} '
        elif lgd_elmt == 'density_func' :
            label += f'DSTYfunc:{DSTYfunc} '
        elif lgd_elmt == 'riming_fraction' :
            label += f'frim:{Frim} '
        elif lgd_elmt == 'diel_func' :
            label += f'DIELfunc:{DSTYfunc} '
        elif lgd_elmt == 'ref' :
            label = f'Ref : AR={ARcnst} ({ARfunc}) Fw={Fw} '+r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0' \
                    + '\n' + f'DSTYfunc: {DSTYfunc} DIELfunc: {DIELfunc}'
    return label


def set_subfig_title(whichRow,SIGBETA,ARfunc,ARcnst,DSTYfunc,Frim,DIELfunc,whichLegend,Fw):
    subfig_title = {'canting_angle' : r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0',
                    'axis_ratio_func' : f'ARfunc: {ARfunc}',
                    'axis_ratio' : f'AR={ARcnst}',
                    'density_func' : f'DSTYfunc: {DSTYfunc}',
                    'riming_fraction' : f'frim={Frim}',
                    'diel_func': f'DIELfunc: {DIELfunc}',
                    'liquid_water_fraction':f'Fw={Fw}',
                    'ref' : 'Reference plot',
                    }
    if whichRow == None :
        subfig_title.pop('ref')
        rowTitle = ''
        if whichLegend != None :
            if type(whichLegend) is str :
                subfig_title.pop(whichLegend)
            elif type(whichLegend) is list :
                for lgd_elmt in whichLegend :
                    subfig_title.pop(lgd_elmt)
        for subfigElement in subfig_title.values() :
            rowTitle += '\n'+subfigElement
    else :
        rowTitle = subfig_title[whichRow]
        
    return rowTitle



def plot_table(table_path:str, h:str,
               delim:str, axeX:str,
               which_dpolVar:list, figAx:str,
               method:str,
               color:str,
               legend:str|None,
               subfigRow:str|None,
               Fw:float,
               ):
    SIGBETA,DIEL,ARfunc,ARcnst,DSTYfunc,Frim,Tcol,ELEV,Fwcol,dpolDict,xaxis = read_table(table_path,delim,axeX,method)
    if (    base_cfg[h]['SIGBETA'] == SIGBETA 
        and base_cfg[h]['ARfunc'] == ARfunc
        and base_cfg[h]['ARvalue'] == ARcnst 
        and base_cfg[h]['DSTYfunc'] == DSTYfunc  
        and base_cfg[h]['DIEL'] == DIEL
        and base_cfg[h]['Frim'] == Frim
        ):
        if h[0:1] == 'w' and subfigRow == None :
            legendItem = legend
        else :
            legendItem = 'ref'
            color='dimgray'
    else :
        legendItem = legend
    
    row = np.where((Tcol == T_dict[h])*(ELEV == 0)*(Fwcol == Fw))
    lgd_label = set_legend(whichLegend=legendItem,
                           SIGBETA=SIGBETA,
                           ARfunc=ARfunc,
                           ARcnst=ARcnst,
                           DSTYfunc=DSTYfunc,
                           Frim=Frim,
                           DIELfunc=DIEL,
                           Fw=Fw,
                           )
    subfig_title = set_subfig_title(whichRow=subfigRow,
                                    SIGBETA=SIGBETA,
                                    ARfunc=ARfunc,
                                    ARcnst=ARcnst,
                                    DSTYfunc=DSTYfunc,
                                    Frim=Frim,
                                    DIELfunc=DIEL,
                                    whichLegend=legend,
                                    Fw=Fw,
                                    )
    for ivar,varName in enumerate(which_dpolVar) :
        ymin,ymax=get_Y_min_max(varName,h,axeX)
        if method == 'Both': listMethod = ['Tmatrix','Rayleigh']
        else : listMethod = [method]
        for m in listMethod :
            yaxis = dpolDict[m][varName]
            figAx[ivar].plot(xaxis[row],yaxis[row],
                             label=lgd_label+f' {m}',
                             linewidth=lw,
                             color=color,
                             ls=ls[m],
                             )
        figAx[ivar].set_xlim(0,dmax_dict[h])
        figAx[ivar].set_ylim(ymin,ymax)
        figAx[ivar].set_ylabel(varName+"("+unit[varName]+")",fontsize = pol_label)
        figAx[ivar].set_xlabel(axeX+" ("+pltunit[axeX]+")",fontsize = pol_label)
        figAx[ivar].tick_params(axis='x', labelsize=pol_label)
        figAx[ivar].tick_params(axis='y', labelsize=pol_label)
        if subfigRow != 'ref' :
            figAx[ivar].set_title(subfig_title,fontsize=pol_suptitle)
        


def analyse_dict(dictParam:dict,hydrometeor:str,combine:list) :
    
    if len(combine)==1:
        print('No legend combination possible with combine =',combine) ; sys.exit()
    if len(dictParam['axis_ratio'])<1 :
        dictParam['axis_ratio'] = [base_cfg[hydrometeor]['ARvalue']]
    if len(dictParam['canting_angle'])<1 :
        dictParam['canting_angle'] += [base_cfg[hydrometeor]['SIGBETA']]
    if len(dictParam['axis_ratio_func'])<1 :
        dictParam['axis_ratio_func'] += [base_cfg[hydrometeor]['ARfunc']]
    if len(dictParam['density_func'])<1 :
        dictParam['density_func'] += [base_cfg[hydrometeor]['DSTYfunc']]
    if len(dictParam['riming_fraction'])<1 :
        dictParam['riming_fraction'] += [base_cfg[hydrometeor]['Frim']]
    if len(dictParam['diel_func'])<1 :
        dictParam['diel_func'] += [base_cfg[hydrometeor]['DIEL']]
    if len(dictParam['liquid_water_fraction'])<1 :
        dictParam['liquid_water_fraction'] += [0.0]
    
    nb_keys = len(dictParam.keys())
    nb_items = np.sum([len(dictParam[key]) for key in dictParam.keys()])
        
    if nb_keys == nb_items :
        nrows = 1
        variation_legend = None
        variation_columns = None
    elif nb_items > nb_keys :
        nb_fixed_param = np.sum([len(item) for _,item in dictParam.items() if len(item)==1])
        if len(combine)>=2 : nb_fixed_param += 1
        if nb_keys-nb_fixed_param == 1 :
            nrows = 1
            if len(combine)>=2 :
                variation_legend = [key for key,item in dictParam.items() if len(item)>1 and key in combine]
                print("Will combine in legend :",variation_legend)
            else :
                variation_legend = [key for key,item in dictParam.items() if len(item)>1][0]
            variation_columns=None
        elif nb_keys-nb_fixed_param == 2 :
            if len(combine)>=2 :
                variation_legend = [key for key,item in dictParam.items() if len(item)>1 and key in combine]
                print("Will combine in legend :",variation_legend)
                variation_columns = [key for key,item in dictParam.items() if len(item)>1 and key not in combine][0]
            else :
                variation_legend = [key for key,item in dictParam.items() if len(item)>1][0]
                variation_columns = [key for key,item in dictParam.items() if len(item)>1][1]
            nrows = len(dictParam[variation_columns])
        else :
            print('')
            print('/!\ ERROR : Only 2 parameters can vary among :',dictParam.keys())
            print('    ------> You need to fix at least',nb_fixed_param-2,'of them.')
            sys.exit()
    #print(nrows,'rows:',variation_columns,'legend:',variation_legend) ; sys.exit()
    return nrows, variation_columns, variation_legend



def get_colors(dictParam:dict,
               asLegend:list|str,
               combinations:tuple,
               tableDir:str,
               tableName:str,
               hydrometeor:str,
               ):
    if type(asLegend) is str :
        nb_colors = len(dictParam[asLegend])
    elif type(asLegend) is list :
        nb_colors = 0
        for arf, ar, cant, dsty, Frim, diel, Fw in combinations :
            if os.path.isfile(f"{tableDir}{hydrometeor}/{arf}_AR{ar}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"):
                nb_colors += 1
    colors = [colorName for colorName in mpl.colors.TABLEAU_COLORS.keys()][0:nb_colors]
    return colors



def get_combinations(dictParam:dict,
                     tableDir:str,
                     tableName:str,
                     hydrometeor:str,
                     ) :
    combinations = list(itertools.product(dictParam['axis_ratio_func'],
                                          dictParam['axis_ratio'],
                                          dictParam['canting_angle'],
                                          dictParam['density_func'],
                                          dictParam['riming_fraction'],
                                          dictParam['diel_func'],
                                          dictParam['liquid_water_fraction']
                                          )
                        )
    final_combinations = []
    for arf, ar, cant, dsty, Frim, diel, Fw in combinations :
        filepath = f"{tableDir}{hydrometeor}/{arf}_AR{ar}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"
        if  os.path.isfile(filepath):
            final_combinations += [(arf, ar, cant, dsty, Frim, diel, Fw)]
    return final_combinations