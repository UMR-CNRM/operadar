# -*- coding: utf-8 -*-


import itertools
from pathlib import Path
from .utils_plot import *
import matplotlib.pyplot as plt
import matplotlib as mpl ; mpl.use('Agg')


  
def sensitivity_test(hydrometeor:str,
                     axeX:str,
                     dictParam:dict,
                     which_dpolVar:list,
                     method:str,
                     band:str,
                     folder_tables:str,
                     folder_figures:str,
                     ref:bool,
                     combine:list):
    # Checkig every scenario
    nrows_plot, asRow, asLegend = analyse_dict(dictParam,hydrometeor,combine) 
    # Plot creation
    fig,axes = plt.subplots(nrows=nrows_plot,
                            ncols=len(which_dpolVar),
                            figsize=(len(which_dpolVar)*4,nrows_plot*4),
                            layout='constrained',
                            )
    # Generic table name
    tableName, delim = get_table_name(axeX,band,hydrometeor)
    
    # SCENARIO 1 : only one value given per argument
    if asRow == None and asLegend == None :
        table_path = f"{folder_tables}{hydrometeor}/{dictParam['axis_ratio_func'][0]}_" \
                    +f"AR{dictParam['axis_ratio'][0]}_CANT{dictParam['canting_angle'][0]}_" \
                    +f"{dictParam['density_func'][0]}_Frim{dictParam['riming_fraction'][0]}_"\
                    +f"{dictParam['diel_func'][0]}/{tableName}" ; print(table_path)
        plot_table(table_path=table_path,
                   h=hydrometeor,
                   delim=delim,
                   color='tab:blue',
                   axeX=axeX,
                   which_dpolVar=which_dpolVar,
                   figAx=axes,
                   method=method,
                   legend=asLegend,
                   subfigRow=asRow,
                   Fw=dictParam['liquid_water_fraction'][0],
                   )
        if ref :
            table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_" \
                        +f"AR{base_cfg[hydrometeor]['ARvalue']}_" \
                        +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_"\
                        +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                        +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                        +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}" ; print(table_path)
            plot_table(table_path=table_path,
                       h=hydrometeor,
                       delim=delim,
                       color='dimgray',
                       axeX=axeX,
                       which_dpolVar=which_dpolVar,
                       figAx=axes,
                       method=method,
                       legend='ref',
                       subfigRow='ref',
                       Fw=dictParam['liquid_water_fraction']
                       )
        dict4naming = dictParam.copy()
        
    # SCENARIO 2 : one argument contains multiple values, other arguments have only one value 
    # ---> the varying parameter is overlaid on each subfigure and a legend is added
    elif nrows_plot==1 :
        
        combinations = get_combinations(dictParam,
                                        folder_tables,
                                        tableName,
                                        hydrometeor,
                                        )
        colors = get_colors(dictParam, asLegend, combinations, folder_tables, tableName, hydrometeor)
        for idxLegend, (arf, ar, cant, dsty, Frim, diel, Fw) in enumerate(combinations):
            table_path = f"{folder_tables}{hydrometeor}/{arf}_AR{ar}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"
            print(table_path)
            plot_table(table_path=table_path,
                       h=hydrometeor,
                       delim=delim,
                       color=colors[idxLegend],
                       axeX=axeX,
                       which_dpolVar=which_dpolVar,
                       figAx=axes,
                       method=method,
                       legend=asLegend,
                       subfigRow=asRow,
                       Fw=Fw,
                       )
        if ref :
            table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_" \
                        +f"AR{base_cfg[hydrometeor]['ARvalue']}_" \
                        +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_" \
                        +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                        +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                        +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}" ; print(table_path)
            plot_table(table_path=table_path,
                       h=hydrometeor,
                       delim=delim,
                       color='dimgray',
                       axeX=axeX,
                       which_dpolVar=which_dpolVar,
                       figAx=axes,
                       method=method,
                       legend='ref',
                       subfigRow='ref',
                       Fw=Fw,
                       )
        dict4naming = dictParam.copy()
        if type(asLegend) is str : dict4naming[asLegend] = ['Multi']
        elif type(asLegend) is list : 
            for lgd_elmt in asLegend : dict4naming[lgd_elmt] = ['Combined']
        
    # SCENARIO 3 :
    # --->
    elif nrows_plot>1 :
        dictCopy = dictParam.copy()
        for idxRow,col in enumerate(dictParam[asRow]) :
            dictCopy[asRow] = [col]
            combinations = get_combinations(dictCopy,
                                            folder_tables,
                                            tableName,
                                            hydrometeor,
                                            )
            colors = get_colors(dictParam, asLegend, combinations, folder_tables, tableName, hydrometeor)
            for idxLegend, (arf, ar, cant, dsty, Frim, diel, Fw) in enumerate(combinations):
                table_path = f"{folder_tables}{hydrometeor}/{arf}_AR{ar}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"
                print(table_path)
                plot_table(table_path=table_path,
                           h=hydrometeor,
                           delim=delim,
                           color=colors[idxLegend],
                           axeX=axeX,
                           which_dpolVar=which_dpolVar,
                           figAx=axes[idxRow,:],
                           method=method,
                           legend=asLegend,
                           subfigRow=asRow,
                           Fw=Fw,
                           )
            if ref :
                table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_"\
                            +f"AR{base_cfg[hydrometeor]['ARvalue']}_"\
                            +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_" \
                            +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                            +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                            +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}" ; print(table_path)
                plot_table(table_path=table_path,
                           h=hydrometeor,
                           delim=delim,
                           color='dimgray',
                           axeX=axeX,
                           which_dpolVar=which_dpolVar,
                           figAx=axes[idxRow,:],
                           method=method,
                           legend='ref',
                           subfigRow='ref',
                           Fw=Fw,
                           )
        dict4naming = dictParam.copy()
        if type(asLegend) is str : dict4naming[asLegend] = ['Multi']
        elif type(asLegend) is list : 
            for lgd_elmt in asLegend : dict4naming[lgd_elmt] = ['Combined']
        dict4naming[asRow] = ['Multi']
    # Add a grid
    axes = axes.flat
    for ax in axes:
        ax.grid(color='darkgray')
        
    # Handle legend
    if asLegend != None :
        ncol = len(combinations)
        if ref : ncol += 1
        if ncol > 5 : ncol = int(np.ceil(ncol/2))
        handles, labels = axes[-1].get_legend_handles_labels()
        fig.legend(handles=handles,
                  labels=labels,
                  loc='upper center',
                  bbox_to_anchor=(0.5, 0),
                  ncol=ncol,
                  )
    else :
        handles, labels = axes[-1].get_legend_handles_labels()
        fig.legend(handles=handles, labels=labels, loc='upper center',
                   bbox_to_anchor=(0.5, 0), ncol=2,
                )
    # Handle title
    fig_title = nameHydrometeors[hydrometeor]+" at T="+str(int(T_dict[hydrometeor]))+u'\u00B0C'
    fig.suptitle(fig_title ,fontsize=15)
    
    # Handle saving
    subfolder = f"{folder_figures}{hydrometeor}_{axeX}xaxis_{band}band/"
    Path(subfolder).mkdir(exist_ok=True, parents=True)
    plotName = f"ARfunc{dict4naming['axis_ratio_func'][0]}_"\
              +f"ARvalue{dict4naming['axis_ratio'][0]}_"\
              +f"Fw{dict4naming['liquid_water_fraction'][0]}_"\
              +f"CANTING{dict4naming['canting_angle'][0]}_"\
              +f"DSTY{dict4naming['density_func'][0]}_"\
              +f"Frim{dict4naming['riming_fraction'][0]}_" \
              +f"DIEL{dict4naming['diel_func'][0]}_"\
              +f"method{method}"
    if ref : plotName += '_withRef'
    fig.savefig(f"{subfolder}{plotName}.png",dpi=120, bbox_inches='tight')