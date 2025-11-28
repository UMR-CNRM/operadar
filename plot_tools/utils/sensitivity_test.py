# -*- coding: utf-8 -*-


from pathlib import Path
from .utils_plot import *
import matplotlib.pyplot as plt
import matplotlib as mpl ; mpl.use('Agg')


  
def sensitivity_test(hydrometeor:str,
                     axeX:str,
                     dictParam:dict,
                     which_dpolVar:list,
                     method:str,
                     microphysics:str,
                     moment:int,
                     band:str,
                     folder_tables:str,
                     folder_figures:str,
                     ref:bool,
                     combine:list,
                     invertColLegend:bool):
    # Checkig every scenario
    nrows_plot, asRow, asLegend = analyse_dict(dictParam,hydrometeor,combine,invertColLegend)
    print('Plot rows (nrows):',asRow,f'({nrows_plot})','\nPlot legend:',asLegend)
    # Plot creation
    fig,axes = plt.subplots(nrows=nrows_plot,
                            ncols=len(which_dpolVar),
                            figsize=(len(which_dpolVar)*3.2,nrows_plot*4),
                            layout='constrained',
                            )
    # Generic table name
    tableName, delim = get_table_name(axeX,band,hydrometeor,microphysics,moment)
    
    # SCENARIO 1 : only one value given per argument
    if asRow == None and asLegend == None :
        table_path = f"{folder_tables}{hydrometeor}/{dictParam['axis_ratio_func'][0]}_" \
                    +f"AR{dictParam['axis_ratio'][0]:.2f}_CANT{dictParam['canting_angle'][0]}_" \
                    +f"{dictParam['density_func'][0]}_Frim{dictParam['riming_fraction'][0]}_"\
                    +f"{dictParam['diel_func'][0]}/{tableName}"
        dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                    delim=delim,
                                                                    axeX=axeX,
                                                                    method=method,
                                                                    )
        Fw = dictParam['liquid_water_fraction'][0]
        Nc = dictParam['number_concentration'][0]
        plot_table(h=hydrometeor, param=dictSettings,
                    color='tab:blue', axeX=axeX, moment=moment,
                    which_dpolVar=which_dpolVar,
                    figAx=axes, method=method,
                    legend=asLegend, subfigRow=asRow,
                    Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                    ELEVcol=ELEVcol, dpolDict=dpolDict, xaxis=xaxis
                    )
        additional_info = add_common_parameterization(param=dictSettings,
                                                      axeX=axeX, hydro=hydrometeor,
                                                      whichLegend=asLegend,
                                                      whichRow=asRow,
                                                      Fw=Fw, Nc=Nc)
        if ref :
            table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_" \
                        +f"AR{base_cfg[hydrometeor]['ARvalue']:.2f}_" \
                        +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_"\
                        +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                        +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                        +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}"
            dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                        delim=delim,
                                                                        axeX=axeX,
                                                                        method=method,
                                                                        )
            plot_table(h=hydrometeor, param=dictSettings,
                        color='k', axeX=axeX, moment=moment,
                        which_dpolVar=which_dpolVar,
                        figAx=axes, method=method,
                        legend=asLegend, subfigRow=asRow,
                        Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                        ELEVcol=ELEVcol, dpolDict=dpolDict,
                        xaxis=xaxis, isRef=True,
                        )
            ref_info = add_common_parameterization(param=dictSettings, isRef=True,
                                                    axeX=axeX, hydro=hydrometeor,
                                                    whichLegend=asLegend,
                                                    whichRow=asRow,
                                                    Fw=Fw, Nc=Nc)
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
        count=0 
        for (arf, ar, cant, dsty, Frim, diel) in combinations:
            table_path = f"{folder_tables}{hydrometeor}/{arf}_AR{ar:.2f}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"
            dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                        delim=delim,
                                                                        axeX=axeX,
                                                                        method=method,
                                                                        )
            for Fw in dictParam['liquid_water_fraction'] :
                for Nc in dictParam['number_concentration'] :
                    plot_table(h=hydrometeor, param=dictSettings,
                                color=colors[count], axeX=axeX,
                                which_dpolVar=which_dpolVar,
                                figAx=axes, method=method, moment=moment,
                                legend=asLegend, subfigRow=asRow,
                                Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                                ELEVcol=ELEVcol, dpolDict=dpolDict, xaxis=xaxis,
                                )
                    count += 1
            additional_info = add_common_parameterization(param=dictSettings,
                                                          axeX=axeX, hydro=hydrometeor,
                                                          whichLegend=asLegend,
                                                          whichRow=asRow,
                                                          Fw=Fw, Nc=Nc)
        if ref :
            table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_" \
                        +f"AR{base_cfg[hydrometeor]['ARvalue']:.2f}_" \
                        +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_" \
                        +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                        +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                        +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}"
            dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                        delim=delim,
                                                                        axeX=axeX,
                                                                        method=method,
                                                                        )
            count_ref = 0
            for Fw in dictParam['liquid_water_fraction'] :
                for Nc in dictParam['number_concentration'] :
                    nbRef = len(dictParam['liquid_water_fraction'])*len(dictParam['number_concentration'])
                    plot_table(h=hydrometeor, param=dictSettings,
                                color=colors[count_ref], axeX=axeX,
                                which_dpolVar=which_dpolVar,
                                figAx=axes, method=method, moment=moment,
                                legend=asLegend, subfigRow=asRow,
                                Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                                ELEVcol=ELEVcol, dpolDict=dpolDict,
                                xaxis=xaxis, isRef=True, nbRef=nbRef,
                                )
                    count_ref += 1
                ref_info = add_common_parameterization(param=dictSettings, isRef=True,
                                                        axeX=axeX, hydro=hydrometeor,
                                                        whichLegend=asLegend,
                                                        whichRow=asRow,
                                                        Fw=Fw, Nc=Nc)
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
            count = 0
            for (arf, ar, cant, dsty, Frim, diel) in combinations:
                table_path = f"{folder_tables}{hydrometeor}/{arf}_AR{ar:.2f}_CANT{cant}_{dsty}_Frim{Frim}_{diel}/{tableName}"
                dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                            delim=delim,
                                                                            axeX=axeX,
                                                                            method=method,
                                                                            )
                for Fw in dictParam['liquid_water_fraction'] :
                    for Nc in dictParam['number_concentration'] :
                        plot_table(h=hydrometeor, param=dictSettings,
                                    color=colors[count], axeX=axeX,
                                    which_dpolVar=which_dpolVar, moment=moment,
                                    figAx=axes[idxRow,:], method=method,
                                    legend=asLegend, subfigRow=asRow,
                                    Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                                    ELEVcol=ELEVcol, dpolDict=dpolDict, xaxis=xaxis,
                                    )
                        count += 1
                additional_info = add_common_parameterization(param=dictSettings,
                                                              axeX=axeX, hydro=hydrometeor,
                                                              whichLegend=asLegend,
                                                              whichRow=asRow,
                                                              Fw=Fw, Nc=Nc)
        if ref :
            table_path = f"{folder_tables}{hydrometeor}/{base_cfg[hydrometeor]['ARfunc']}_" \
                        +f"AR{base_cfg[hydrometeor]['ARvalue']:.2f}_" \
                        +f"CANT{base_cfg[hydrometeor]['SIGBETA']}_" \
                        +f"{base_cfg[hydrometeor]['DSTYfunc']}_" \
                        +f"Frim{base_cfg[hydrometeor]['Frim']}_" \
                        +f"{base_cfg[hydrometeor]['DIEL']}/{tableName}"
            dictSettings,Tcol,ELEVcol,P3col,dpolDict,xaxis = read_table(table_path=table_path,
                                                                        delim=delim,
                                                                        axeX=axeX,
                                                                        method=method,
                                                                        )
            for idxRow,col in enumerate(dictParam[asRow]) :
                count_ref = 0
                for Fw in dictParam['liquid_water_fraction'] :
                    for Nc in dictParam['number_concentration'] :
                        nbRef = len(dictParam['liquid_water_fraction'])*len(dictParam['number_concentration'])
                        plot_table(h=hydrometeor, param=dictSettings,
                                    color=colors[count_ref], axeX=axeX,
                                    which_dpolVar=which_dpolVar, moment=moment,
                                    figAx=axes[idxRow,:], method=method,
                                    legend=asLegend, subfigRow=asRow,
                                    Fw=Fw, Nc=Nc, Tcol=Tcol, P3col=P3col,
                                    ELEVcol=ELEVcol, dpolDict=dpolDict,
                                    xaxis=xaxis, isRef=True, nbRef=nbRef,
                                    )
                        count_ref += 1
                ref_info = add_common_parameterization(param=dictSettings, isRef=True,
                                                        axeX=axeX, hydro=hydrometeor,
                                                        whichLegend=asLegend,
                                                        whichRow=asRow,
                                                        Fw=Fw, Nc=Nc)
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
        ncol = count
        if ref : ncol += count_ref
        if ncol > 6 : ncol = int(np.ceil(ncol/2))
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
    if axeX == 'M' : fig_title += f' (with {microphysics} microphysics)'
    fig.suptitle(fig_title ,fontsize=15)
    
    #Handle additional info about parametrizations used
    if ref :
        if asRow is not None :
            additional_info = additional_info + '\n'*6+ref_info
        else :
            additional_info = additional_info + '\n'*3+ref_info
    fig.text(1.01,0.5,additional_info,fontsize=12,ha='left',va='center')
    
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
    if axeX == 'M' :
        plotName += f"_{microphysics}_Nc{dict4naming['number_concentration'][0]}"
    if ref : plotName += '_withRef'
    fig.savefig(f"{subfolder}{plotName}.png",dpi=110, bbox_inches='tight')
    print('Saved :',f"{subfolder}{plotName}.png")
    fig.savefig(f"{subfolder}{plotName}.pdf",dpi=110, format='pdf', bbox_inches='tight')
    print('Saved :',f"{subfolder}{plotName}.pdf")
