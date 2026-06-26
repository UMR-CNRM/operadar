# -- coding: utf-8 --
"""
Created on Tue Dec 1 09:23:04 2020
@author: borderiesm and augrosc
Trace distributions D-Zhh et M-Zhh (ou Zdr, Kdp, Rhohv ou Ah)
pour bandes de frequence specifiee dans band_list
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
from utils.utils_plot import sample_cmap

# --- Configuration ---
micro = "ICE3"
TmatOption={'S':'improved','C':'vertical','K':'default','Ka':'default','Ku':'default','W':'default','L':'David2026PhD'}
moments = {
    "ICE3": {"rr": "1M", "ss": "1M", "gg": "1M", "wg": "1M", "cl": "1M", "cs": "1M", "ii": "1M"},
    "ICJW": {"rr": "1M", "ss": "1M", "gg": "1M", "wg": "1M", "cl": "1M", "cs": "1M", "ii": "1M"},
    "LIMA": {"rr": "2M", "ss": "1M", "gg": "1M", "wg": "1M", "cl": "2M", "cs": "2M", "ii": "2M"},
    "LIMC": {"rr": "2M", "ss": "1M", "gg": "1M", "wg": "1M", "cl": "2M", "cs": "2M", "ii": "1M"}
}
plotR = True
band_list = ['C', 'Ku','K', 'Ka', 'W']
typeh_list = ['rr','ss', 'gg', 'cl', 'ii', 'wg']
listplot=['M'] #'D','M']

Path_tables = "../tables_generator/tables/"
dir_fig = "IMG/"

# Bandes et types

pltunit={'D':'mm','M':r'kg m$^{-3}$'}
pltX = {'D': 'Deq', 'M': 'M'}
listvar = ['Zh']
unit={'Zh':'dBZ','Zdr':'dB','Kdp':u'\u00B0'+r' km$^{-1}$','Rhohv':'/','Ah':'dB'+r' km$^{-1}$','Av':'dB'+r' km$^{-1}$'}
typeName = {
    'ii': 'Pristine ice', 'ss': 'Dry Snow', 'gg': 'Dry Graupel', 'cl': 'Cloud Water',
    'cs': 'Cloud Water', 'rr': 'Rain', 'wg': 'Wet Graupel', 'hh': 'Dry Hail', 'wh': 'Wet Hail'
}

# Limites Y
ymin_dict, ymax_dict = {}, {}
for var in ['Zh', 'Zdr', 'Rhohv', 'Kdp', 'Ah', 'Av']:
    ymin_dict[var], ymax_dict[var] = {}, {}

ymin_dict["Zh"] = {'ii': -20, 'ss': -20, 'gg': -20, 'cl': -20, 'cs': -20, 'rr': -20, 'wg': -20, 'hh': -20, 'wh': -20}
ymax_dict["Zh"] = {'ii': 70, 'ss': 70, 'gg': 70, 'cl': 70, 'cs': 70, 'rr': 70, 'wg': 70, 'hh': 70, 'wh': 70}
ymin_dict["Zdr"] = {'ii': 0, 'ss': -2, 'gg': -2, 'cl': 0, 'cs': 0, 'rr': -4, 'wg': -4, 'hh': -2, 'wh': -4}
ymax_dict["Zdr"] = {'ii': 6, 'ss': 2, 'gg': 2, 'cl': 1, 'cs': 1, 'rr': 10, 'wg': 10, 'hh': 2, 'wh': 10}
ymin_dict["Kdp"] = {'ii': 0, 'ss': 0, 'gg': -0.2, 'cl': 0, 'cs': 0, 'rr': -2, 'wg': -4, 'hh': -30, 'wh': -30}
ymax_dict["Kdp"] = {'ii': 0.2, 'ss': 0.2, 'gg': 0.2, 'cl': 1, 'cs': 1, 'rr': 5, 'wg': 4, 'hh': 20, 'wh': 20}
ymin_dict["Ah"] = {'ii': 0, 'ss': 0, 'gg': 0, 'cl': 0, 'rr': 0, 'wg': 0, 'hh': 0, 'wh': 0}
ymax_dict["Ah"] = {'ii': 10, 'ss': 10, 'gg': 10, 'cl': 10, 'rr': 10, 'wg': 10, 'hh': 10, 'wh': 10}
ymin_dict["Av"] = {'ii': 0, 'ss': 0, 'gg': 0, 'cl': 0, 'rr': 0, 'wg': 0, 'hh': 0, 'wh': 0}
ymax_dict["Av"] = {'ii': 10, 'ss': 10, 'gg': 10, 'cl': 10, 'rr': 10, 'wg': 10, 'hh': 10, 'wh': 10}

dmax_dict = {'ii': 10, 'ss': 20, 'gg': 50, 'cl': 2, 'cs': 2, 'rr': 10, 'wg': 50, 'hh': 100, 'wh': 100}
Fwsel = 0
Fw_list, Fw_ls = [0.0, 0.1, 0.6, 1.0], ['-.', ':', '--', '-']
ELEVsel = 90
Nii = 800
expN_list, N_ls = [3], ['-.', '-', '--']
T_dict = {'ii': -30, 'ss': -10, 'gg': 0, 'cl': 10, 'cs': 10, 'rr': 10, 'wg': 10, 'hh': 1, 'wh': 10}
delim = {'D': r"\s+", 'M': ';'}

# Palette de couleurs
palette = sample_cmap('viridis', n=len(band_list), as_hex=True)
color = {lab: palette[i] for i, lab in enumerate(band_list)}

# Style de ligne par bande
style = {b: '-' for b in band_list}

# Boucle sur variables et plots
for var in listvar:
    for plot in listplot:
        print("Plotting " + var + "-" + plot)
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(24, 12))
        ax = axes.flat

        for itypeh, typeh in enumerate(typeh_list):
            print(f"type : {typeh}")
            begtable={'D':'TmatCoefDiff_','M':'TmatCoefInt_'+micro+'_'+moments[micro][typeh]+'_'}
            
            for iband, band in enumerate(band_list):
                print(f"band : {band}")
                table = f"{Path_tables}/{TmatOption[band]}/{begtable[plot]}{band}{typeh}"

                # Lecture des paramètres
                df_param = pd.read_csv(table, sep=delim[plot], nrows=1, engine='python')
                SIGBETA = int(df_param['SIGBETA'][0])
                ARfunc = df_param['ARfunc'][0]
                ARcnst = df_param['ARcnst'][0]

                # Lecture données
                df = pd.read_csv(table, sep=delim[plot], skiprows=2, engine='python')
                Tcol = df['Tc'].to_numpy()
                ELEV = df['ELEV'].to_numpy()

                if plot == 'M':
                    P3col = df['P3'].to_numpy()
                else:
                    P3col = df['Fw'].to_numpy()

                vn = {
                    'Zh': df['zhh'].to_numpy(),
                    'Zdr': df['zdr'].to_numpy(),
                    'Rhohv': df['rhohv'].to_numpy(),
                    'Kdp': df['kdp'].to_numpy(),
                    'Ah': df['Ah'].to_numpy(),
                    'Av': df['Av'].to_numpy()
                }
                vn[var]=np.array(vn[var], dtype=float)

                vn_R = {
                    'Zh': df['zhhR'].to_numpy(),
                    'Zdr': df['zdrR'].to_numpy(),
                    'Rhohv': df['rhohvR'].to_numpy(),
                    'Kdp': df['kdpR'].to_numpy(),
                    'Ah': df['AhR'].to_numpy(),
                    'Av': df['AvR'].to_numpy()
                }
                vn_R[var]=np.array(vn_R[var], dtype=float)

                x = df[pltX[plot]].to_numpy()

                # Gestion des cas spéciaux
                if typeh in ['wg', 'wh']:
                    for iFw, Fw in enumerate(Fw_list):
                        ind = np.where((Tcol == T_dict[typeh]) & (ELEV == ELEVsel) & (P3col == Fw))
                        label = f"Fw={Fw}"
                        ax[itypeh].plot(x[ind], vn[var][ind], label=label, color=color[band], ls=Fw_ls[iFw], linewidth=3)
                elif (typeh == 'ii' and plot == 'M') or (plot == 'M' and moments[micro][typeh] == "2M"):
                    expNcol = np.log10(P3col)
                    for iN, expN in enumerate(expN_list):
                        id = np.argmin(np.abs(expNcol - expN))
                        expNsel = expNcol[id]
                        ind = np.where((Tcol == T_dict[typeh]) & (ELEV == ELEVsel) & (expNcol == expNsel))
                        label = f" expN={expNsel:.1f}"
                        ax[itypeh].plot(x[ind], vn[var][ind], label=label, color=color[band], ls=N_ls[iN], linewidth=3)
                else:
                    ind = np.where((Tcol == T_dict[typeh]) & (ELEV == ELEVsel) & (P3col == Fwsel))
                    ax[itypeh].plot(x[ind], vn[var][ind], label=band, color=color[band], linewidth=3)
                    if plotR:
                        ax[itypeh].plot(x[ind], vn_R[var][ind], label=f"{band}, Rayleigh", color=color[band], ls='--', linewidth=3)

                # Titre et labels
                title1 = f"{typeName[typeh]} T={T_dict[typeh]}°C\n"
                title2 = r'$\sigma_{\beta}=$'+str(SIGBETA)+u'\u00B0'
                title3 = f" AR={ARcnst} elev={ELEVsel}°"
                if typeh == "rr":
                    title3 = f" AR={ARfunc}"
                ax[itypeh].set_title(title1 + title2 + title3, fontsize=22)
                ax[itypeh].set_ylabel(f"{var}({unit[var]})", fontsize=20)
                ax[itypeh].set_xlabel(f"{plot}({pltunit[plot]})", fontsize=20)
                ax[itypeh].tick_params(axis='x', labelsize=20)
                ax[itypeh].tick_params(axis='y', labelsize=20)
                ax[itypeh].set_ylim(ymin_dict[var][typeh], ymax_dict[var][typeh])
                if plot == "D":
                    ax[itypeh].set_xlim(0, dmax_dict[typeh])
                elif plot == "M":
                    ax[itypeh].set_xlim(1e-5, 1e-2)
                    ax[itypeh].set_xscale('log')

                if (iband==0):
                    if (typeh in ['wg', 'wh']):
                        ax[itypeh].legend(loc='best', fontsize=16)
                    if ((typeh == 'ii' and plot == 'M') or (plot == 'M' and moments[micro][typeh] == "2M")):
                        ax[itypeh].legend(loc='best', fontsize=16)
                
                    
                ax[itypeh].grid()

        # Légende commune sous les figures
        handles, labels = ax[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.02), ncol=len(band_list), fontsize=18)

        # Titre global
        figtitle = f"{var}({unit[var]}) as a function of {plot}({pltunit[plot]}) - ELEV={ELEVsel}° - {micro}"
        fig.suptitle(figtitle, fontsize=30)
        fig.tight_layout(rect=[0, 0.08, 1, 0.95])

        # Sauvegarde
        bandtit = ''.join(band_list)
        micron = '' if plot == 'D' else micro
        nomfig = f"{dir_fig}DistTmat{bandtit}{micron}_{plot}{var}{TmatOption[band_list[0]]}"
        if plotR:
            nomfig += "_R"

        if not os.path.exists(dir_fig):
            os.makedirs(dir_fig)

        fig.savefig(nomfig + ".png", dpi=200, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        print(f"Figure saved in: {nomfig}")
