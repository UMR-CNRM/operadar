#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:50:19 2026

@author: augros
"""

import os
import sys
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import glob
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ------------------------------------------------------------------
# Variable metadata (global so that plot_map can see it)
# ------------------------------------------------------------------
VAR_DICT = {
    "Zh": {"name": "Reflectivity", "min": 8, "max": 64, "step": 4, "unit": "dBZ"},
    "Zh_att": {"name": "Attenuated Reflectivity", "min": 8, "max": 64, "step": 4, "unit": "dBZ"},
    "Zdr": {"name": "Differential Reflectivity", "min": 0, "max": 6, "step": 0.5, "unit": "dB"},
    "Kdp": {"name": "Specific Differential Phase", "min": 0, "max": 6, "step": 0.5, "unit": "°/km"},
}

# ------------------------------------------------------------------
# Mapping routine
# ------------------------------------------------------------------
def plot_map(var, lon, lat, data, bounds, cmap, lev,
             out_dir, fname_base, model, micro, time,
             lon_min, lon_max, lat_min, lat_max):
    """
    Plot a 2-D filled-contour map of *data* and save to PNG.
    """
    fig = plt.figure(figsize=(13, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    cont=ax.contourf(lon, lat, data, levels=bounds, cmap=cmap, extend="both")

    # colorbar
    cbar = plt.colorbar(cont, orientation="vertical",
                        ticks=bounds, shrink=0.7, pad=0.02)
    cbar.set_label(f"{var} ({VAR_DICT[var]['unit']})", fontsize=16)

    # map decorations
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":", edgecolor="gray")
    ax.add_feature(cfeature.STATES.with_scale("10m"),
                   linewidth=1, linestyle="-", edgecolor="gray")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = gl.right_labels = False

    # title & output file name
    day = pd.to_datetime(time.values).strftime("%Y%m%d")
    hour = pd.to_datetime(time.values).strftime("%H%M")
    if lev == -1:
        title = f"Max {VAR_DICT[var]['name']} – {model} {micro} – {day} {hour} UTC"
        png_name = fname_base.replace("dpolvar", var).replace(".nc", "_max.png")
    else:
        title = f"{VAR_DICT[var]['name']} at level {lev} – {model} {micro} – {day} {hour} UTC"
        png_name = fname_base.replace("dpolvar", var).replace(".nc", f"_lev{lev}.png")

    plt.title(title, fontsize=16)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.savefig(os.path.join(out_dir, png_name),
                bbox_inches="tight", pad_inches=0, dpi=100)
    plt.close(fig)
