#!/usr/bin/env python3
"""
Create PNG maps of polarimetric radar variables from operadar NetCDF files.
For real case only (latitude and longitude are needed) => operadar configuration file needs
real_case=True

Select the variables to plot, the directory containing operadar files, 
the output directory for the images to save and the domain 

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

import epygram
epygram.init_env()

from utils.plot_map import plot_map, VAR_DICT


# === Directories, variables and model levels to plot ===
dataDir="../modelFiles/AROME/"
outputDir="./IMG/"
variables=["Zh","Zdr"] #["Zh", "Zdr", "Kdp", "Ah", "Rhohv", "Zh_att"]
levels=[89] 


# === Domain ===
lon_min, lon_max, lat_min, lat_max = -5., 10., 41., 52.
# SESAR domain ?
#lon_min, lon_max, lat_min, lat_max = -12.0, 16., 37.5, 55.4

# === Colormap ===
epygram.util.load_cmap("radar")
cmap = plt.get_cmap("radar")
cmap.set_under("white")
cmap.set_over("deeppink")

# === Loop over variables to plot ===
for fpath in glob.glob(dataDir + "//*.nc"):
    if os.path.isfile(fpath):
        fname=os.path.basename(fpath)
        print(f"Processing: {fpath}")
        print(f"Filename: {fname}")
        ds = xr.open_dataset(fpath)

        for var in variables:
            if var not in VAR_DICT:
                print(f"  Warning: {var} not in VAR_DICT – skipped")
                continue
            vmin = VAR_DICT[var]["min"]
            vmax = VAR_DICT[var]["max"]
            step = VAR_DICT[var]["step"]
            bounds = np.arange(vmin, vmax + step, step)

            # individual levels
            for lev in levels:
                print(f"  Plot {var} at level {lev}")
                fld = ds[var].sel(level=lev)
                plot_map(var, fld.lon, fld.lat, fld, bounds, cmap, lev,
                         outputDir, fname, ds.model, ds.microphysics,
                         ds.time, lon_min, lon_max, lat_min, lat_max)

            # column maximum
            print(f"  Plot {var} max")
            fld_max = ds[var].max(dim="level")
            plot_map(var, fld_max.lon, fld_max.lat, fld_max, bounds,
                     cmap, -1, outputDir, fname,
                     ds.model, ds.microphysics, ds.time,
                     lon_min, lon_max, lat_min, lat_max)

        ds.close()