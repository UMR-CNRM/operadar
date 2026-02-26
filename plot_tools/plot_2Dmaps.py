#!/usr/bin/env python3
"""
Create PNG maps of polarimetric radar variables from operadar NetCDF files.
For real case only (latitude and longitude are needed) => operadar configuration file needs
real_case=True


Example
-------
python3 plot_2Dmaps.py --dataDir ../modelFiles/AROME/20250831/ --vars Zh Zdr Kdp --levels 89 80 75
"""
import os
import sys
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import epygram
epygram.init_env()

# ------------------------------------------------------------------
# Variable metadata (global so that plot_map can see it)
# ------------------------------------------------------------------
VAR_DICT = {
    "Zh": {"name": "Reflectivity", "min": 8, "max": 64, "step": 4, "unit": "dBZ"},
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
    ax.contourf(lon, lat, data, levels=bounds, cmap=cmap, extend="both")

    # colorbar
    cbar = plt.colorbar(ax.collections[0], orientation="vertical",
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

# ------------------------------------------------------------------
# Processing loop
# ------------------------------------------------------------------
def main(data_dir, variables, levels):
    """Loop over NetCDF files and create maps."""

    # domain
    lon_min, lon_max, lat_min, lat_max = -5., 10., 41., 52.
    # SESAR domain ?
    #lon_min, lon_max, lat_min, lat_max = -12.0, 16., 37.5, 55.4

    # colormap
    epygram.util.load_cmap("radar")
    cmap = plt.get_cmap("radar")
    cmap.set_under("white")
    cmap.set_over("deeppink")

    for root, _, files in os.walk(data_dir):
        for fname in files:
            if not fname.lower().endswith(".nc"):
                continue
            fpath = os.path.join(root, fname)
            print(f"Processing: {fpath}")
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
                             root, fname, ds.model, ds.microphysics, ds.time,
                             lon_min, lon_max, lat_min, lat_max)

                # column maximum
                print(f"  Plot {var} max")
                fld_max = ds[var].max(dim="level")
                plot_map(var, fld_max.lon, fld_max.lat, fld_max, bounds,
                         cmap, -1, root, fname,
                         ds.model, ds.microphysics, ds.time,
                         lon_min, lon_max, lat_min, lat_max)

            ds.close()

# ------------------------------------------------------------------
# Command-line interface
# ------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate PNG maps of radar variables from AROME NetCDF files.")
    parser.add_argument("--dataDir", required=True,
                        help="Directory containing NetCDF files (walked recursively)")
    parser.add_argument("--vars", nargs="+", default=["Zh"],
                        help="Variables to plot (space-separated)")
    parser.add_argument("--levels", nargs="+", type=int,
                        default=[89],
                        help="Model levels to plot (space-separated)")
    args = parser.parse_args()

    if not os.path.isdir(args.dataDir):
        sys.exit(f"Error: directory '{args.dataDir}' does not exist.")

    main(args.dataDir, args.vars, args.levels)