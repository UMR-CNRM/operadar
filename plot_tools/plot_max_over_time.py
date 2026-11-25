#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Maximum over several times of a radar variable (storm swaths), at a given height or over the column, with the
contour Zh = --contour dBZ of each time step (tracks of the cells). Works with operadar outputs (one file per
time) and with gridded files with a time dimension (e.g. observations). Replaces draft/plot_2Dcumulmaps.py.

Maps are drawn in lon/lat when the files have lat/lon (with coastlines if cartopy is installed), else in km.

Examples (from plot_tools/):
  python plot_max_over_time.py "../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20220818_0*.nc" --var Zh
  python plot_max_over_time.py "../operadarFiles/AROME/*.nc" --var Zdr --height 2 --contour 40
  python plot_max_over_time.py OBS/20220818_obs.nc --var Kdp --box 4 11 41 45
"""
import argparse
import os

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

from utils.operadar_output import (RADAR_CMAP, HeightGrid, axis_label, column_max, describe, expand_files,
                                   horizontal_km, open_output, time_steps, time_string, var_info, vertical_dim)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', help='operadar output files or gridded files (shell patterns allowed)')
    p.add_argument('--var', default='Zh', help='variable (default Zh)')
    p.add_argument('--height', type=float, default=None,
                   help='height in km (default: column maximum); above the lowest model level (agl)')
    p.add_argument('--contour', type=float, default=40., help='Zh contour of each time (dBZ; <= -100: none)')
    p.add_argument('--box', nargs=4, type=float, default=None, metavar=('LONMIN', 'LONMAX', 'LATMIN', 'LATMAX'))
    p.add_argument('--bounds', nargs=3, type=float, default=None, metavar=('MIN', 'MAX', 'STEP'),
                   help='colour bounds (default: from the variable)')
    p.add_argument('--out', default='IMG/maxtime/', help='output folder (default IMG/maxtime/)')
    p.add_argument('--dpi', type=int, default=120)
    return p.parse_args()


def field_2d(ds, name, height):
    if height is None:
        return column_max(ds[name], vertical_dim(ds))
    grid = HeightGrid(ds, [height * 1000.], 'agl')
    return grid(ds[name])[0].reshape(ds.sizes['y'], ds.sizes['x'])


def main():
    args = parse_args()
    files = expand_files(args.files)
    vmax = None
    times, coords, desc = [], None, ''
    zh_list = []
    for path in files:
        ds_file = open_output(path, args.box)
        for ds in time_steps(ds_file):
            f = field_2d(ds, args.var, args.height)
            vmax = f if vmax is None else np.fmax(vmax, f)
            if args.contour > -100:
                zh_list.append(field_2d(ds, 'Zh', args.height))
            times.append(time_string(ds))
            if coords is None:
                desc = describe(ds) or os.path.basename(path)
                if 'lon' in ds.coords and ds.lon.ndim == 2:
                    coords = (ds.lon.values, ds.lat.values, True)
                else:
                    x, y = horizontal_km(ds)
                    coords = (x, y, False)
        ds_file.close()
        print(f'  {os.path.basename(path)}')

    X, Y, geo = coords
    lo, hi, step = args.bounds or var_info(args.var)['bounds']
    norm = BoundaryNorm(np.arange(lo, hi + step / 2, step), RADAR_CMAP.N, extend='both')
    ax = None
    if geo:
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            fig = plt.figure(figsize=(10, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
            gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.top_labels = gl.right_labels = False
            kw = dict(transform=ccrs.PlateCarree())
        except ImportError:
            ax = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
        kw = {}
        ax.set_xlabel('longitude' if geo else 'x (km)')
        ax.set_ylabel('latitude' if geo else 'y (km)')
        ax.set_aspect(1. / np.cos(np.deg2rad(np.nanmean(Y))) if geo else 'equal')
    im = ax.pcolormesh(X, Y, vmax, cmap=RADAR_CMAP, norm=norm, shading='auto', **kw)
    for zh in zh_list:
        if np.nanmax(zh) >= args.contour:
            ax.contour(X, Y, np.nan_to_num(zh, nan=-999.), levels=[args.contour], colors='k', linewidths=0.7, **kw)
    plt.colorbar(im, ax=ax, label=axis_label(args.var), shrink=0.7, extend='both')
    where = f'at {args.height:g} km' if args.height is not None else 'column max'
    period = f'{times[0]} to {times[-1]}' if len(times) > 1 else times[0]
    cont = f', contours Zh = {args.contour:g} dBZ of each time' if args.contour > -100 else ''
    ax.set_title(f'{desc}: maximum of {args.var} ({where}) over {len(times)} time(s)\n{period}{cont}', fontsize=10)
    os.makedirs(args.out, exist_ok=True)
    tag = f'{args.height:g}km' if args.height is not None else 'colmax'
    stem = os.path.splitext(os.path.basename(files[0]))[0]
    path = os.path.join(args.out, f'{stem}_{len(times)}times_{args.var}_{tag}.png')
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


if __name__ == '__main__':
    main()
