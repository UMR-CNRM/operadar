#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vertical cross-sections of an operadar output (or of a gridded file): radar variables with the isotherms, and
optionally the hydrometeor contents. A small map of the column maximum of Zh shows where the section is.

Where is the section (default: west-east line through the maximum of Zh):
  --x KM                      south-north line at x = KM (km from the west edge of the domain)
  --y KM                      west-east line at y = KM (km from the south edge)
  --points LON1 LAT1 LON2 LAT2  between two points (needs lat/lon, real_case=True)
  --ij I1 J1 I2 J2            between two grid points (x index, y index)

Heights are altitudes above sea level (the bottom of the section follows the model orography).

Examples (from plot_tools/):
  python plot_cross_section.py ../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
  python plot_cross_section.py FILE --points 2.3 45.2 2.8 45.9 --vars Zh Zdr Kdp Rhohv --contents
  python plot_cross_section.py FILE --y 60 --hmax 12 --isotherms 0 -10 -20 -40
"""
import argparse
import os

import numpy as np
import xarray as xr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm

from utils.operadar_output import (HYDRO_NAMES, RADAR_CMAP, axis_label, column_max, describe, horizontal_km,
                                   open_output, time_string, var_info, vertical_dim)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='operadar output file (or gridded file)')
    where = p.add_mutually_exclusive_group()
    where.add_argument('--x', type=float, help='south-north section at x = X km')
    where.add_argument('--y', type=float, help='west-east section at y = Y km')
    where.add_argument('--points', nargs=4, type=float, metavar=('LON1', 'LAT1', 'LON2', 'LAT2'))
    where.add_argument('--ij', nargs=4, type=int, metavar=('I1', 'J1', 'I2', 'J2'))
    p.add_argument('--vars', nargs='+', default=None, help='radar variables (default: Zh Zdr Kdp Rhohv found)')
    p.add_argument('--contents', action='store_true', help='also plot the hydrometeor contents')
    p.add_argument('--hmax', type=float, default=15., help='top of the section (km)')
    p.add_argument('--isotherms', nargs='+', type=float, default=[0., -20., -40.], help='isotherms (degC)')
    p.add_argument('--time', type=int, default=0, help='time index (files with a time dimension)')
    p.add_argument('--out', default='IMG/sections/', help='output folder (default IMG/sections/)')
    p.add_argument('--dpi', type=int, default=130)
    return p.parse_args()


def nearest_index(ds, lon, lat):
    d2 = (ds.lon.values - lon) ** 2 * np.cos(np.deg2rad(lat)) ** 2 + (ds.lat.values - lat) ** 2
    j, i = np.unravel_index(np.argmin(d2), d2.shape)
    return int(i), int(j)


def section_indices(args, ds):
    """Grid indices (i along x, j along y) of the points of the section, and a description."""
    nx, ny = ds.sizes['x'], ds.sizes['y']
    x, y = horizontal_km(ds)
    if args.x is not None:
        i = int(np.argmin(np.abs(x - args.x)))
        return (i, 0), (i, ny - 1), f'x = {x[i]:.0f} km'
    if args.y is not None:
        j = int(np.argmin(np.abs(y - args.y)))
        return (0, j), (nx - 1, j), f'y = {y[j]:.0f} km'
    if args.points is not None:
        if 'lon' not in ds.coords:
            raise SystemExit('--points needs lat/lon in the file (real_case=True); use --ij, --x or --y')
        lon1, lat1, lon2, lat2 = args.points
        return (nearest_index(ds, lon1, lat1), nearest_index(ds, lon2, lat2),
                f'({lon1:.2f}E {lat1:.2f}N) to ({lon2:.2f}E {lat2:.2f}N)')
    if args.ij is not None:
        i1, j1, i2, j2 = args.ij
        return (i1, j1), (i2, j2), f'grid points ({i1},{j1}) to ({i2},{j2})'
    zmax = np.nan_to_num(column_max(ds.Zh, vertical_dim(ds)), nan=-999.)
    j, i = np.unravel_index(np.argmax(zmax), zmax.shape)
    return (0, int(j)), (nx - 1, int(j)), f'west-east through the Zh maximum (y = {y[j]:.0f} km)'


def sample(ds, p1, p2):
    """Indices of the grid points along the segment p1-p2 (nearest neighbour) and the distance (km)."""
    n = int(np.ceil(np.hypot(p2[0] - p1[0], p2[1] - p1[1]))) + 1
    ii = np.rint(np.linspace(p1[0], p2[0], n)).astype(int)
    jj = np.rint(np.linspace(p1[1], p2[1], n)).astype(int)
    res = float(ds.attrs.get('horizontal_resolution', 1000.)) / 1000.
    dist = np.hypot(ii - ii[0], jj - jj[0]) * res
    return ii, jj, dist


def vertical_slice(ds, da, ii, jj):
    s = da.isel(x=xr.DataArray(ii, dims='s'), y=xr.DataArray(jj, dims='s'))
    vdim = vertical_dim(ds)
    return s.transpose(vdim, 's').values


def heights_2d(ds, ii, jj):
    if 'Alt' in ds:
        return vertical_slice(ds, ds.Alt, ii, jj) / 1000.
    return np.repeat(ds.z.values[:, None] / 1000., len(ii), axis=1)


def draw_isotherms(ax, dist2, alt2, T, levels):
    if T is None:
        return
    cs = ax.contour(dist2, alt2, T, levels=sorted(levels), colors='k', linewidths=0.8, linestyles='--')
    ax.clabel(cs, fmt='%g°C', fontsize=7)


def locator_map(ax, ds, ii, jj):
    x, y = horizontal_km(ds)
    zmax = column_max(ds.Zh, vertical_dim(ds))
    lo, hi, step = var_info('Zh')['bounds']
    norm = BoundaryNorm(np.arange(lo, hi + step / 2, step), RADAR_CMAP.N, extend='both')
    ax.pcolormesh(x, y, zmax, cmap=RADAR_CMAP, norm=norm, shading='auto')
    ax.plot(x[ii], y[jj], 'k-', lw=1.5)
    ax.plot(x[ii[0]], y[jj[0]], 'ko', ms=4)
    ax.text(x[ii[0]], y[jj[0]], ' A', fontsize=9, fontweight='bold')
    ax.text(x[ii[-1]], y[jj[-1]], ' B', fontsize=9, fontweight='bold')
    ax.set_aspect('equal')
    ax.set_title('Zh column max', fontsize=9)
    ax.set_xlabel('x (km)', fontsize=8)
    ax.set_ylabel('y (km)', fontsize=8)
    ax.tick_params(labelsize=7)


def main():
    args = parse_args()
    ds = open_output(args.file)
    if 'time' in ds.dims:
        ds = ds.isel(time=args.time)
    radar_vars = args.vars or [v for v in ['Zh', 'Zdr', 'Kdp', 'Rhohv'] if v in ds]
    p1, p2, where = section_indices(args, ds)
    ii, jj, dist = sample(ds, p1, p2)
    alt2 = heights_2d(ds, ii, jj)
    dist2 = np.broadcast_to(dist, alt2.shape)
    T = vertical_slice(ds, ds['T'], ii, jj) if 'T' in ds else None
    base = os.path.splitext(os.path.basename(args.file))[0]
    title = f'{describe(ds)}  {time_string(ds)}  section {where}'
    os.makedirs(args.out, exist_ok=True)
    tag = '_'.join(map(str, (*p1, *p2)))

    # Radar variables
    n = len(radar_vars)
    fig = plt.figure(figsize=(13, 2.6 * n + 0.8))
    gs = fig.add_gridspec(n, 2, width_ratios=[1, 3.2], wspace=0.15)
    locator_map(fig.add_subplot(gs[0:min(2, n), 0]), ds, ii, jj)
    for k, v in enumerate(radar_vars):
        ax = fig.add_subplot(gs[k, 1])
        lo, hi, step = var_info(v)['bounds']
        norm = BoundaryNorm(np.arange(lo, hi + step / 2, step), RADAR_CMAP.N, extend='both')
        im = ax.pcolormesh(dist2, alt2, vertical_slice(ds, ds[v], ii, jj), cmap=RADAR_CMAP, norm=norm,
                           shading='nearest')
        draw_isotherms(ax, dist2, alt2, T, args.isotherms)
        fig.colorbar(im, ax=ax, label=axis_label(v), pad=0.01, extend='both')
        ax.set_ylim(max(0., np.nanmin(alt2) - 0.2), args.hmax)
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel('Distance from A (km)' if k == n - 1 else '')
        ax.set_title(v, fontsize=10, loc='left')
    fig.suptitle(title, fontsize=11)
    path = os.path.join(args.out, f'{base}_section_{tag}.png')
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)

    # Contents
    if args.contents and 'Contents' in ds:
        hyd = [str(h) for h in ds.hydrometeor.values]
        ncol = 2
        nrow = int(np.ceil(len(hyd) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(13, 2.6 * nrow + 0.8), squeeze=False, sharex=True, sharey=True)
        norm = LogNorm(vmin=0.01, vmax=10.)
        cmap = mpl.colormaps['viridis'].copy()
        cmap.set_under('white')
        for ax, h in zip(axs.ravel(), hyd):
            m = vertical_slice(ds, ds.Contents.sel(hydrometeor=h), ii, jj)
            im = ax.pcolormesh(dist2, alt2, np.where(m > 0, m, np.nan), cmap=cmap, norm=norm, shading='nearest')
            draw_isotherms(ax, dist2, alt2, T, args.isotherms)
            ax.set_title(f'{HYDRO_NAMES.get(h, h)} ({h}), max {np.nanmax(m):.2f} g/m3', fontsize=10, loc='left')
            ax.set_ylim(max(0., np.nanmin(alt2) - 0.2), args.hmax)
        for ax in axs.ravel()[len(hyd):]:
            ax.set_visible(False)
        for ax in axs[:, 0]:
            ax.set_ylabel('Altitude (km)')
        for ax in axs[-1]:
            ax.set_xlabel('Distance from A (km)')
        fig.colorbar(im, ax=axs, label=r'content (g m$^{-3}$)', pad=0.01, extend='both')
        fig.suptitle(title, fontsize=11)
        path = os.path.join(args.out, f'{base}_section_{tag}_contents.png')
        fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
        plt.close(fig)
        print('  ->', path)


if __name__ == '__main__':
    main()
