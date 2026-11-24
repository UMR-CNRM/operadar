#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare two operadar outputs computed on the same model file (e.g. old and new lookup tables, two
configurations, two radar bands, T-matrix and Rayleigh).

Prints the statistics of B - A for each radar variable (number of points, mean difference, RMSE, correlation),
for all points and by temperature range (rain T >= 4 degC, melting layer 0-4 degC, ice T < 0 degC), and draws:
  <name>_scatter.png      B versus A (2D histogram, log colour scale) with the 1:1 line
  <name>_maps.png         column maximum of A, B and B - A for each variable
  <name>_temperature.png  median and Q5-Q95 of B - A as a function of temperature

Points are compared where both values exist and max(Zh_A, Zh_B) >= --zh-min.
If one file is a horizontal subdomain of the other (operadar subDomain), the larger one is cut to it
(located by lat/lon); --box LONMIN LONMAX LATMIN LATMAX compares only a lon/lat box.

Examples (from plot_tools/):
  python plot_compare_outputs.py old/dpolvar_Arome_ICE3_Cband_20250831_1800.nc \\
                                 new/dpolvar_Arome_ICE3_Cband_20250831_1800.nc --labels "old tables" "new tables"
  python plot_compare_outputs.py FILE_C FILE_X --labels C X --vars Zh Zdr --zh-min 10
"""
import argparse
import os
import re

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm, TwoSlopeNorm

from utils.operadar_output import (RADAR_CMAP, RADAR_VARS, axis_label, column_max, describe, horizontal_km,
                                   open_output, var_info, vertical_dim)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file_a', help='reference file A')
    p.add_argument('file_b', help='file B')
    p.add_argument('--labels', nargs=2, default=['A', 'B'])
    p.add_argument('--vars', nargs='+', default=None, help='radar variables (default: those in both files)')
    p.add_argument('--zh-min', type=float, default=0., help='compare where max(Zh_A, Zh_B) >= zh-min (dBZ)')
    p.add_argument('--time', type=int, default=0, help='time index (files with a time dimension)')
    p.add_argument('--box', nargs=4, type=float, default=None, metavar=('LONMIN', 'LONMAX', 'LATMIN', 'LATMAX'),
                   help='compare only this lon/lat box')
    p.add_argument('--out', default='IMG/compare/', help='output folder (default IMG/compare/)')
    p.add_argument('--name', default=None, help='prefix of the figures (default: from the labels)')
    p.add_argument('--dpi', type=int, default=120)
    return p.parse_args()


def stats(a, b):
    d = b - a
    if d.size < 2:
        return dict(n=d.size, bias=np.nan, rmse=np.nan, r=np.nan, p99=np.nan)
    r = np.corrcoef(a, b)[0, 1] if np.std(a) > 0 and np.std(b) > 0 else np.nan
    return dict(n=d.size, bias=d.mean(), rmse=np.sqrt(np.mean(d ** 2)), r=r, p99=np.percentile(np.abs(d), 99))


def align_subdomain(A, B):
    """If one file is a horizontal subdomain of the other (e.g. an operadar run with subDomain and a run on the
    whole model domain), cut the larger one to the grid of the smaller one, located by its lat/lon."""
    if (A.sizes.get('y'), A.sizes.get('x')) == (B.sizes.get('y'), B.sizes.get('x')):
        return A, B
    small, large = (A, B) if A.sizes['y'] * A.sizes['x'] <= B.sizes['y'] * B.sizes['x'] else (B, A)
    if 'lat' not in small or small.lat.ndim != 2 or 'lat' not in large or large.lat.ndim != 2:
        return A, B
    dist = np.abs(large.lat.values - float(small.lat[0, 0])) + np.abs(large.lon.values - float(small.lon[0, 0]))
    j, i = np.unravel_index(np.argmin(dist), dist.shape)
    ny, nx = small.sizes['y'], small.sizes['x']
    if dist[j, i] > 1e-4 or j + ny > large.sizes['y'] or i + nx > large.sizes['x']:
        return A, B
    cut = large.isel(y=slice(j, j + ny), x=slice(i, i + nx))
    if not (np.allclose(cut.lat.values, small.lat.values, atol=1e-4)
            and np.allclose(cut.lon.values, small.lon.values, atol=1e-4)):
        return A, B
    print(f'One file is a subdomain of the other: compared on its {ny} x {nx} points (offset y={j}, x={i})')
    return (cut, B) if large is A else (A, cut)


def main():
    args = parse_args()
    la, lb = args.labels
    A, B = open_output(args.file_a, args.box), open_output(args.file_b, args.box)
    A, B = [d.isel(time=args.time) if 'time' in d.dims else d for d in (A, B)]
    A, B = align_subdomain(A, B)
    for dim in ('x', 'y', 'level', 'z'):
        if A.sizes.get(dim) != B.sizes.get(dim):
            raise SystemExit(f'the files are not on the same grid ({dim}: {A.sizes.get(dim)} / {B.sizes.get(dim)})')
    if 'Contents' in A and 'Contents' in B:
        dm = float(np.nanmax(np.abs(A.Contents.values - B.Contents.values)))
        if dm > 1e-6:
            print(f'/!\\ the contents differ (max difference {dm:.3g} g/m3): not the same model file?')
    radar_vars = args.vars or [v for v in RADAR_VARS if v in A and v in B]
    vdim = vertical_dim(A)
    name = args.name or 'compare_' + '_vs_'.join(re.sub(r'\W+', '', lab) for lab in (lb, la))
    os.makedirs(args.out, exist_ok=True)
    print(f'A = {la}: {args.file_a} ({describe(A)})\nB = {lb}: {args.file_b} ({describe(B)})')

    zha, zhb = A.Zh.values.ravel(), B.Zh.values.ravel()
    with np.errstate(invalid='ignore'):
        strong = np.fmax(zha, zhb) >= args.zh_min
    T = A['T'].values.ravel() if 'T' in A else None
    regimes = [('all', np.ones_like(strong))]
    if T is not None:
        regimes += [('T >= 4', T >= 4), ('0 <= T < 4', (T >= 0) & (T < 4)), ('T < 0', T < 0)]

    lines = [f'B - A, points where both exist and max(Zh_A, Zh_B) >= {args.zh_min:g} dBZ',
             f"{'variable':8s} {'points':12s} {'n':>10s} {'mean':>9s} {'RMSE':>9s} {'|d| P99':>9s} {'corr':>7s}"]
    data = {}
    for v in radar_vars:
        a, b = A[v].values.ravel(), B[v].values.ravel()
        ok = np.isfinite(a) & np.isfinite(b) & strong
        only_a, only_b = int(np.sum(np.isfinite(a) & ~np.isfinite(b))), int(np.sum(~np.isfinite(a) & np.isfinite(b)))
        data[v] = (a, b, ok)
        for rname, sel in regimes:
            s = stats(a[ok & sel], b[ok & sel])
            lines.append(f"{v:8s} {rname:12s} {s['n']:10d} {s['bias']:9.3f} {s['rmse']:9.3f} {s['p99']:9.3f} "
                         f"{s['r']:7.4f}")
        if only_a or only_b:
            lines.append(f"{'':8s} values only in A: {only_a}, only in B: {only_b}")
    print('\n'.join(lines))
    with open(os.path.join(args.out, f'{name}_stats.txt'), 'w') as f:
        f.write(f'A = {la}: {args.file_a}\nB = {lb}: {args.file_b}\n' + '\n'.join(lines) + '\n')

    # Scatter
    n = len(radar_vars)
    fig, axs = plt.subplots(1, n, figsize=(4.4 * n, 4.2), squeeze=False)
    for ax, v in zip(axs[0], radar_vars):
        a, b, ok = data[v]
        if ok.sum() == 0:
            continue
        lo, hi = np.percentile(np.concatenate([a[ok], b[ok]]), [0.1, 99.9])
        lo, hi = lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo)
        bins = np.linspace(lo, hi, 150)
        h, _, _ = np.histogram2d(a[ok], b[ok], bins=[bins, bins])
        im = ax.pcolormesh(bins, bins, np.where(h > 0, h, np.nan).T, norm=LogNorm(), cmap='viridis')
        ax.plot([lo, hi], [lo, hi], 'k--', lw=0.8)
        s = stats(a[ok], b[ok])
        ax.text(0.03, 0.97, f"n = {s['n']}\nmean B-A = {s['bias']:.3f}\nRMSE = {s['rmse']:.3f}\nr = {s['r']:.4f}",
                transform=ax.transAxes, va='top', fontsize=8, bbox=dict(fc='w', alpha=0.8, lw=0))
        ax.set_xlabel(f'{la}: {axis_label(v)}')
        ax.set_ylabel(f'{lb}: {axis_label(v)}')
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, label='points', shrink=0.8)
    fig.suptitle(f'{lb} vs {la}  ({describe(A)})')
    fig.tight_layout()
    path = os.path.join(args.out, f'{name}_scatter.png')
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)

    # Maps of the column maximum
    x, y = horizontal_km(A)
    fig, axs = plt.subplots(n, 3, figsize=(13, 3.9 * n), squeeze=False, layout='constrained')
    for k, v in enumerate(radar_vars):
        ma, mb = column_max(A[v], vdim), column_max(B[v], vdim)
        lo, hi, step = var_info(v)['bounds']
        norm = BoundaryNorm(np.arange(lo, hi + step / 2, step), RADAR_CMAP.N, extend='both')
        for j, (m, lab) in enumerate([(ma, la), (mb, lb)]):
            im = axs[k, j].pcolormesh(x, y, m, cmap=RADAR_CMAP, norm=norm, shading='auto')
            axs[k, j].set_title(f'{lab}: {v} column max', fontsize=10)
        fig.colorbar(im, ax=axs[k, :2], label=axis_label(v), pad=0.01, extend='both')
        d = mb - ma
        amax = np.nanpercentile(np.abs(d), 99.5) if np.isfinite(d).any() else 1.
        amax = amax if amax > 0 else 1e-3
        imd = axs[k, 2].pcolormesh(x, y, d, cmap='RdBu_r', norm=TwoSlopeNorm(0., -amax, amax), shading='auto')
        axs[k, 2].set_title(f'{lb} - {la}', fontsize=10)
        fig.colorbar(imd, ax=axs[k, 2], pad=0.01, extend='both')
        for ax in axs[k]:
            ax.set_aspect('equal')
            ax.set_xlabel('x (km)', fontsize=8)
            ax.set_ylabel('y (km)', fontsize=8)
            ax.tick_params(labelsize=7)
    fig.suptitle(f'{describe(A)}: column maximum')
    path = os.path.join(args.out, f'{name}_maps.png')
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)

    # Difference as a function of temperature
    if T is not None:
        tb = np.arange(np.floor(np.nanmin(T) / 2) * 2, np.nanmax(T) + 2, 2.)
        tc = 0.5 * (tb[1:] + tb[:-1])
        fig, axs = plt.subplots(1, n, figsize=(3.4 * n, 5), squeeze=False, sharey=True)
        for ax, v in zip(axs[0], radar_vars):
            a, b, ok = data[v]
            order = np.argsort(T[ok])
            ts, ds_ = T[ok][order], (b - a)[ok][order]
            cuts = np.searchsorted(ts, tb)
            q = np.full((len(tc), 3), np.nan)
            for k in range(len(tc)):
                vals = ds_[cuts[k]:cuts[k + 1]]
                if vals.size >= 20:
                    q[k] = np.percentile(vals, [5, 50, 95])
            ax.fill_betweenx(tc, q[:, 0], q[:, 2], color='tab:red', alpha=0.2, lw=0, label='Q5-Q95')
            ax.plot(q[:, 1], tc, color='tab:red', lw=2, label='median')
            ax.axvline(0, color='k', lw=0.8)
            ax.axhline(0, color='k', lw=0.6, ls='--')
            ax.set_xlabel(f'{lb} - {la}: {axis_label(v)}', fontsize=9)
            ax.grid(ls='--', lw=0.4, alpha=0.6)
        axs[0, 0].set_ylabel('Temperature (°C)')
        axs[0, 0].invert_yaxis()
        axs[0, 0].legend(fontsize=8, loc='lower right')
        fig.suptitle(f'{lb} - {la} as a function of temperature (2 °C bins)')
        fig.tight_layout()
        path = os.path.join(args.out, f'{name}_temperature.png')
        fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
        plt.close(fig)
        print('  ->', path)


if __name__ == '__main__':
    main()
