#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick look and sanity checks of one operadar output file (any model, real or idealised case; no cartopy needed).

Prints (and writes to <out>/<file>_check.txt) a summary of each variable (valid points, NaN fraction,
percentiles) and warnings on suspicious values, then draws:
  <file>_maps.png      column maximum of each radar variable (x/y in km)
  <file>_contents.png  column maximum of each hydrometeor content
  <file>_joint.png     joint distributions Zh-Zdr, Zh-Kdp (Zh-Rhohv) in rain (T >= --t-rain) and in ice
                       (T <= --t-ice), with the median of the y variable for each Zh bin
  <file>_temperature.png  median and Q5-Q95 of each radar variable as a function of temperature (melting
                       layer, temperature dependence of the tables)

Checks:
  - very low Zh (below the default sensitivity of the band: -30 dBZ L-X, -40 Ku-Ka, -50 W): radar
    sensitivity not applied (operadar option sensitivity_threshold)
  - points with exactly Zh = 0 dBZ (false echoes of older operadar versions)
  - points with precipitation (contents without cloud water >= 0.1 g/m3) but no radar value: out of the tables?
  - Zdr < -1 dB or > 8 dB, Kdp < -1 deg/km (not at Ka and W band, where it is expected), Rhohv > 1 or < 0.7

Examples (from plot_tools/):
  python plot_output_check.py ../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
  python plot_output_check.py "../operadarFiles/*.nc" --figs maps joint --out IMG/check/
"""
import argparse
import os

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm

from utils.operadar_output import (HYDRO_NAMES, RADAR_CMAP, RADAR_VARS, axis_label, column_max, describe,
                                   expand_files, horizontal_km, open_output, time_string, var_info,
                                   vertical_dim, wet_fraction)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', help='operadar output file(s)')
    p.add_argument('--figs', nargs='+', default=['maps', 'contents', 'joint', 'temperature'],
                   choices=['maps', 'contents', 'joint', 'temperature', 'none'])
    p.add_argument('--zh-min', type=float, default=0., help='joint and temperature plots: Zh >= zh-min (dBZ)')
    p.add_argument('--t-rain', type=float, default=4., help='rain: T >= t-rain (degC)')
    p.add_argument('--t-ice', type=float, default=-10., help='ice: T <= t-ice (degC)')
    p.add_argument('--max-points', type=float, default=5e6,
                   help='random subsample for the percentiles and joint/temperature plots (default 5e6)')
    p.add_argument('--box', nargs=4, type=float, default=None, metavar=('LONMIN', 'LONMAX', 'LATMIN', 'LATMAX'))
    p.add_argument('--out', default='IMG/check/', help='output folder (default IMG/check/)')
    p.add_argument('--dpi', type=int, default=120)
    return p.parse_args()


# ====== Summary and checks ======
def subsample(n, max_points, seed=0):
    """Indices of a random subsample of n points (all points if n <= max_points)."""
    if n <= max_points:
        return slice(None)
    return np.sort(np.random.default_rng(seed).choice(n, int(max_points), replace=False))


# Default minimum detectable reflectivity at 1 km per band, as in operadar/radar/sensitivity.py
DEFAULT_ZMIN_1KM = {"L": -30.0, "S": -30.0, "C": -30.0, "X": -30.0,
                    "Ku": -40.0, "K": -40.0, "Ka": -40.0, "W": -50.0}


def summary(ds, radar_vars, log, max_points):
    vdim = vertical_dim(ds)
    log(f"{describe(ds)}  {time_string(ds)}  grid {dict(ds.sizes)}")
    log("attributes: " + ', '.join(f'{k}={v}' for k, v in ds.attrs.items() if not k.startswith('_')))
    if 'T' in ds:
        log(f"T: {float(ds.T.min()):.1f} to {float(ds.T.max()):.1f} degC")
    log(f"\n{'variable':10s} {'valid':>10s} {'NaN %':>7s} {'min':>9s} {'P1':>9s} {'P50':>9s} {'P99':>9s} {'max':>9s}")
    values = {}
    for v in radar_vars:
        a = ds[v].values.ravel()
        ok = np.isfinite(a)
        values[v] = a
        if ok.any():
            valid = a[ok]
            p = np.percentile(valid[subsample(valid.size, max_points)], [1, 50, 99])
            log(f"{v:10s} {ok.sum():10d} {100 * (1 - ok.mean()):7.1f} {a[ok].min():9.3f} {p[0]:9.3f} {p[1]:9.3f} "
                f"{p[2]:9.3f} {a[ok].max():9.3f}")
        else:
            log(f"{v:10s} {0:10d} {100.:7.1f}")
    total = None
    if 'Contents' in ds:
        log(f"\n{'content':14s} {'max (g/m3)':>11s} {'% > 0':>7s}")
        for h in ds.hydrometeor.values:
            m = ds.Contents.sel(hydrometeor=h).values
            if str(h) not in ('cc', 'cl', 'cs'):  # precipitation only (cloud water alone is below the radar sensitivity)
                total = m.copy() if total is None else total + m
            log(f"{HYDRO_NAMES.get(str(h), str(h)):14s} {np.nanmax(m):11.3f} {100 * np.mean(m > 0):7.2f}")
        fw, why = wet_fraction(ds)
        if fw is not None:
            f = fw.values[np.isfinite(fw.values)]
            log(f"Fw of the wet species: {f.size} points, median {np.median(f) if f.size else np.nan:.2f}")
        else:
            log(f"Fw: {why}")

    log('\nChecks:')
    n_warn = 0

    def warn(msg):
        nonlocal n_warn
        n_warn += 1
        log('  /!\\ ' + msg)

    band = str(ds.attrs.get('radar_band', ''))
    if 'Zh' in values:
        zh = values['Zh']
        zmin = DEFAULT_ZMIN_1KM.get(band, -30.)
        low = np.sum(zh < zmin - 0.01)
        if low:
            warn(f'{low} points with Zh < {zmin:.0f} dBZ (min {np.nanmin(zh):.0f}): radar sensitivity not applied '
                 '(operadar option sensitivity_threshold, or older file), Zmin_1km set lower in the '
                 'configuration, or ranges below 1 km')
        zero = int(np.sum(zh == 0.))
        if zero:
            warn(f'{zero} points with exactly Zh = 0 dBZ: false echoes of operadar versions before Nov. 2026 '
                 '(linear Z = 0 written as 0 dBZ), recompute the file')
        if total is not None:
            precip = total.ravel() >= 0.1
            nozh = precip & ~np.isfinite(zh)
            if 'Ah' in values and nozh.sum():
                # Ah is not masked by the sensitivity threshold: where it exists, the point was inside the tables
                masked = nozh & np.isfinite(values['Ah'])
                if masked.sum():
                    log(f'  (info) {masked.sum()} points ({100 * masked.sum() / max(precip.sum(), 1):.2f} % of the '
                        'points with contents >= 0.1 g/m3) have Ah but no Zh: removed by the radar sensitivity '
                        '(Zh_att below the threshold after attenuation)')
                nozh = nozh & ~masked
            if nozh.sum():
                warn(f'{nozh.sum()} points ({100 * nozh.sum() / max(precip.sum(), 1):.2f} % of the points with '
                     'contents >= 0.1 g/m3) have no Zh: outside the distance range, or contents/temperatures '
                     'outside the tables?')
    for v, cond, text in [('Zdr', lambda a: a < -1, 'Zdr < -1 dB'), ('Zdr', lambda a: a > 8, 'Zdr > 8 dB'),
                          ('Kdp', lambda a: a < -1, 'Kdp < -1 deg/km'), ('Rhohv', lambda a: a > 1.0001, 'Rhohv > 1'),
                          ('Rhohv', lambda a: a < 0.7, 'Rhohv < 0.7')]:
        if v in values:
            with np.errstate(invalid='ignore'):
                n = int(np.sum(cond(values[v])))
            if n:
                if text.startswith('Kdp <') and band in ('Ka', 'W'):
                    log(f'  (info) {n} points with {text}: expected at {band} band (Mie resonance of the '
                        'raindrops: Kdp of rain and melting particles becomes negative)')
                else:
                    warn(f'{n} points with {text}')
    if n_warn == 0:
        log('  nothing suspicious')
    return vdim


# ====== Figures ======
def binned_percentiles(x, y, edges, percents, nmin=10):
    """Percentiles of y in each bin of x (NaN for bins with fewer than nmin values)."""
    order = np.argsort(x)
    x, y = x[order], y[order]
    cuts = np.searchsorted(x, edges)
    out = np.full((len(edges) - 1, len(percents)), np.nan)
    for k in range(len(edges) - 1):
        vals = y[cuts[k]:cuts[k + 1]]
        if vals.size >= nmin:
            out[k] = np.percentile(vals, percents)
    return out


def grid_axes(n, size=4.2):
    ncol = min(n, 3)
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(size * ncol + 1, size * nrow), squeeze=False)
    for ax in axs.ravel()[n:]:
        ax.set_visible(False)
    return fig, axs.ravel()


def map_axes(ds):
    if 'lon' in ds.coords and 'lat' in ds.coords and ds.lon.ndim == 2:
        return ds.lon.values, ds.lat.values, 'longitude', 'latitude'
    x, y = horizontal_km(ds)
    return x, y, 'x (km)', 'y (km)'


def fig_maps(ds, radar_vars, vdim, title, path, dpi):
    X, Y, xl, yl = map_axes(ds)
    fig, axs = grid_axes(len(radar_vars))
    for ax, v in zip(axs, radar_vars):
        lo, hi, step = var_info(v)['bounds']
        bounds = np.arange(lo, hi + step / 2, step)
        norm = BoundaryNorm(bounds, RADAR_CMAP.N, extend='both')
        im = ax.pcolormesh(X, Y, column_max(ds[v], vdim), cmap=RADAR_CMAP, norm=norm, shading='auto')
        fig.colorbar(im, ax=ax, label=axis_label(v), shrink=0.85, extend='both')
        ax.set_title(f'{v} column max')
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_aspect('equal' if xl.startswith('x') else 1. / np.cos(np.deg2rad(np.nanmean(Y))))
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


def fig_contents(ds, vdim, title, path, dpi):
    X, Y, xl, yl = map_axes(ds)
    hyd = [str(h) for h in ds.hydrometeor.values]
    fig, axs = grid_axes(len(hyd), size=3.8)
    norm = LogNorm(vmin=0.01, vmax=10.)
    cmap = mpl.colormaps['viridis'].copy()
    cmap.set_under('white')
    for ax, h in zip(axs, hyd):
        m = column_max(ds.Contents.sel(hydrometeor=h), vdim)
        im = ax.pcolormesh(X, Y, np.where(m > 0, m, np.nan), cmap=cmap, norm=norm, shading='auto')
        fig.colorbar(im, ax=ax, label=r'g m$^{-3}$', shrink=0.85, extend='both')
        ax.set_title(f'{HYDRO_NAMES.get(h, h)} ({h}) column max, max {np.nanmax(m):.2f}', fontsize=10)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_aspect('equal' if xl.startswith('x') else 1. / np.cos(np.deg2rad(np.nanmean(Y))))
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


def fig_joint(ds, radar_vars, args, title, path, dpi):
    ys = [v for v in ['Zdr', 'Kdp', 'Rhohv'] if v in radar_vars]
    if 'Zh' not in ds or not ys or 'T' not in ds:
        return
    sub = subsample(ds.Zh.size, args.max_points)
    zh, T = ds.Zh.values.ravel()[sub], ds['T'].values.ravel()[sub]
    regimes = [(f'rain (T ≥ {args.t_rain:g} °C)', T >= args.t_rain), (f'ice (T ≤ {args.t_ice:g} °C)', T <= args.t_ice)]
    fig, axs = plt.subplots(len(regimes), len(ys), figsize=(4.3 * len(ys) + 1, 3.8 * len(regimes)), squeeze=False)
    xb = np.arange(args.zh_min, 72., 1.)
    for i, (name, sel) in enumerate(regimes):
        for j, v in enumerate(ys):
            ax = axs[i, j]
            y = ds[v].values.ravel()[sub]
            ok = sel & np.isfinite(zh) & np.isfinite(y) & (zh >= args.zh_min)
            lim = {'Zdr': (-1., 6.), 'Kdp': (-0.5, 6.), 'Rhohv': (0.9, 1.001)}[v]
            yb = np.linspace(lim[0], lim[1], 120)
            if ok.sum() == 0:
                ax.text(0.5, 0.5, 'no data', ha='center', transform=ax.transAxes)
            else:
                h, _, _ = np.histogram2d(zh[ok], np.clip(y[ok], lim[0], lim[1] - 1e-9), bins=[xb, yb])
                im = ax.pcolormesh(xb, yb, np.where(h > 0, h, np.nan).T, norm=LogNorm(), cmap='viridis')
                fig.colorbar(im, ax=ax, label='number of points', shrink=0.9)
                med = binned_percentiles(zh[ok], y[ok], xb, [50], nmin=10)[:, 0]
                ax.plot(0.5 * (xb[1:] + xb[:-1]), med, color='r', lw=1.5, label='median')
                ax.legend(loc='upper left', fontsize=8)
            ax.set_xlabel(axis_label('Zh'))
            ax.set_ylabel(axis_label(v))
            ax.set_ylim(*lim)
            ax.set_title(f'{name}: {ok.sum()} points' + (' (subsample)' if isinstance(sub, np.ndarray) else ''),
                         fontsize=10)
            ax.grid(ls='--', lw=0.4, alpha=0.6)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


def fig_temperature(ds, radar_vars, args, title, path, dpi):
    if 'T' not in ds or 'Zh' not in ds:
        return
    sub = subsample(ds.Zh.size, args.max_points)
    T, zh = ds['T'].values.ravel()[sub], ds.Zh.values.ravel()[sub]
    tb = np.arange(np.floor(np.nanmin(T) / 2) * 2, np.nanmax(T) + 2, 2.)
    fig, axs = plt.subplots(1, len(radar_vars), figsize=(3.4 * len(radar_vars), 5), squeeze=False, sharey=True)
    for ax, v in zip(axs[0], radar_vars):
        a = ds[v].values.ravel()[sub]
        ok = np.isfinite(a) & np.isfinite(zh) & (zh >= args.zh_min)
        q = binned_percentiles(T[ok], a[ok], tb, [5, 50, 95], nmin=20)
        tc = 0.5 * (tb[1:] + tb[:-1])
        ax.fill_betweenx(tc, q[:, 0], q[:, 2], color='tab:blue', alpha=0.2, lw=0, label='Q5-Q95')
        ax.plot(q[:, 1], tc, color='tab:blue', lw=2, label='median')
        ax.axhline(0, color='k', lw=0.8, ls='--')
        ax.set_xlabel(axis_label(v))
        ax.grid(ls='--', lw=0.4, alpha=0.6)
    axs[0, 0].set_ylabel('Temperature (°C)')
    axs[0, 0].invert_yaxis()
    axs[0, 0].legend(fontsize=8, loc='lower right')
    fig.suptitle(f'{title}\npoints with Zh ≥ {args.zh_min:g} dBZ, 2 °C bins')
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    for path in expand_files(args.files):
        ds = open_output(path, args.box)
        base = os.path.splitext(os.path.basename(path))[0]
        lines = []

        def log(msg=''):
            print(msg)
            lines.append(msg)

        log(f'=== {path}')
        radar_vars = [v for v in RADAR_VARS if v in ds]
        vdim = summary(ds, radar_vars, log, args.max_points)
        with open(os.path.join(args.out, f'{base}_check.txt'), 'w') as f:
            f.write('\n'.join(lines) + '\n')
        title = f'{base}\n{describe(ds)}  {time_string(ds)}'
        figs = set(args.figs)
        if 'maps' in figs and radar_vars:
            fig_maps(ds, radar_vars, vdim, title, os.path.join(args.out, f'{base}_maps.png'), args.dpi)
        if 'contents' in figs and 'Contents' in ds:
            fig_contents(ds, vdim, title, os.path.join(args.out, f'{base}_contents.png'), args.dpi)
        if 'joint' in figs:
            fig_joint(ds, radar_vars, args, title, os.path.join(args.out, f'{base}_joint.png'), args.dpi)
        if 'temperature' in figs:
            fig_temperature(ds, [v for v in radar_vars if v != 'Zh_att'], args, title,
                            os.path.join(args.out, f'{base}_temperature.png'), args.dpi)
        ds.close()


if __name__ == '__main__':
    main()
