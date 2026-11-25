#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the CFADs computed by compute_cfad.py (one statistics file per dataset: observations, simulations...).

Figures (--figs, default: summary):
  summary   layout of David et al. (2025, AMT, Fig. 6): radar variables of all datasets (a, b, c...), contents of
            each simulation (one panel per dataset) and Fw of the wet graupel
  radar     radar variables only: median (line) and Q25-Q75 (shaded) of each dataset, one panel per variable
  contents  hydrometeor contents, one panel per dataset
  fw        liquid water fraction of the wet graupel/hail, median, Q25-Q75 and Q95 (dotted)
  freq      frequency CFADs (2D histograms, Yuter and Houze 1995), one figure per variable, one panel per
            dataset, and with --diff the difference with the first dataset (e.g. the observations)
  all       all of them

Statistics files with the same label are added (e.g. cases computed separately).

Examples (from plot_tools/):
  python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc
  python plot_cfad.py CFAD/*_core.nc --figs radar fw --vars Zh Zdr Kdp Rhohv --extremes
  python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc --figs freq --diff --vars Zh Zdr
  python plot_cfad.py CFAD/ICE3_core.nc --figs contents --contents-log --xlim M 0.001 10
"""
import argparse
import os
import re

import numpy as np
import xarray as xr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

from utils.operadar_output import HYDRO_COLORS, HYDRO_NAMES, axis_label, var_info
from utils import cfad

DATASET_COLORS = ['tab:orange', 'tab:green', 'tab:blue', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']
OBS_COLOR = 'dimgray'
HYDRO_ORDER = ['rr', 'cc', 'cl', 'cs', 'ii', 'ss', 'gg', 'wg', 'hh', 'wh']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stats', nargs='+', help='statistics files of compute_cfad.py')
    p.add_argument('--labels', nargs='+', default=None, help='dataset names (default: label of each file)')
    p.add_argument('--colors', nargs='+', default=None, help='one colour per dataset')
    p.add_argument('--figs', nargs='+', default=['summary'],
                   choices=['summary', 'radar', 'contents', 'fw', 'freq', 'all'])
    p.add_argument('--vars', nargs='+', default=None, help='radar variables (default: Zh Zdr Kdp found in the files)')
    p.add_argument('--hydrometeors', nargs='+', default=None, help='contents to plot (default: all)')
    p.add_argument('--extremes', action='store_true', help='add Q5 (dashed) and Q95 (dotted) to the radar CFADs')
    p.add_argument('--min-count', type=int, default=20, help='hide the heights with fewer values (default 20)')
    p.add_argument('--hmin', type=float, default=0., help='bottom of the plots (km)')
    p.add_argument('--hmax', type=float, default=12., help='top of the plots (km)')
    p.add_argument('--xlim', nargs=3, action='append', default=[], metavar=('VAR', 'MIN', 'MAX'),
                   help='x range of a variable (M for the contents), e.g. --xlim Zh 0 65 (repeatable)')
    p.add_argument('--contents-log', action='store_true', help='log x axis for the contents')
    p.add_argument('--freq-norm', default='level', choices=['level', 'total', 'count'],
                   help='frequency CFADs: %% of the values of each level (default), %% of all values, or counts')
    p.add_argument('--bin', nargs=2, action='append', default=[], metavar=('VAR', 'WIDTH'),
                   help='bin width of the frequency CFADs, e.g. --bin Zh 1 (repeatable)')
    p.add_argument('--log', action='store_true', help='logarithmic colour scale for the frequency CFADs')
    p.add_argument('--diff', action='store_true', help='frequency CFADs: add the difference with the first dataset')
    p.add_argument('--out', default='IMG/cfad/', help='output folder (default IMG/cfad/)')
    p.add_argument('--name', default=None, help='prefix of the figure names (default: from the labels and mask)')
    p.add_argument('--dpi', type=int, default=150)
    return p.parse_args()


# ====== Datasets ======
def load_datasets(args):
    """Ordered dict label -> merged statistics."""
    groups = {}
    for i, path in enumerate(args.stats):
        st = xr.load_dataset(path)
        label = args.labels[i] if args.labels and i < len(args.labels) else st.attrs.get('label', path)
        groups.setdefault(label, []).append(st)
    datasets = {lab: cfad.merge(lst) for lab, lst in groups.items()}
    for lab, lst in groups.items():
        if len(lst) > 1:
            print(f'{lab}: {len(lst)} statistics files added')
    return datasets


def dataset_colors(labels, user_colors):
    colors, k = {}, 0
    for i, lab in enumerate(labels):
        if user_colors and i < len(user_colors):
            colors[lab] = user_colors[i]
        elif re.search(r'obs', lab, re.IGNORECASE):
            colors[lab] = OBS_COLOR
        else:
            colors[lab] = DATASET_COLORS[k % len(DATASET_COLORS)]
            k += 1
    return colors


def mask_title(st, value=True) -> str:
    a = st.attrs
    m = a.get('mask', '')
    if m == 'core':
        if not value:
            return f"inside storm cores (column max {a.get('zh_var', 'Zh')} ≥"
        return f"inside storm cores (column max {a.get('zh_var', 'Zh')} ≥ {float(a.get('core_dbz', 40)):g} dBZ)"
    if m == 'envelope':
        return f"inside storm cells (column max {a.get('zh_var', 'Zh')} ≥ {float(a.get('envelope_dbz', 25)):g} dBZ)"
    if m == 'zh':
        return f"where {a.get('zh_var', 'Zh')} ≥ {float(a.get('zh_local', 40)):g} dBZ"
    if m == 'var':
        return f"inside {a.get('mask_var', 'mask')}"
    return 'all grid points'


def mask_text(datasets) -> str:
    """Description of the mask for the titles; different core thresholds between datasets (e.g. 48 dBZ for the
    observations, 40 dBZ for the simulations as in the paper) are listed."""
    st0 = next(iter(datasets.values()))
    cores = {lab: float(st.attrs.get('core_dbz')) for lab, st in datasets.items() if st.attrs.get('mask') == 'core'}
    masks = {st.attrs.get('mask') for st in datasets.values()}
    if len(masks) > 1:
        return 'for different masks (see the legend)'
    if len(set(cores.values())) > 1:
        return mask_title(st0, value=False) + ' ' + ', '.join(f'{v:g} dBZ for {lab}' for lab, v in cores.items()) + ')'
    return mask_title(st0)


def xlim_of(args, name):
    key = 'M' if name.startswith('M_') else name
    for v, lo, hi in args.xlim:
        if v == key:
            return float(lo), float(hi)
    return var_info(name)['lim']


def heights_km(st):
    return st.height.values / 1000.


# ====== Quantile CFADs (David et al. 2025) ======
def profile(ax, st, name, color, args, label=None, q95=False, extremes=False, lw=2.):
    """Median (line) and Q25-Q75 (shaded) of a variable; returns the line (None if absent)."""
    if f'count_{name}' not in st:
        return None
    c, e = cfad.get_for_quantiles(st, name)
    q = cfad.quantiles(c, e, [5, 25, 50, 75, 95], args.min_count)
    z = heights_km(st)
    line, = ax.plot(q[:, 2], z, color=color, lw=lw, label=label)
    ax.fill_betweenx(z, q[:, 1], q[:, 3], color=color, alpha=0.2, lw=0)
    if extremes:
        ax.plot(q[:, 0], z, color=color, lw=1, ls='--')
    if extremes or q95:
        ax.plot(q[:, 4], z, color=color, lw=1.2, ls=':')
    return line


def finish_axis(ax, args, xlabel, xlim=None, ylabel=True, letter=None, title=None):
    ax.set_ylim(args.hmin, args.hmax)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel('Height (km)')
    ax.grid(ls='--', lw=0.5, alpha=0.7)
    if title:
        ax.set_title(title, fontsize=11)
    if letter:
        ax.text(0.02, 0.98, f'({letter})', transform=ax.transAxes, va='top', ha='left', fontsize=11,
                fontweight='bold')


def panel_radar(ax, datasets, colors, var, args, letter=None, ylabel=True, legend=True):
    lines = []
    for lab, st in datasets.items():
        ln = profile(ax, st, var, colors[lab], args, label=lab, extremes=args.extremes)
        if ln is not None:
            lines.append(ln)
    finish_axis(ax, args, axis_label(var), xlim_of(args, var), ylabel, letter)
    if legend and lines:
        ax.legend(handles=lines, loc='upper right', fontsize=9)


def panel_contents(ax, st, lab, args, letter=None, ylabel=True, legend=True):
    hyds = [h for h in HYDRO_ORDER if f'count_M_{h}' in st]
    hyds += [n[2:] for n in cfad.names_in(st) if n.startswith('M_') and n[2:] not in hyds]
    if args.hydrometeors:
        hyds = [h for h in hyds if h in args.hydrometeors]
    lines = []
    for h in hyds:
        ln = profile(ax, st, f'M_{h}', HYDRO_COLORS.get(h, 'k'), args, label=HYDRO_NAMES.get(h, h))
        lines.append(ln)
    if args.contents_log:
        ax.set_xscale('log')
    lim = xlim_of(args, 'M_rr')
    if args.contents_log and lim[0] <= 0:
        lim = (1e-3, max(lim[1], 10.))
    finish_axis(ax, args, axis_label('M_rr').replace('rain content', 'Content'), lim, ylabel, letter, lab)
    if legend and lines:
        ax.legend(handles=lines, loc='upper right', fontsize=8)


def panel_fw(ax, datasets, colors, args, letter=None, ylabel=True):
    lines = []
    for lab, st in datasets.items():
        ln = profile(ax, st, 'Fw', colors[lab], args, label=lab, q95=True)
        if ln is not None:
            lines.append(ln)
    finish_axis(ax, args, axis_label('Fw'), xlim_of(args, 'Fw'), ylabel, letter)
    if lines:
        ax.legend(handles=lines, loc='upper right', fontsize=9)


def save(fig, args, suffix):
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, f'{args.name}_{suffix}.png')
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print('  ->', path)


def fig_summary(datasets, colors, variables, args):
    with_contents = [lab for lab, st in datasets.items() if any(n.startswith('M_') for n in cfad.names_in(st))]
    with_fw = any('count_Fw' in st for st in datasets.values())
    row2 = len(with_contents) + (1 if with_fw else 0)
    ncol = max(len(variables), row2, 1)
    nrow = 2 if row2 else 1
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 5.2 * nrow), squeeze=False, sharey=True)
    letters = iter('abcdefghijklmnop')
    for j, var in enumerate(variables):
        panel_radar(axs[0, j], datasets, colors, var, args, next(letters), ylabel=(j == 0))
    for j in range(len(variables), ncol):
        axs[0, j].set_visible(False)
    if row2:
        for j, lab in enumerate(with_contents):
            panel_contents(axs[1, j], datasets[lab], lab, args, next(letters), ylabel=(j == 0))
        if with_fw:
            panel_fw(axs[1, len(with_contents)], datasets, colors, args, next(letters),
                     ylabel=(len(with_contents) == 0))
        for j in range(row2, ncol):
            axs[1, j].set_visible(False)
    fig.suptitle(f'CFADs {mask_text(datasets)}: median and Q25-Q75', fontsize=13)
    fig.tight_layout()
    save(fig, args, 'summary')


def fig_radar(datasets, colors, variables, args):
    fig, axs = plt.subplots(1, len(variables), figsize=(4.2 * len(variables), 5.2), squeeze=False, sharey=True)
    for j, var in enumerate(variables):
        panel_radar(axs[0, j], datasets, colors, var, args, 'abcdefgh'[j], ylabel=(j == 0))
    ext = ', Q5 (dashed), Q95 (dotted)' if args.extremes else ''
    fig.suptitle(f'CFADs {mask_text(datasets)}: median, Q25-Q75{ext}', fontsize=12)
    fig.tight_layout()
    save(fig, args, 'radar')


def fig_contents(datasets, args):
    labs = [lab for lab, st in datasets.items() if any(n.startswith('M_') for n in cfad.names_in(st))]
    if not labs:
        print('  no contents in the statistics (observations only?)')
        return
    fig, axs = plt.subplots(1, len(labs), figsize=(4.5 * len(labs), 5.5), squeeze=False, sharey=True)
    for j, lab in enumerate(labs):
        panel_contents(axs[0, j], datasets[lab], lab, args, 'abcdefgh'[j], ylabel=(j == 0))
    fig.suptitle(f'Hydrometeor contents {mask_text({lab: datasets[lab] for lab in labs})}: median and Q25-Q75',
                 fontsize=12)
    fig.tight_layout()
    save(fig, args, 'contents')


def fig_fw(datasets, colors, args):
    if not any('count_Fw' in st for st in datasets.values()):
        print('  no Fw in the statistics (needs Fw_posg and wet species)')
        return
    fig, ax = plt.subplots(figsize=(4.8, 5.5))
    panel_fw(ax, datasets, colors, args)
    ax.set_title(f'Fw {mask_text(datasets)}\nmedian, Q25-Q75, Q95 (dotted)', fontsize=9)
    fig.tight_layout()
    save(fig, args, 'fw')


# ====== Frequency CFADs ======
def bin_width(args, var):
    for v, w in args.bin:
        if v == var or (v == 'M' and var.startswith('M_')):
            return float(w)
    return var_info(var)['freq_bin']


def frequency(st, var, args):
    """(frequency (height, bin), coarse edges, counts per level)."""
    c, e = cfad.get(st, var)
    total = c.sum(axis=1)
    lim = xlim_of(args, var)
    cc, ce = cfad.rebin(c, e, bin_width(args, var), lim)
    if args.freq_norm == 'level':
        with np.errstate(invalid='ignore', divide='ignore'):
            f = 100. * cc / total[:, None]
    elif args.freq_norm == 'total':
        f = 100. * cc / max(total.sum(), 1)
    else:
        f = cc.astype(float)
    f[total < args.min_count] = np.nan
    return f, ce, total


def fig_frequency(datasets, colors, var, args):
    labs = [lab for lab, st in datasets.items() if f'count_{var}' in st]
    if not labs:
        return
    data = {lab: frequency(datasets[lab], var, args) for lab in labs}
    diff = args.diff and len(labs) > 1
    nrow = 2 if diff else 1
    fig, axs = plt.subplots(nrow, len(labs), figsize=(max(4.3 * len(labs), 6.5), 4.8 * nrow + 0.6), squeeze=False,
                            sharey=True, sharex=True, layout='constrained')
    unit = {'level': '% of the values of the level', 'total': '% of all values', 'count': 'number of values'}
    vmax = np.nanmax([np.nanmax(f) for f, _, _ in data.values()])
    if args.log:
        positive = [np.nanmin(np.where(f > 0, f, np.nan)) for f, _, _ in data.values()]
        norm = LogNorm(vmin=max(np.nanmin(positive), vmax * 1e-4), vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cmap = mpl.colormaps['turbo'].copy()
    cmap.set_bad('white')
    for j, lab in enumerate(labs):
        f, ce, _ = data[lab]
        st = datasets[lab]
        z = heights_km(st)
        dz = z[1] - z[0] if len(z) > 1 else 0.5
        zedges = np.append(z - dz / 2, z[-1] + dz / 2)
        ax = axs[0, j]
        im = ax.pcolormesh(ce, zedges, np.where(f > 0, f, np.nan), cmap=cmap, norm=norm)
        c, e = cfad.get_for_quantiles(st, var)
        q = cfad.quantiles(c, e, [25, 50, 75], args.min_count)
        ax.plot(q[:, 1], z, color='k', lw=1.8)
        ax.plot(q[:, 0], z, color='k', lw=1, ls='--')
        ax.plot(q[:, 2], z, color='k', lw=1, ls='--')
        finish_axis(ax, args, axis_label(var) if not diff else '', (ce[0], ce[-1]), j == 0, title=lab)
    fig.colorbar(im, ax=axs[0, :], label=unit[args.freq_norm], pad=0.01, shrink=0.9)
    if diff:
        ref = labs[0]
        fref = data[ref][0]
        diffs = {lab: data[lab][0] - fref for lab in labs[1:]}
        amax = np.nanmax([np.nanmax(np.abs(d)) for d in diffs.values()]) or 1.
        dnorm = TwoSlopeNorm(vcenter=0., vmin=-amax, vmax=amax)
        axs[1, 0].set_visible(False)
        axs[0, 0].tick_params(labelbottom=True)
        axs[0, 0].set_xlabel(axis_label(var))
        axs[1, 1].tick_params(labelleft=True)
        for j, lab in enumerate(labs[1:], start=1):
            ce = data[lab][1]
            z = heights_km(datasets[lab])
            dz = z[1] - z[0] if len(z) > 1 else 0.5
            zedges = np.append(z - dz / 2, z[-1] + dz / 2)
            imd = axs[1, j].pcolormesh(ce, zedges, diffs[lab], cmap='RdBu_r', norm=dnorm)
            finish_axis(axs[1, j], args, axis_label(var), (ce[0], ce[-1]), j == 1, title=f'{lab} - {ref}')
        fig.colorbar(imd, ax=axs[1, :], label='difference (' + unit[args.freq_norm].replace('number of', '') + ')',
                     pad=0.01, shrink=0.9)
    fig.suptitle(f'Frequency CFAD of {var} {mask_text(datasets)}\n'
                 'median (line) and Q25-Q75 (dashed)', fontsize=12)
    save(fig, args, f'freq_{var}')


def main():
    args = parse_args()
    datasets = load_datasets(args)
    colors = dataset_colors(list(datasets), args.colors)
    present = set().union(*[cfad.names_in(st) for st in datasets.values()])
    variables = args.vars or [v for v in ['Zh', 'Zdr', 'Kdp'] if v in present]
    figs = {'summary', 'radar', 'contents', 'fw', 'freq'} if 'all' in args.figs else set(args.figs)
    if args.name is None:
        st0 = next(iter(datasets.values()))
        args.name = 'CFAD_' + '_'.join(re.sub(r'\W+', '', lab) for lab in datasets) + f"_{st0.attrs.get('mask', '')}"
    for lab, st in datasets.items():
        print(f"{lab}: {st.attrs.get('n_files', '?')} file(s), {st.attrs.get('n_times', '?')} time step(s), "
              f"{st.attrs.get('n_mask_columns', '?')} columns, mask '{st.attrs.get('mask', '?')}'")
    if 'summary' in figs:
        fig_summary(datasets, colors, variables, args)
    if 'radar' in figs:
        fig_radar(datasets, colors, variables, args)
    if 'contents' in figs:
        fig_contents(datasets, args)
    if 'fw' in figs:
        fig_fw(datasets, colors, args)
    if 'freq' in figs:
        for var in variables + (['Fw'] if 'Fw' in present and args.vars is None else []):
            fig_frequency(datasets, colors, var, args)


if __name__ == '__main__':
    main()
