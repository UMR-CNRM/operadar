#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistics for the CFADs of operadar outputs (and of gridded radar observations), as in
David et al. (2025, AMT, Fig. 6 and 8): distribution of Zh, Zdr, Kdp (Rhohv...), of the hydrometeor contents and
of the liquid water fraction of the wet graupel, at each height, inside storm cells.

For each height (regular grid, 0 to 15 km every 500 m by default) the values of all files and time steps are
accumulated in histograms with fine bins, and saved in a small NetCDF file. Plot it with plot_cfad.py
(quantile CFADs as in the paper, and frequency CFADs).

Where the values are taken (--mask):
  core      storm cores (default): contiguous areas where the column maximum of Zh >= --core-dbz (40 dBZ)
            with at least --min-core (5) grid points, inside an envelope;
            the paper uses 40 dBZ for simulations and 48 dBZ for observations
  envelope  storm envelopes: column max of Zh >= --envelope-dbz (25 dBZ), at least --min-envelope (20) points
  zh        every grid point where Zh >= --zh-local (3D criterion)
  all       every grid point
  var       where the variable --mask-var of the files is > 0 (2D or 3D), e.g. a cell or ZDR column mask
            computed elsewhere (tracking)
With --cells-from, the cores and envelopes are those of other files on the same grid (same storm cells for
several radar bands).
The radar variables are only used where Zh >= --zh-min (0 dBZ by default), the contents everywhere in the mask
(zeros included, as in the paper), Fw where there is wet graupel/hail.

Examples (from plot_tools/):
  python compute_cfad.py "../operadarFiles/AROME_ICE3/dpolvar_*.nc" --label "AROME ICE3" --out CFAD/ICE3_core.nc
  python compute_cfad.py "../operadarFiles/AROME_LIMA/*.nc" --label "AROME LIMA" --out CFAD/LIMA_core.nc
  python compute_cfad.py "OBS/*.nc" --label obs --core-dbz 48 --out CFAD/obs_core.nc
  python compute_cfad.py "../operadarFiles/*.nc" --mask zh --zh-local 40 --label ICE3 --out CFAD/ICE3_zh40.nc
  python compute_cfad.py W.nc --cells-from C.nc --label "W band, C band cores" --out CFAD/W_core.nc
"""
import argparse
import os
import sys
import time

import numpy as np

from utils.operadar_output import (RADAR_VARS, HeightGrid, cell_masks, column_max, describe, expand_files,
                                   open_output, time_steps, vertical_dim, wet_fraction)
from utils.cfad import CfadAccumulator


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', help='operadar output files or gridded files (shell patterns allowed)')
    p.add_argument('--out', required=True, help='output NetCDF file of statistics')
    p.add_argument('--label', default=None, help='name of the dataset in the legends (default: model micro)')
    p.add_argument('--mask', default='core', choices=['core', 'envelope', 'zh', 'all', 'var'])
    p.add_argument('--cells-from', nargs='+', default=None, metavar='FILE',
                   help='take the cores/envelopes (--mask core or envelope) from these files instead (same grid, '
                        'one per input file, shell patterns allowed), e.g. the C band output of the same run '
                        'for a W band output, whose Zh rarely reaches 40 dBZ')
    p.add_argument('--core-dbz', type=float, default=40.)
    p.add_argument('--envelope-dbz', type=float, default=25.)
    p.add_argument('--min-core', type=int, default=5)
    p.add_argument('--min-envelope', type=int, default=20)
    p.add_argument('--zh-local', type=float, default=40., help='threshold of --mask zh (dBZ)')
    p.add_argument('--mask-var', default=None, help='variable of the files used with --mask var')
    p.add_argument('--zh-var', default='Zh', help='reflectivity used for the cells and --zh-min (Zh or Zh_att)')
    p.add_argument('--zh-min', type=float, default=0., help='radar variables only where Zh >= zh-min (dBZ)')
    p.add_argument('--vars', nargs='+', default=None,
                   help=f'radar variables (default: those of {RADAR_VARS[:4]} found in the files)')
    p.add_argument('--no-contents', action='store_true', help='do not compute the CFADs of the contents')
    p.add_argument('--no-fw', action='store_true', help='do not compute the CFAD of Fw')
    p.add_argument('--hmin', type=float, default=0., help='lowest height (km)')
    p.add_argument('--hmax', type=float, default=15., help='highest height (km)')
    p.add_argument('--dz', type=float, default=0.5, help='height step (km)')
    p.add_argument('--height', default='agl', choices=['agl', 'asl'],
                   help='agl: above the lowest model level (~ground, as in the paper); asl: above sea level')
    p.add_argument('--box', nargs=4, type=float, default=None, metavar=('LONMIN', 'LONMAX', 'LATMIN', 'LATMAX'),
                   help='use only this lon/lat box')
    return p.parse_args()


def selected_columns(args, ds, vdim, zh, zh_cells=None):
    """2D mask of the columns to use and number of cells (or None). zh_cells: reflectivity defining the cells
    (default zh)."""
    if args.mask in ('core', 'envelope'):
        if zh_cells is not None:
            zh = zh_cells
        envelope, core, ncore = cell_masks(column_max(zh, vdim), args.envelope_dbz, args.core_dbz,
                                           args.min_envelope, args.min_core)
        return (core if args.mask == 'core' else envelope), ncore
    if args.mask == 'zh':
        return np.nan_to_num(column_max(zh, vdim), nan=-999.) >= args.zh_local, None
    if args.mask == 'var':
        m = ds[args.mask_var]
        if vdim in m.dims:
            return (m.fillna(0) > 0).any(dim=vdim).values, None
        return m.fillna(0).values > 0, None
    return np.ones(zh.shape[-2:], bool), None


def point_mask(args, ds, vdim, grid, columns):
    """3D mask (n_heights, n_columns) inside the selected columns."""
    nz, ncol = len(grid.heights), int(columns.sum())
    if args.mask == 'zh':
        return np.nan_to_num(grid(ds[args.zh_var]), nan=-999.) >= args.zh_local
    if args.mask == 'var' and vdim in ds[args.mask_var].dims:
        return np.nan_to_num(grid(ds[args.mask_var])) > 0
    return np.ones((nz, ncol), bool)


def main():
    args = parse_args()
    if args.mask == 'var' and not args.mask_var:
        sys.exit('--mask var needs --mask-var NAME')
    files = expand_files(args.files)
    cell_files = expand_files(args.cells_from) if args.cells_from else None
    if cell_files is not None:
        if args.mask not in ('core', 'envelope'):
            sys.exit('--cells-from needs --mask core or envelope')
        if len(cell_files) != len(files):
            sys.exit(f'--cells-from: {len(cell_files)} file(s) for {len(files)} input file(s)')
    heights = np.arange(args.hmin * 1000., args.hmax * 1000. + 1., args.dz * 1000.)

    # Variables available in the first file
    ds0 = open_output(files[0], args.box)
    radar_vars = args.vars or [v for v in RADAR_VARS[:4] if v in ds0]
    missing = [v for v in radar_vars + [args.zh_var] if v not in ds0]
    if missing:
        sys.exit(f'{files[0]}: missing variables {missing}')
    hydros = [] if args.no_contents or 'Contents' not in ds0 else [str(h) for h in ds0.hydrometeor.values]
    fw0, why = (None, '') if args.no_fw else wet_fraction(ds0)
    do_fw = fw0 is not None
    if not args.no_fw and not do_fw:
        print(f'No CFAD of Fw: {why}')
    names = radar_vars + [f'M_{h}' for h in hydros] + (['Fw'] if do_fw else [])
    attrs = dict(label=args.label or describe(ds0) or os.path.basename(files[0]),
                 description=describe(ds0), source_kind=ds0.attrs['_kind'],
                 mask=args.mask, zh_var=args.zh_var, zh_min=args.zh_min, height_reference=args.height)
    if args.mask in ('core', 'envelope'):
        attrs.update(core_dbz=args.core_dbz, envelope_dbz=args.envelope_dbz,
                     min_core=args.min_core, min_envelope=args.min_envelope)
        if cell_files:
            attrs['cells_from'] = ' '.join(os.path.basename(f) for f in cell_files)[:20000]
    elif args.mask == 'zh':
        attrs.update(zh_local=args.zh_local)
    elif args.mask == 'var':
        attrs.update(mask_var=args.mask_var)
    for k in ('model', 'microphysics', 'radar_band', 'mixed_phase_type', 'scattering_method'):
        if k in ds0.attrs:
            attrs[k] = str(ds0.attrs[k])
    ds0.close()
    print(f"{attrs['label']}: {len(files)} file(s), mask '{args.mask}', variables {names}")

    acc = CfadAccumulator(heights, names)
    n_times = n_columns = n_cells = 0
    t0 = time.time()
    for ifile, path in enumerate(files):
        ds_file = open_output(path, args.box)
        vdim = vertical_dim(ds_file)
        ds_cells = open_output(cell_files[ifile], args.box) if cell_files else None
        if ds_cells is not None and ds_cells[args.zh_var].shape != ds_file[args.zh_var].shape:
            sys.exit(f'--cells-from: {cell_files[ifile]} is not on the grid of {path}')
        cells_steps = time_steps(ds_cells) if ds_cells is not None else None
        for ds in time_steps(ds_file):
            n_times += 1
            zh = ds[args.zh_var]
            zh_cells = next(cells_steps)[args.zh_var] if cells_steps is not None else None
            columns, ncell = selected_columns(args, ds, vdim, zh, zh_cells)
            n_cells += ncell or 0
            ncol = int(columns.sum())
            if ncol == 0:
                continue
            n_columns += ncol
            grid = HeightGrid(ds, heights, args.height, columns)
            inside = point_mask(args, ds, vdim, grid, columns)
            zh_h = np.nan_to_num(grid(zh), nan=-999.)
            radar_ok = inside & (zh_h >= args.zh_min)
            for v in radar_vars:
                vals = zh_h if v == args.zh_var else grid(ds[v])
                acc.add(v, np.where(radar_ok, vals, np.nan))
            for h in hydros:
                acc.add(f'M_{h}', np.where(inside, grid(ds.Contents.sel(hydrometeor=h)), np.nan))
            if do_fw:
                fw, _ = wet_fraction(ds)
                acc.add('Fw', np.where(inside, grid(fw), np.nan))
        ds_file.close()
        if ds_cells is not None:
            ds_cells.close()
        print(f'  [{ifile + 1}/{len(files)}] {os.path.basename(path)}: {n_columns} columns so far '
              f'({time.time() - t0:.0f} s)')

    attrs.update(n_files=len(files), n_times=n_times, n_mask_columns=n_columns,
                 files=' '.join(os.path.basename(f) for f in files)[:20000])
    if args.mask in ('core', 'envelope'):
        attrs['n_cores'] = n_cells
    if n_columns == 0:
        print('/!\\ no column selected: check the thresholds (--mask, --core-dbz...)')
    out = acc.to_dataset(attrs)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_netcdf(args.out, encoding={v: {'zlib': True, 'complevel': 4} for v in out.data_vars})
    total = {n: int(acc.counts[n].sum() + acc.zeros[n].sum()) for n in names[:4]}
    print(f'Saved {args.out}: {n_times} time step(s), {n_columns} columns, values per variable {total}')


if __name__ == '__main__':
    main()
