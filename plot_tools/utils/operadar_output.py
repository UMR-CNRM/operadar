#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read operadar output files (and gridded radar observations) for the analysis scripts of plot_tools:
compute_cfad.py, plot_cfad.py, plot_output_check.py, plot_cross_section.py, plot_compare_outputs.py and
plot_max_over_time.py.

Two kinds of files are read:
- "model": operadar NetCDF outputs (dpolvar_<model>_<micro>_<band>band_<date>.nc), with dimensions
  (level, y, x), the altitude of each grid point in `Alt` (m) and the contents in `Contents` (g/m3, one value per
  hydrometeor). Model levels are interpolated on regular heights with `HeightGrid`.
- "grid": data already on regular heights, e.g. radar observations interpolated on a Cartesian grid (Py-ART),
  with a vertical dimension `z` (or `height`, `alt`, `altitude`) in metres, and possibly a `time` dimension.
  Lower-case names (zh, zdr, kdp, rhohv) and `M` for the contents are accepted.

Only numpy, xarray and scipy are needed (no cartopy, no epygram).
"""
import glob
import os
import sys

import numpy as np
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage


# ====== Variables ======
RADAR_VARS = ['Zh', 'Zdr', 'Kdp', 'Rhohv', 'Ah', 'Av', 'Zh_att']

# Other names accepted in the input files (observations, older files)
ALIASES = {'zh': 'Zh', 'ZH': 'Zh', 'dbz': 'Zh', 'DBZ': 'Zh',
           'zdr': 'Zdr', 'ZDR': 'Zdr',
           'kdp': 'Kdp', 'KDP': 'Kdp',
           'rhohv': 'Rhohv', 'RHOHV': 'Rhohv', 'rhv': 'Rhohv',
           'M': 'Contents',
           }
HEIGHT_DIMS = ['z', 'height', 'alt', 'altitude']

# label, unit, histogram (min, max, bin width; contents: 0 then 500 log bins from 1e-5 to 30 g/m3, i.e. 3 %)
# used by compute_cfad.py, default plot range,
# bin width of the frequency CFADs, colour bounds of the maps (min, max, step)
VAR_INFO = {
    'Zh':     dict(label=r'$Z_H$', unit='dBZ', hist=(-30., 80., 0.1), lim=(10., 60.), freq_bin=2.,
                   bounds=(8., 64., 4.)),
    'Zh_att': dict(label=r'$Z_H$ attenuated', unit='dBZ', hist=(-30., 80., 0.1), lim=(10., 60.), freq_bin=2.,
                   bounds=(8., 64., 4.)),
    'Zdr':    dict(label=r'$Z_{DR}$', unit='dB', hist=(-3., 10., 0.002), lim=(-0.5, 3.), freq_bin=0.25,
                   bounds=(-0.5, 5., 0.5)),
    'Kdp':    dict(label=r'$K_{DP}$', unit=r'° km$^{-1}$', hist=(-3., 20., 0.002), lim=(-0.5, 3.), freq_bin=0.25,
                   bounds=(0., 6., 0.5)),
    'Rhohv':  dict(label=r'$\rho_{HV}$', unit='-', hist=(0.5, 1.002, 0.0005), lim=(0.9, 1.0), freq_bin=0.005,
                   bounds=(0.9, 1., 0.01)),
    'Ah':     dict(label=r'$A_H$', unit=r'dB km$^{-1}$', hist=(0., 20., 0.001), lim=(0., 0.5), freq_bin=0.02,
                   bounds=(0., 1., 0.1)),
    'Av':     dict(label=r'$A_V$', unit=r'dB km$^{-1}$', hist=(0., 20., 0.001), lim=(0., 0.5), freq_bin=0.02,
                   bounds=(0., 1., 0.1)),
    'M':      dict(label='Content', unit=r'g m$^{-3}$', hist=('log', 1e-5, 30., 500), lim=(0., 3.), freq_bin=0.1,
                   bounds=(0., 5., 0.5)),
    'Fw':     dict(label=r'$F_w$ (wet graupel/hail)', unit='-', hist=(0., 1.0001, 0.001), lim=(0., 1.),
                   freq_bin=0.05, bounds=(0., 1., 0.1)),
}

HYDRO_NAMES = {'rr': 'rain', 'cc': 'cloud water', 'cl': 'cloud water', 'cs': 'cloud water (sea)',
               'ii': 'pristine ice', 'ss': 'snow', 'gg': 'graupel', 'wg': 'wet graupel',
               'hh': 'hail', 'wh': 'wet hail'}
# Colours close to the CFADs of David et al. (2025, AMT)
HYDRO_COLORS = {'rr': 'navy', 'cc': 'tab:cyan', 'cl': 'tab:cyan', 'cs': 'tab:cyan', 'ii': 'tab:olive',
                'ss': 'tab:pink', 'gg': 'tab:gray', 'wg': 'tab:purple', 'hh': 'tab:brown', 'wh': 'tab:red'}


def var_info(name: str) -> dict:
    """Metadata of a variable: radar variable, content 'M_<hydrometeor>' or 'Fw'."""
    if name.startswith('M_'):
        info = dict(VAR_INFO['M'])
        info['label'] = f"{HYDRO_NAMES.get(name[2:], name[2:])} content"
        return info
    return VAR_INFO.get(name, dict(label=name, unit='', hist=(-100., 100., 0.1), lim=None, freq_bin=1.,
                                   bounds=(0., 10., 1.)))


def axis_label(name: str) -> str:
    info = var_info(name)
    return f"{info['label']} ({info['unit']})" if info['unit'] not in ('', '-') else info['label']


# Discrete radar-like colour map (no epygram needed); white under the lowest bound, pink over the highest
RADAR_CMAP = LinearSegmentedColormap.from_list(
    'radar_operadar', ['#a6dcff', '#3c97e0', '#0b4fbf', '#08a14a', '#3fd13a', '#b5e837', '#ffe733',
                       '#ffb02e', '#f5671f', '#d7191c', '#99000d', '#c51b8a', '#7a0177'])
RADAR_CMAP.set_under('white')
RADAR_CMAP.set_over('deeppink')
RADAR_CMAP.set_bad('white')


# ====== Files ======
def expand_files(patterns: list[str]) -> list[str]:
    """Expand the shell patterns (quoted globs are allowed) into a sorted list of existing files."""
    files = []
    for p in patterns:
        found = sorted(glob.glob(os.path.expanduser(p)))
        if not found and os.path.isfile(p):
            found = [p]
        if not found:
            print(f"/!\\ no file matches {p}")
        files += [f for f in found if f not in files]
    if not files:
        sys.exit("No input file found")
    return files


def open_output(path: str, box: list[float] | None = None) -> xr.Dataset:
    """Open an operadar output or a gridded file, harmonise the names, optionally cut a lon/lat box.

    The kind of file ('model' or 'grid') is stored in ds.attrs['_kind'].
    box = [lon_min, lon_max, lat_min, lat_max] (needs the lat/lon coordinates, i.e. real_case=True).
    """
    ds = xr.open_dataset(path)
    rename = {k: v for k, v in ALIASES.items() if k in ds.variables and v not in ds.variables}
    for d in HEIGHT_DIMS[1:]:
        if d in ds.dims and 'z' not in ds.dims:
            rename[d] = 'z'
    if rename:
        ds = ds.rename(rename)
    if 'level' in ds.dims and 'Alt' in ds:
        ds.attrs['_kind'] = 'model'
    elif 'z' in ds.dims:
        ds.attrs['_kind'] = 'grid'
        if float(ds.z.max()) < 100.:  # heights given in km
            ds = ds.assign_coords(z=ds.z * 1000.)
    else:
        raise ValueError(f"{path}: no 'level' dimension with 'Alt', and no height dimension {HEIGHT_DIMS}")
    if box is not None:
        ds = cut_box(ds, box)
    return ds


def cut_box(ds: xr.Dataset, box: list[float]) -> xr.Dataset:
    if 'lat' not in ds.variables or 'lon' not in ds.variables:
        sys.exit("--box needs lat/lon in the files (operadar configuration with real_case=True)")
    lon_min, lon_max, lat_min, lat_max = box
    lat, lon = ds.lat.values, ds.lon.values
    if lat.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    inside = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    if not inside.any():
        sys.exit(f"--box {box}: no grid point inside")
    iy, ix = np.nonzero(inside)
    ydim, xdim = ds.lat.dims[-2:] if ds.lat.ndim == 2 else ('y', 'x')
    return ds.isel({ydim: slice(iy.min(), iy.max() + 1), xdim: slice(ix.min(), ix.max() + 1)})


def vertical_dim(ds) -> str:
    return 'level' if ds.attrs.get('_kind') == 'model' else 'z'


def time_steps(ds: xr.Dataset):
    """Yield the dataset at each time step (operadar files have a single, scalar time)."""
    if 'time' in ds.dims:
        for i in range(ds.sizes['time']):
            yield ds.isel(time=i)
    else:
        yield ds


def time_string(ds) -> str:
    if 'time' in ds.coords and ds.time.size == 1:
        try:
            return np.datetime_as_string(ds.time.values, unit='m').replace('T', ' ')
        except (TypeError, ValueError):
            return str(ds.time.values)
    return ''


def describe(ds) -> str:
    """Short description: model, microphysics, band (from the operadar attributes)."""
    parts = [str(ds.attrs[k]) for k in ('model', 'microphysics') if k in ds.attrs]
    if 'radar_band' in ds.attrs:
        parts.append(f"{ds.attrs['radar_band']} band")
    return ' '.join(parts)


def horizontal_km(ds):
    """1D x and y coordinates in km (grid index times the resolution if x/y are missing)."""
    res = float(ds.attrs.get('horizontal_resolution', 1000.))
    x = ds.x.values / 1000. if 'x' in ds.coords else np.arange(ds.sizes['x']) * res / 1000.
    y = ds.y.values / 1000. if 'y' in ds.coords else np.arange(ds.sizes['y']) * res / 1000.
    return x - x[0], y - y[0]


# ====== Vertical interpolation ======
class HeightGrid:
    """Values of 3D fields on regular heights, for a selection of columns.

    - model files: linear interpolation between the model levels, using the altitude `Alt` of each grid point;
      reference='agl' gives heights above the lowest model level (close to the height above ground),
      reference='asl' gives altitudes above sea level.
    - grid files: nearest height of the `z` coordinate (within half a height step).

    columns: 2D boolean array (y, x) of the columns to keep (None: all). Calling the object on a 3D DataArray
    returns an array (n_heights, n_columns), NaN where the height is outside the column.
    """

    def __init__(self, ds: xr.Dataset, heights: np.ndarray, reference: str = 'agl', columns=None):
        self.heights = np.asarray(heights, dtype='f4')
        self.kind = ds.attrs.get('_kind')
        self.columns = None if columns is None else np.asarray(columns, bool)
        if self.kind == 'model':
            alt = self._columns(ds['Alt'].values.astype('f4'))
            self.flip = alt.shape[0] > 1 and np.nanmean(alt[0]) > np.nanmean(alt[-1])  # AROME: level 0 = top
            if self.flip:
                alt = alt[::-1]
            if reference == 'agl':
                alt = alt - alt[0:1]
            nl = alt.shape[0]
            nz, ncol = len(self.heights), alt.shape[1]
            self.k = np.zeros((nz, ncol), np.intp)
            self.w = np.zeros((nz, ncol), 'f4')
            self.valid = np.zeros((nz, ncol), bool)
            for i, h in enumerate(self.heights):
                below = (alt <= h).sum(axis=0) - 1
                valid = (below >= 0) & (below < nl - 1)
                kc = np.clip(below, 0, nl - 2)
                a0 = np.take_along_axis(alt, kc[None], 0)[0]
                a1 = np.take_along_axis(alt, kc[None] + 1, 0)[0]
                with np.errstate(divide='ignore', invalid='ignore'):
                    self.w[i] = np.where(a1 > a0, (h - a0) / (a1 - a0), 0.)
                self.k[i], self.valid[i] = kc, valid
        else:
            z = ds.z.values.astype('f8')
            dz = np.median(np.diff(np.sort(z))) if z.size > 1 else 1e9
            self.zidx = np.array([np.argmin(np.abs(z - h)) for h in self.heights])
            self.zok = np.abs(z[self.zidx] - self.heights) <= dz / 2 + 1e-3

    def _columns(self, field: np.ndarray) -> np.ndarray:
        """(nv, y, x) -> (nv, n_columns)."""
        field = field.reshape(field.shape[0], -1)
        if self.columns is not None:
            field = field[:, self.columns.ravel()]
        return field

    def __call__(self, da: xr.DataArray) -> np.ndarray:
        vdim = 'level' if self.kind == 'model' else 'z'
        da = da.transpose(vdim, ...)
        if self.kind == 'model':
            f = self._columns(da.values.astype('f4'))
            if self.flip:
                f = f[::-1]
            f0 = np.take_along_axis(f, self.k, 0)
            f1 = np.take_along_axis(f, self.k + 1, 0)
            out = f0 + (f1 - f0) * self.w
            out[~self.valid] = np.nan
            return out
        out = self._columns(da.values.astype('f4'))[self.zidx]
        out[~self.zok] = np.nan
        return out


def column_max(da: xr.DataArray, vdim: str) -> np.ndarray:
    """Maximum over the vertical (NaN ignored)."""
    return da.max(dim=vdim, skipna=True).values


# ====== Storm cells ======
def cell_masks(zmax: np.ndarray, envelope_dbz=25., core_dbz=40., min_envelope=20, min_core=5):
    """Storm cell envelopes and cores from the column maximum of Zh, as in David et al. (2025, AMT, Sect. 3.3):
    envelopes = contiguous areas of Zmax >= envelope_dbz with at least min_envelope grid points,
    cores = contiguous areas of Zmax >= core_dbz with at least min_core grid points, inside an envelope.
    Returns two 2D boolean arrays (envelope, core) and the number of cores.
    """
    structure = np.ones((3, 3), bool)  # 8-connectivity
    zmax = np.nan_to_num(zmax, nan=-999.)

    def keep_large(binary, nmin):
        lab, n = ndimage.label(binary, structure=structure)
        if n == 0:
            return np.zeros_like(binary), 0
        sizes = np.bincount(lab.ravel())
        keep = sizes >= nmin
        keep[0] = False
        return keep[lab], int(keep.sum())

    envelope, _ = keep_large(zmax >= envelope_dbz, min_envelope)
    core, ncore = keep_large((zmax >= core_dbz) & envelope, min_core)
    return envelope, core, ncore


# ====== Liquid water fraction ======
def wet_fraction(ds: xr.Dataset):
    """Liquid water fraction Fw of the wet graupel (and wet hail), rebuilt from the saved contents.

    With mixed_phase_parametrization='Fw_posg', the rain is kept and the wet graupel/hail holds the former
    graupel/hail content in the melting layer, so Fw = rr / (rr + wg + wh) where wg + wh > 0.
    Returns (DataArray or None, message).
    """
    if 'Contents' not in ds:
        return None, 'no Contents in the file'
    mp = ds.attrs.get('mixed_phase_type', '')
    hyd = [str(h) for h in ds.hydrometeor.values]
    wets = [h for h in ('wg', 'wh') if h in hyd]
    if not wets or 'rr' not in hyd:
        return None, 'no wet species in the file'
    if mp != 'Fw_posg':
        return None, (f"mixed phase '{mp}': Fw can only be rebuilt from the contents with 'Fw_posg' "
                      "(with 'Fw_pos' the rain of the melting layer is moved into the wet species)")
    rr = ds.Contents.sel(hydrometeor='rr')
    wet = sum(ds.Contents.sel(hydrometeor=h) for h in wets)
    return (rr / (rr + wet)).where(wet > 0), ''
