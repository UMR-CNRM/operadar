#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reading and writing of the lookup tables (NetCDF format).

The Fortran programs write text tables:
  * TmatCoefDiff_<band><type>             scattering coefficients of one particle (Tmat),
                                          rows ordered Tc > ELEV > Fw > D, columns separated by spaces
  * TmatCoefInt_<micro>_<n>M_<band><type> integrated over the size distribution (TmatInt),
                                          rows ordered Tc > ELEV > P3 > M, columns separated by ';'
Both have 3 header lines: parameter names, parameter values, column names.

generate_tables.sh converts them to NetCDF (same name + '.nc') and deletes the text files.
In the NetCDF file each column is a variable on the regular grid of the table:
  * TmatCoefDiff: dims (Tc, ELEV, Fw, D)
  * TmatCoefInt : dims (Tc, ELEV, Fw, M) for 1-moment species and (Tc, ELEV, N, M) for 2-moment species
    (and pristine ice, whose third parameter is a concentration), M in kg/m3, N in m-3.
The header parameters (Tcmin, ELEVstep, expMmax, LAM, canting, ...) are global attributes, and the three
original header lines are kept (attributes header_line1/2/3) so that a text file can be written back
(needed by TmatInt, which reads TmatCoefDiff text files).

Examples
--------
>>> import xarray as xr
>>> ds = xr.open_dataset("TmatCoefInt_ICE3_1M_Crr.nc")
>>> ds.zdr.sel(Tc=10, ELEV=0, Fw=0).plot(x="M", xscale="log")

Command line:
    python table_io.py to-netcdf TABLE [TABLE ...] [--keep-text]
    python table_io.py to-text TABLE.nc [OUTPUT]      (TmatCoefDiff only)
"""

import os
import sys
import numpy as np
import pandas as pd

# The Fortran programs write 'Infinity' (rhohv of cells with an extremely small signal, see the wiki):
# read as NaN, as operadar always did
_TEXT_OPTIONS = dict(skipinitialspace=True, na_values=["Infinity", "-Infinity", "NaN"])


# ----------------------------------------------------------------------------------------------
# Paths and header
# ----------------------------------------------------------------------------------------------
def table_path(path: str) -> str:
    """Return the existing file for a table given with or without '.nc'
    (the NetCDF file is preferred, the text file is the fallback for tables not converted)."""
    base = path[:-3] if path.endswith(".nc") else path
    for candidate in (base + ".nc", base):
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError(f"Lookup table not found: {base}(.nc)")


def is_netcdf(path: str) -> bool:
    with open(path, "rb") as f:
        magic = f.read(4)
    return magic[:3] == b"CDF" or magic == b"\x89HDF"


def _convert(value: str):
    value = value.strip()
    try:
        return float(value)
    except ValueError:
        return value


def _read_text_header(path: str):
    with open(path) as f:
        lines = [f.readline().rstrip("\n") for _ in range(3)]
    sep = ";" if ";" in lines[0] else None
    names = lines[0].split(sep)
    values = lines[1].split(sep)
    params = {n.strip(): _convert(v) for n, v in zip(names, values)}
    columns = [c.strip() for c in lines[2].split(sep) if c.strip()]
    kind = "Int" if sep == ";" else "Diff"
    return kind, sep, lines, params, columns


def _third_parameter_name(path: str, params: dict) -> str:
    """Name of the third dimension of an integrated table: Fw (1 moment) or N (2 moments, pristine ice)."""
    name = os.path.basename(path)
    name = name[:-3] if name.endswith(".nc") else name
    if "_2M_" in name or name.endswith("ii"):
        return "N"
    return "Fw"


# ----------------------------------------------------------------------------------------------
# Text -> NetCDF
# ----------------------------------------------------------------------------------------------
def _grid_shape(kind: str, path: str, params: dict) -> tuple[list[str], list[int]]:
    """Dimensions and size of the table grid from the header parameters."""
    def n(prefix):
        return int(round((params[prefix + "max"] - params[prefix + "min"]) / params[prefix + "step"])) + 1
    if kind == "Diff":
        return ["Tc", "ELEV", "Fw", "D"], [n("Tc"), n("ELEV"), n("Fw"), n("expD")]
    p3 = _third_parameter_name(path, params)
    return ["Tc", "ELEV", p3, "M"], [n("Tc"), n("ELEV"), n("Fw") if p3 == "Fw" else n("expCC"), n("expM")]


def text_to_netcdf(path: str, ncpath: str | None = None, complevel: int = 4) -> str:
    """Convert a text table (TmatCoefDiff or TmatCoefInt) to NetCDF. Returns the NetCDF path.
    The file is read by blocks of temperatures, so the memory used is about the size of the NetCDF
    arrays (4 bytes per value)."""
    import xarray as xr

    kind, sep, lines, params, columns = _read_text_header(path)
    ncpath = ncpath or path + ".nc"
    dims, shape = _grid_shape(kind, path, params)
    grid = ["Tc", "ELEV", "Fw", "D"] if kind == "Diff" else ["Tc", "ELEV", "P3", "M"]
    block = int(np.prod(shape[1:]))  # rows of one temperature
    values = {c: np.empty(shape, dtype=np.float32) for c in columns if c not in grid}
    coords = {d: None for d in dims}
    reader = pd.read_csv(path, sep=r"\s+" if sep is None else sep, skiprows=3, header=None, names=columns,
                         dtype=np.float64, engine="c", chunksize=block, **_TEXT_OPTIONS)
    it = -1
    for it, chunk in enumerate(reader):
        if it >= shape[0] or len(chunk) != block:
            raise ValueError(f"{path}: more rows than the {' x '.join(map(str, shape))} grid of the header, "
                             "or a temperature with missing rows (incomplete table?)")
        sub = shape[1:]
        # grid columns: the rows must be ordered Tc > ELEV > P3 > M (or D) on the grid of the header
        tc = chunk[grid[0]].to_numpy()
        if not np.all(tc == tc[0]):
            raise ValueError(f"{path}: rows are not ordered by temperature blocks")
        if it == 0:
            coords[dims[0]] = []
            for axis, (col, dim) in enumerate(zip(grid[1:], dims[1:])):
                index = [0, 0, 0]
                index[axis] = slice(None)
                coords[dim] = chunk[col].to_numpy().reshape(sub)[tuple(index)].copy()
        for axis, (col, dim) in enumerate(zip(grid[1:], dims[1:])):
            expected = np.expand_dims(coords[dim], [a for a in range(3) if a != axis])
            if not np.array_equal(chunk[col].to_numpy().reshape(sub), np.broadcast_to(expected, sub)):
                raise ValueError(f"{path}: rows are not ordered {' > '.join(grid)} on a regular grid")
        coords[dims[0]].append(tc[0])
        for c in values:
            v = chunk[c].to_numpy().reshape(sub)
            v32 = v.astype(np.float32)
            lost = ((v32 == 0) & (v != 0)) | (np.isinf(v32) & np.isfinite(v))
            if lost.any() and values[c].dtype == np.float32:
                values[c] = values[c].astype(np.float64)  # float32 cannot store a value: keep float64
            values[c][it] = v if values[c].dtype == np.float64 else v32
    if it + 1 != shape[0]:
        raise ValueError(f"{path}: {it + 1} temperatures, {shape[0]} expected from the header (incomplete table?)")
    coords[dims[0]] = np.array(coords[dims[0]])

    encoding = {c: {"zlib": True, "complevel": complevel, "shuffle": True} for c in values}
    attrs = dict(params)
    attrs.update(header_line1=lines[0], header_line2=lines[1], header_line3=lines[2],
                 table_kind="TmatCoef" + kind, columns=" ".join(columns),
                 source_text_file=os.path.basename(path))
    ds = xr.Dataset({c: (dims, v) for c, v in values.items()}, coords=coords, attrs=attrs)
    units = {"Tc": "degC", "ELEV": "degree", "Fw": "1", "N": "m-3", "M": "kg m-3", "D": "mm"}
    for dim in dims:
        ds[dim].attrs["units"] = units[dim]
    tmp = ncpath + ".part"
    ds.to_netcdf(tmp, encoding=encoding, format="NETCDF4")
    os.replace(tmp, ncpath)
    return ncpath


# ----------------------------------------------------------------------------------------------
# NetCDF -> text (TmatCoefDiff, read by TmatInt with fixed formats)
# ----------------------------------------------------------------------------------------------
def _efmt(width: int, digits: int):
    def fmt(x):
        if not np.isfinite(x):
            return ("NaN" if np.isnan(x) else ("Infinity" if x > 0 else "-Infinity")).rjust(width)
        return f"{x:{width}.{digits - 1}E}"
    return fmt


def netcdf_to_text(ncpath: str, path: str | None = None) -> str:
    """Write a TmatCoefDiff NetCDF table back to the text format read by TmatInt
    (Fortran format 5002: F5.1,F4.1,F4.2, 6 x E10.4, 20 x E12.5, separated by 2 spaces)."""
    import xarray as xr

    path = path or (ncpath[:-3] if ncpath.endswith(".nc") else ncpath + ".txt")
    with xr.open_dataset(ncpath) as ds:
        if ds.attrs.get("table_kind") != "TmatCoefDiff":
            raise ValueError("netcdf_to_text only writes TmatCoefDiff tables")
        columns = ds.attrs["columns"].split()
        shape = [ds.sizes[d] for d in ("Tc", "ELEV", "Fw", "D")]
        grids = np.meshgrid(ds.Tc.values, ds.ELEV.values, ds.Fw.values, ds.D.values, indexing="ij")
        cols = dict(zip(["Tc", "ELEV", "Fw", "D"], [g.ravel() for g in grids]))
        for c in columns:
            if c not in cols:
                cols[c] = ds[c].values.astype(np.float64).ravel()
        header = [ds.attrs["header_line1"], ds.attrs["header_line2"], ds.attrs["header_line3"]]
    n = int(np.prod(shape))
    e10, e12 = _efmt(10, 4), _efmt(12, 5)
    fmts = [lambda x: f"{x:5.1f}", lambda x: f"{x:4.1f}", lambda x: f"{x:4.2f}"] + [e10] * 6 + [e12] * (len(columns) - 9)
    tmp = path + ".part"
    with open(tmp, "w") as f:
        f.write("\n".join(header) + "\n")
        arrays = [cols[c] for c in columns]
        for i in range(n):
            f.write("  ".join(fmt(a[i]) for fmt, a in zip(fmts, arrays)) + "\n")
    os.replace(tmp, path)
    return path


# ----------------------------------------------------------------------------------------------
# Reading (operadar and plot tools)
# ----------------------------------------------------------------------------------------------
def read_table_params(path: str) -> dict:
    """Header parameters of a table (Tcmin, Tcstep, ..., LAM, canting ...), NetCDF or text."""
    path = table_path(path)
    if is_netcdf(path):
        import xarray as xr
        with xr.open_dataset(path) as ds:
            return {k: (v.item() if hasattr(v, "item") else v) for k, v in ds.attrs.items()
                    if not k.startswith("header_line") and k not in ("columns", "table_kind", "source_text_file")}
    return _read_text_header(path)[3]


def read_table_columns(path: str, columns: list[str] | None = None, inf_to_nan: bool = True) -> dict:
    """Read the requested columns of a table as 1D float arrays in the row order of the text
    tables (Tc > ELEV > third parameter > M or D). Only these columns are read from a NetCDF file.
    The grid columns can be requested by their text name: 'Tc', 'ELEV', 'P3' (or 'Fw'/'N'), 'M', 'D'.
    columns=None reads all the columns."""
    path = table_path(path)
    out = {}
    if is_netcdf(path):
        import xarray as xr
        with xr.open_dataset(path) as ds:
            dims = [d for d in ("Tc", "ELEV", "Fw", "N", "D", "M") if d in ds.dims]
            if columns is None:
                columns = ds.attrs["columns"].split()
            grid_alias = {"P3": dims[2]}
            needed_grid = [c for c in columns if grid_alias.get(c, c) in dims]
            if needed_grid:
                grids = np.meshgrid(*[ds[d].values for d in dims], indexing="ij")
                grids = dict(zip(dims, grids))
            for c in columns:
                name = grid_alias.get(c, c)
                if name in dims:
                    out[c] = grids[name].ravel()
                elif name in ds:
                    out[c] = ds[name].transpose(*dims).values.ravel()
                else:
                    raise KeyError(f"Column '{c}' not in {path}")
    else:
        kind, sep, lines, params, all_columns = _read_text_header(path)
        if columns is None:
            columns = all_columns
        df = pd.read_csv(path, sep=r"\s+" if sep is None else sep, skiprows=3, header=None,
                         names=all_columns, usecols=list(dict.fromkeys(columns)), dtype=np.float64,
                         **_TEXT_OPTIONS)
        out = {c: df[c].to_numpy() for c in columns}
    if inf_to_nan:
        for c, v in out.items():
            if v.dtype.kind == "f":
                bad = np.isinf(v)
                if bad.any():
                    if not v.flags.writeable:
                        v = out[c] = v.copy()
                    v[bad] = np.nan
    return out


def read_table_dataframe(path: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Table as a pandas DataFrame with the same columns as the text tables (for the plot tools)."""
    return pd.DataFrame(read_table_columns(path, columns))


def open_table(path: str):
    """Table as an xarray Dataset (NetCDF tables only)."""
    import xarray as xr
    return xr.open_dataset(table_path(path))


# ----------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------
def _main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Convert lookup tables between text and NetCDF.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("to-netcdf", help="text table(s) -> TABLE.nc (deletes the text file unless --keep-text)")
    p1.add_argument("tables", nargs="+")
    p1.add_argument("--keep-text", action="store_true")
    p2 = sub.add_parser("to-text", help="TmatCoefDiff NetCDF -> text file read by TmatInt")
    p2.add_argument("table")
    p2.add_argument("output", nargs="?")
    args = parser.parse_args(argv)
    if args.cmd == "to-netcdf":
        for t in args.tables:
            nc = text_to_netcdf(t)
            if not args.keep_text:
                os.remove(t)
            print(f"{t} -> {nc}")
    else:
        print(netcdf_to_text(args.table, args.output))


if __name__ == "__main__":
    _main(sys.argv[1:])
