# -*- coding: utf-8 -*-
"""
2D maps of the lookup tables, covering their full content (all temperatures, contents, liquid water
fractions and concentrations), as a complement to the curves of plot_D_M_varpol_*.py.

For each hydrometeor type, the figures are:
  * Int  M-Tc : variables as a function of the content M (x, log) and the temperature (y),
                at a fixed P3 (Fw for wet species, N for 2-moment species and pristine ice);
  * Int  M-P3 : (wet species, 2-moment species, ice) variables as a function of M and of Fw or N,
                at a fixed temperature;
  * Diff D-Tc : single-particle variables as a function of the diameter D (x, log) and the temperature
                (dry species), or D-Fw at a fixed temperature (wet species).
The Int maps also show Mfrac, the fraction of the content inside the diameter range of the table:
where Mfrac < 0.99 (hatched) the size distribution is cut by the table and the values are too low.
"Zh-ZhR" is the difference between the T-matrix and the Rayleigh reflectivities (Mie effects).

Examples (from the plot_tools folder):
    python plot_tables_2D.py --band C
    python plot_tables_2D.py --band C --micro LIMA --types rr cl ii
    python plot_tables_2D.py --band W --config David2025AMT --types ss --vars Zh Zdr --Tc -20
Figures are written to IMG/tables2D/<config>/<band>/.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, "..", "tables_generator", "tables")

# Number of moments of each type in each microphysics scheme (as in generate_tables.sh / operadar)
MOMENTS = {
    "ICE3": {"rr": 1, "ss": 1, "gg": 1, "wg": 1, "cl": 1, "cs": 1, "ii": 1, "hh": 1, "wh": 1},
    "ICJW": {"rr": 1, "ss": 1, "gg": 1, "wg": 1, "cl": 1, "cs": 1, "ii": 1, "hh": 1, "wh": 1},
    "LIMA": {"rr": 2, "ss": 1, "gg": 1, "wg": 1, "cl": 2, "cs": 2, "ii": 2, "hh": 1, "wh": 1},
    "LIMC": {"rr": 2, "ss": 1, "gg": 1, "wg": 1, "cl": 2, "cs": 2, "ii": 1, "hh": 1, "wh": 1},
}
TYPE_NAME = {"rr": "rain", "cl": "cloud water (land)", "cs": "cloud water (sea)", "ii": "pristine ice",
             "ss": "snow", "gg": "graupel", "hh": "hail", "wg": "wet graupel", "wh": "wet hail"}
WET = ("wg", "wh")

# Default slice values: temperature of the M-P3 / D-Fw maps, Fw of the M-Tc maps of wet species,
# N (/m3) of the M-Tc maps of 2-moment species and ice (the closest grid value is used)
TC_SLICE = {"rr": 10, "cl": 10, "cs": 10, "ii": -30, "ss": -10, "gg": -10, "hh": -10, "wg": 0, "wh": 0}
FW_SLICE = 0.2
N_SLICE = {"rr": 1e3, "cl": 1e8, "cs": 1e8, "ii": 1e4}

# Variables: table column(s), label, unit, colour scale
VARS = {
    "Zh":     dict(col="zhh", label="Zh", unit="dBZ"),
    "Zdr":    dict(col="zdr", label="Zdr", unit="dB"),
    "Kdp":    dict(col="kdp", label="Kdp", unit="° km$^{-1}$"),
    "Rhohv":  dict(col="rhohv", label=r"$\rho_{hv}$", unit=""),
    "Ah":     dict(col="Ah", label="Ah", unit="dB km$^{-1}$"),
    "Zh-ZhR": dict(col=("zhh", "zhhR"), label="Zh - Zh Rayleigh", unit="dB"),
}
# A variable whose range is below this tolerance (in the cells covered by the table) is shown as a
# constant (e.g. Zdr = 0 and rhohv = 1 for spheres, Kdp = 0 for spheres)
CONST_TOL = {"Zh": 1e-3, "Zdr": 1e-3, "Kdp": 1e-6, "Rhohv": 1e-5, "Ah": 1e-9, "Zh-ZhR": 1e-3}
# Diameter axis of the single-particle maps: the grid diameter D of the table (equivolume diameter of the
# particle before melting; for wet species the column Deq, after melting, changes with Fw)
DLABEL = "D (mm, equivolume, before melting)"
DEFAULT_VARS = ["Zh", "Zdr", "Kdp", "Rhohv", "Ah", "Zh-ZhR"]


def read_table(path, kind, elev):
    """Header parameters and data (rows at the elevation closest to elev) of a Diff or Int table."""
    sep = r"\s+" if kind == "Diff" else ";"
    header = pd.read_csv(path, sep=sep, nrows=1, engine="python")
    df = pd.read_csv(path, sep=sep, skiprows=2, engine="c" if kind == "Int" else "python")
    df.columns = [c.strip() for c in df.columns]
    elevs = np.unique(df["ELEV"])
    e = elevs[np.argmin(np.abs(elevs - elev))]
    return header.iloc[0], df[df["ELEV"] == e], e


def closest(values, target, log=False):
    grid = np.unique(values)
    d = np.abs(np.log10(np.maximum(grid, 1e-30)) - np.log10(target)) if log else np.abs(grid - target)
    return grid[np.argmin(d)]


def values_of(df, var):
    col = VARS[var]["col"]
    if isinstance(col, tuple):
        return df[col[0]].to_numpy(float) - df[col[1]].to_numpy(float)
    return df[col].to_numpy(float)


def grid2d(df, xcol, ycol, z):
    """(x, y, Z[y, x]) from the table rows (x and y take all the grid values)."""
    x, xi = np.unique(df[xcol].to_numpy(float), return_inverse=True)
    y, yi = np.unique(df[ycol].to_numpy(float), return_inverse=True)
    Z = np.full((len(y), len(x)), np.nan)
    Z[yi, xi] = z
    return x, y, Z


def norm_and_cmap(var, Z):
    """Colour scale of one variable (robust to outliers, centred on 0 where the sign matters)."""
    finite = Z[np.isfinite(Z)]
    if finite.size == 0:
        return mcolors.Normalize(0, 1), "viridis"
    if var == "Zh":
        return mcolors.Normalize(max(np.nanmin(finite), -40), np.nanmax(finite)), "viridis"
    if var == "Rhohv":
        return mcolors.Normalize(max(np.percentile(finite, 1), 0.8), 1.0), "magma"
    if var in ("Zdr", "Zh-ZhR"):
        v = max(np.percentile(np.abs(finite), 99.5), 1e-3)
        return mcolors.TwoSlopeNorm(0, -v, v), "RdBu_r"
    if var == "Kdp":  # spans many decades and can be negative (hail, vertical incidence)
        v = max(np.nanmax(np.abs(finite)), 1e-4)
        return mcolors.SymLogNorm(linthresh=max(v * 1e-4, 1e-5), vmin=-v, vmax=v), "RdBu_r"
    if var == "Ah":
        pos = finite[finite > 0]
        vmax = pos.max() if pos.size else 1.0
        return mcolors.LogNorm(max(vmax * 1e-6, pos.min() if pos.size else 1e-9), vmax), "viridis"
    return mcolors.Normalize(np.nanmin(finite), np.nanmax(finite)), "viridis"


def panel_axes(ax, lab, xlabel, ylabel, logx, logy):
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(lab.split(" (")[0], fontsize=11)


def draw_maps(df, xcol, ycol, varlist, title, xlabel, ylabel, outfile, logx=True, logy=False, mfrac=True):
    panels = list(varlist) + (["Mfrac"] if mfrac and "Mfrac" in df.columns else [])
    ncols = 3
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.0 * nrows), squeeze=False,
                             layout="constrained")
    cut = None
    if mfrac and "Mfrac" in df.columns:
        _, _, cut = grid2d(df, xcol, ycol, df["Mfrac"].to_numpy(float))
    for ax, var in zip(axes.flat, panels):
        if var == "Mfrac":
            x, y, Z = grid2d(df, xcol, ycol, df["Mfrac"].to_numpy(float))
            norm, cmap, lab = mcolors.Normalize(0, 1), "cividis", "Mfrac (fraction of M inside the table)"
        else:
            x, y, Z = grid2d(df, xcol, ycol, values_of(df, var))
            # colour scale from the cells covered by the table only: where the size distribution is cut
            # (Mfrac < 0.99, hatched) the values can be meaningless (e.g. Zdr of -700 dB for N > 1e7 /m3 ice)
            norm, cmap = norm_and_cmap(var, Z if cut is None else np.where(cut >= 0.99, Z, np.nan))
            u = VARS[var]["unit"]
            lab = VARS[var]["label"] + (f" ({u})" if u else "")
        covered = Z if cut is None or var == "Mfrac" else np.where(cut >= 0.99, Z, np.nan)
        finite = covered[np.isfinite(covered)]
        if finite.size and np.ptp(finite) <= CONST_TOL.get(var, 1e-6) * max(1.0, np.abs(finite).max()):
            # constant field (e.g. rhohv = 1 for a single particle without canting): no colour scale
            ax.pcolormesh(x, y, np.zeros_like(Z), cmap="Greys", vmin=0, vmax=4, shading="nearest", rasterized=True)
            where = "" if cut is None or not np.any(cut < 0.99) else "\n(outside the hatched area)"
            ax.text(0.5, 0.5, f"constant: {np.median(finite):.3g}{where}", transform=ax.transAxes,
                    ha="center", va="center", bbox=dict(fc="white", ec="none", alpha=0.8))
            if cut is not None and np.any(cut < 0.99):
                ax.contourf(x, y, cut, levels=[-1, 0.99], colors="none", hatches=["////"])
                ax.contour(x, y, cut, levels=[0.99], colors="k", linewidths=0.8)
            panel_axes(ax, lab, xlabel, ylabel, logx, logy)
            continue
        pm = ax.pcolormesh(x, y, Z, norm=norm, cmap=cmap, shading="nearest", rasterized=True)
        cb = fig.colorbar(pm, ax=ax, pad=0.01, extend="both" if cut is not None and np.any(cut < 0.99) else "neither")
        cb.set_label(lab, fontsize=10)
        if isinstance(norm, mcolors.SymLogNorm):  # decades only, away from the linear zone around 0
            v, lt = norm.vmax, norm.linthresh
            dec = [10.0 ** k for k in range(int(np.ceil(np.log10(lt))) + 1, int(np.floor(np.log10(v))) + 1)][-4:]
            cb.set_ticks([-d for d in dec[::-1]] + [0] + dec)
        if var == "Zh":  # a few iso-reflectivity lines help reading the values
            levels = [l for l in (0, 10, 20, 30, 40, 50, 60) if np.nanmin(Z) < l < np.nanmax(Z)]
            if levels:
                cs = ax.contour(x, y, Z, levels=levels, colors="white", linewidths=0.6, alpha=0.8)
                ax.clabel(cs, fontsize=7, fmt="%d", inline_spacing=2, use_clabeltext=True)
        if cut is not None and np.any(cut < 0.99):
            ax.contourf(x, y, cut, levels=[-1, 0.99], colors="none", hatches=["////"])
            ax.contour(x, y, cut, levels=[0.99], colors="k", linewidths=0.8)
        panel_axes(ax, lab, xlabel, ylabel, logx, logy)
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    note = "  -  hatched: Mfrac < 0.99 (size distribution cut by the table)" if cut is not None and np.any(cut < 0.99) else ""
    fig.suptitle(title + note, fontsize=13)
    fig.savefig(outfile, dpi=130)
    plt.close(fig)
    print("  saved", outfile)


def plot_type(h, band, config, micro, varlist, elev, outdir, tc_user, fw_user, n_user):
    folder = os.path.join(TABLES, config)
    nmom = MOMENTS[micro][h]
    tc0 = TC_SLICE[h] if tc_user is None else tc_user
    name = TYPE_NAME[h]

    # --- Integrated table (as a function of the content M)
    path = os.path.join(folder, f"TmatCoefInt_{micro}_{nmom}M_{band}{h}")
    if os.path.exists(path):
        hd, df, e = read_table(path, "Int", elev)
        base = f"{band} band, {name}, {micro} {nmom}M ({config}), elev {e:g}°"
        if h in WET:
            p3, p3name, logp3 = closest(df["P3"], fw_user if fw_user is not None else FW_SLICE), "Fw", False
        elif nmom == 2 or h == "ii":
            target = n_user if n_user is not None else N_SLICE.get(h, 1e3)
            p3, p3name, logp3 = closest(df["P3"], target, log=True), "N", True
        else:
            p3, p3name, logp3 = None, None, False
        sub = df if p3 is None else df[df["P3"] == p3]
        sl = "" if p3 is None else f", {p3name} = {p3:g}" + (" m$^{-3}$" if p3name == "N" else "")
        draw_maps(sub, "M", "Tc", varlist, f"{base}{sl}", "M (kg m$^{-3}$)", "Tc (°C)",
                  os.path.join(outdir, f"Int_{micro}_{nmom}M_{band}{h}_M-Tc.png"))
        if p3 is not None:
            t = closest(df["Tc"], tc0)
            sub = df[df["Tc"] == t]
            ylabel = "Fw" if p3name == "Fw" else "N (m$^{-3}$)"
            draw_maps(sub, "M", "P3", varlist, f"{base}, Tc = {t:g} °C", "M (kg m$^{-3}$)", ylabel,
                      os.path.join(outdir, f"Int_{micro}_{nmom}M_{band}{h}_M-{p3name}.png"), logy=logp3)
    else:
        print("  missing", path)

    # --- Single-particle table (as a function of the diameter D)
    path = os.path.join(folder, f"TmatCoefDiff_{band}{h}")
    if os.path.exists(path):
        hd, df, e = read_table(path, "Diff", elev)
        base = f"{band} band, {name}, single particle ({config}), elev {e:g}°"
        if h in WET:
            t = closest(df["Tc"], tc0)
            draw_maps(df[df["Tc"] == t], "D", "Fw", varlist, f"{base}, Tc = {t:g} °C",
                      DLABEL, "Fw", os.path.join(outdir, f"Diff_{band}{h}_D-Fw.png"), mfrac=False)
            fw = closest(df["Fw"], fw_user if fw_user is not None else FW_SLICE)
            draw_maps(df[df["Fw"] == fw], "D", "Tc", varlist, f"{base}, Fw = {fw:g}",
                      DLABEL, "Tc (°C)", os.path.join(outdir, f"Diff_{band}{h}_D-Tc.png"), mfrac=False)
        else:
            draw_maps(df, "D", "Tc", varlist, base, DLABEL, "Tc (°C)",
                      os.path.join(outdir, f"Diff_{band}{h}_D-Tc.png"), mfrac=False)
    else:
        print("  missing", path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--band", default="C", help="radar band (L, S, C, X, Ku, K, Ka, W). Default: C")
    parser.add_argument("--config", default="David2026PhD", help="table folder: David2026PhD (default; L, S, C), David2025AMT (all bands), vertical...")
    parser.add_argument("--micro", default="ICE3", choices=list(MOMENTS), help="microphysics scheme. Default: ICE3")
    parser.add_argument("--types", nargs="*", default=["rr", "cl", "cs", "ii", "ss", "gg", "hh", "wg", "wh"],
                        help="hydrometeor types (missing tables are skipped)")
    parser.add_argument("--vars", nargs="*", default=DEFAULT_VARS, choices=list(VARS), help="variables to map")
    parser.add_argument("--elev", type=float, default=0.0, help="elevation (closest table value). Default: 0")
    parser.add_argument("--Tc", type=float, default=None, help="temperature of the M-Fw, M-N and D-Fw maps")
    parser.add_argument("--Fw", type=float, default=None, help=f"Fw of the M-Tc maps of wet species (default {FW_SLICE})")
    parser.add_argument("--N", type=float, default=None, help="N (/m3) of the M-Tc maps of 2-moment species and ice")
    parser.add_argument("--out", default=os.path.join(HERE, "IMG", "tables2D"), help="output folder")
    a = parser.parse_args()

    outdir = os.path.join(a.out, a.config, a.band)
    os.makedirs(outdir, exist_ok=True)
    for h in a.types:
        print(f"{a.band} {h} ({a.micro})")
        plot_type(h, a.band, a.config, a.micro, a.vars, a.elev, outdir, a.Tc, a.Fw, a.N)
