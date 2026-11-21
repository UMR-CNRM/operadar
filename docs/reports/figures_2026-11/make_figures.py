"""Figures 6-8 of the report 2026-11_lookup_tables_update.md (run from the repository root).

    python docs/reports/figures_2026-11/make_figures.py

Figures 1-5 come from plot_tools (plot_tables_2D.py, plot_D_M_varpol_alltypes.py).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAB = "tables_generator/tables"
OUT = os.path.dirname(os.path.abspath(__file__))
COLS = ["Tc", "ELEV", "P3", "M", "zhh", "zdr", "kdp", "rhohv", "zhhR", "zdrR", "kdpR", "Mfrac"]
COLORS = ["#2a6fdb", "#e3711b", "#1a9e77", "#c2366b", "#7a5cc7", "#8c6d1f", "#4c4c4c", "#d4a017"]


def nearest(values, x):
    values = np.unique(values)
    return values[np.argmin(np.abs(values - x))]


def load(config, name, Tc, P3=None, elev=None):
    """Rows of TmatCoefInt_<name> at the grid values closest to Tc and P3 (and elev)."""
    path = f"{TAB}/{config}/TmatCoefInt_{name}"
    keep = []
    for chunk in pd.read_csv(path, sep=";", skiprows=2, usecols=COLS, chunksize=2_000_000,
                             dtype={"rhohv": str}):
        keep.append(chunk[np.abs(chunk.Tc - Tc) < 1e-3])
    df = pd.concat(keep)
    # rhohv is NaN/Infinity where Zh < -150 dBZ (see the wiki page): read it as a number
    df["rhohv"] = pd.to_numeric(df.rhohv, errors="coerce")
    if df.empty:
        raise ValueError(f"{path}: Tc={Tc} not on the grid")
    if P3 is not None:
        df = df[df.P3 == nearest(df.P3.values, P3)]
    if elev is not None:
        df = df[np.abs(df.ELEV - elev) < 1e-3]
    return df


def at_M(df, M):
    m = nearest(df.M.values, M)
    return df[df.M == m].sort_values("ELEV")


def fig_elevation():
    """Radar variables as a function of the elevation, David2026PhD, C band."""
    cases = [
        ("rain, 10 °C, 1 g/m³", "ICE3_1M_Crr", 10, None, 1e-3),
        ("rain, 10 °C, 5 g/m³", "ICE3_1M_Crr", 10, None, 5e-3),
        ("wet graupel, 0 °C, Fw 0.5, 5 g/m³", "ICE3_1M_Cwg", 0, 0.5, 5e-3),
        ("snow, -10 °C, 1 g/m³", "ICE3_1M_Css", -10, None, 1e-3),
        ("pristine ice, -30 °C, 0.1 g/m³, N 1000 m⁻³", "LIMA_2M_Cii", -30, 1000, 1e-4),
    ]
    variables = [("zhh", "Zh (dBZ)"), ("zdr", "Zdr (dB)"), ("kdp", "Kdp (°/km)"), ("rhohv", "ρhv")]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.8))
    for i, (label, name, Tc, P3, M) in enumerate(cases):
        d = at_M(load("David2026PhD", name, Tc, P3), M)
        for ax, (v, vlabel) in zip(axes, variables):
            ax.plot(d.ELEV, d[v], "-o", color=COLORS[i], lw=2, ms=4, label=label)
            if name.endswith("rr") and v in ("zdr", "kdp"):
                ax.plot(d.ELEV, d[v + "R"], "--", color=COLORS[i], lw=1.2,
                        label=label + " (Rayleigh)" if v == "zdr" else None)
    for ax, (v, vlabel) in zip(axes, variables):
        ax.set_xlabel("elevation (°)")
        ax.set_title(vlabel)
        ax.set_xticks(range(0, 91, 15))
        ax.grid(alpha=0.3)
    axes[3].set_ylim(0.95, 1.001)
    axes[0].legend(fontsize=8, loc="lower left")
    axes[1].legend(*axes[1].get_legend_handles_labels(), fontsize=7, loc="upper right")
    fig.suptitle("David2026PhD, C band (canting orient): dependence on the elevation", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/elevation_C_David2026PhD.png", dpi=110)
    plt.close(fig)


def fig_configurations():
    """David2025AMT (solid) vs David2026PhD (dashed), C band, 0 deg elevation."""
    cases = [
        ("rain 10 °C", "ICE3_1M_Crr", 10, None),
        ("pristine ice -30 °C, N 1000 m⁻³", "LIMA_2M_Cii", -30, 1000),
        ("snow -10 °C", "ICE3_1M_Css", -10, None),
        ("graupel 0 °C", "ICE3_1M_Cgg", 0, None),
        ("wet graupel 0 °C, Fw 0.5", "ICE3_1M_Cwg", 0, 0.5),
    ]
    variables = [("zhh", "Zh (dBZ)"), ("zdr", "Zdr (dB)"), ("kdp", "Kdp (°/km)"), ("rhohv", "ρhv")]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.8))
    for i, (label, name, Tc, P3) in enumerate(cases):
        for config, ls in (("David2025AMT", "-"), ("David2026PhD", "--")):
            d = load(config, name, Tc, P3, elev=0).sort_values("M")
            d = d[(d.Mfrac > 0.99) & (d.M <= 1e-2)]
            for ax, (v, vlabel) in zip(axes, variables):
                ax.plot(d.M * 1e3, d[v], ls, color=COLORS[i], lw=2,
                        label=f"{label}" if config == "David2025AMT" else None)
    for ax, (v, vlabel) in zip(axes, variables):
        ax.set_xscale("log")
        ax.set_xlabel("content M (g/m³)")
        ax.set_title(vlabel)
        ax.grid(alpha=0.3)
    axes[2].set_yscale("symlog", linthresh=0.01)
    axes[3].set_ylim(0.9, 1.002)
    axes[0].plot([], [], "k-", label="David2025AMT")
    axes[0].plot([], [], "k--", label="David2026PhD")
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("C band, 0° elevation: David2025AMT (solid) vs David2026PhD (dashed)", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/configurations_C.png", dpi=110)
    plt.close(fig)


def fig_bands():
    """David2025AMT rain and snow at the 8 bands, 0 deg elevation."""
    bands = ["L", "S", "C", "X", "Ku", "K", "Ka", "W"]
    cases = [("rain 10 °C", "rr", 10), ("snow -10 °C", "ss", -10)]
    variables = [("zhh", "Zh (dBZ)"), ("zdr", "Zdr (dB)"), ("kdp", "Kdp (°/km)")]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for row, (label, t, Tc) in enumerate(cases):
        for i, b in enumerate(bands):
            d = load("David2025AMT", f"ICE3_1M_{b}{t}", Tc, elev=0).sort_values("M")
            d = d[(d.Mfrac > 0.99) & (d.M <= 1e-2) & (d.M >= 1e-5)]
            for ax, (v, vlabel) in zip(axes[row], variables):
                ax.plot(d.M * 1e3, d[v], "-", color=COLORS[i], lw=2, label=b)
                if v == "zhh":
                    ax.plot(d.M * 1e3, d[v + "R"], ":", color=COLORS[i], lw=1)
        for ax, (v, vlabel) in zip(axes[row], variables):
            ax.set_xscale("log")
            ax.set_xlabel("content M (g/m³)")
            ax.set_title(f"{label}: {vlabel}")
            ax.grid(alpha=0.3)
    axes[0][0].legend(ncol=2, fontsize=9, title="band (dotted: Rayleigh)")
    fig.suptitle("David2025AMT, ICE3, 0° elevation: all bands", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/bands_David2025AMT.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    fig_elevation()
    fig_configurations()
    fig_bands()
    print("figures written to", OUT)
