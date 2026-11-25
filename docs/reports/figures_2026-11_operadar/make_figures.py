"""Verification of operadar outputs (AROME, ICE3, C band): figures and statistics of the report
docs/reports/2026-11_operadar_verification.md. Run from the repository root:

    python docs/reports/figures_2026-11_operadar/make_figures.py REF=path/to/ref.nc LABEL=path/to/run.nc ...

The first run is the reference; every other run is compared with it. CFADs follow David et al. (2025, AMT,
Fig. 6): median and 25th-75th percentiles inside the convective cores (columns where max_z Zh >= 40 dBZ),
plus frequency CFADs (Yuter and Houze, 1995) normalised at each height.
"""
import os
import sys
import json
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

OUT = os.path.dirname(os.path.abspath(__file__))
VARS = ["Zh", "Zdr", "Kdp", "Ah", "Zh_att"]
UNITS = {"Zh": "dBZ", "Zdr": "dB", "Kdp": "° km⁻¹", "Ah": "dB km⁻¹", "Zh_att": "dBZ", "Rhohv": ""}
XLIM = {"Zh": (10, 60), "Zdr": (-0.5, 4), "Kdp": (-0.25, 3), "Ah": (0, 0.3), "Zh_att": (10, 60)}
BINS = {"Zh": np.arange(-10, 71, 2.0), "Zdr": np.arange(-1, 5.01, 0.1), "Kdp": np.arange(-0.5, 5.01, 0.1),
        "Ah": np.arange(0, 0.5, 0.01), "Zh_att": np.arange(-10, 71, 2.0)}
HYDRO_NAMES = {"cc": "cloud water", "rr": "rain", "ss": "snow", "gg": "graupel", "ii": "ice",
               "wg": "wet graupel", "hh": "hail", "wh": "wet hail"}
HYDRO_COLORS = {"cc": "#17becf", "rr": "#1f77b4", "ss": "#e7298a", "gg": "#7b1c1c", "ii": "#d62728",
                "wg": "#8c6d1f", "hh": "#555555", "wh": "#bcbd22"}
RUN_COLORS = ["#4c4c4c", "#e3711b", "#2a6fdb", "#1a9e77", "#c2366b", "#7a5cc7"]
DZ = 0.5          # height bins (km)
HMAX = 13.
ZCORE = 40.       # convective core threshold on max_z Zh (dBZ)
ZMIN = -30.       # detectability threshold applied to every file (dBZ)


def load(path):
    ds = xr.open_dataset(path)
    alt = ds["Alt"].values.astype(float)
    ground = alt[-1]  # lowest model level (a few metres above the ground)
    out = {"ds": ds, "H": (alt - ground[None]) / 1000., "T": ds["T"].values}
    for v in VARS + ["Rhohv"]:
        if v in ds:
            x = ds[v].values.astype(float)
            x[~np.isfinite(x)] = np.nan
            out[v] = x
    # Same detectability rule for every file (threshold of the C band at 1 km, -30 dBZ), and removal of the
    # Zh = 0 dBZ written by operadar before 2026-11-23 where the linear reflectivity was 0 (no signal)
    z = out["Zh"]
    out["n_zero_dbz"] = int((z == 0).sum())
    out["n_below_zmin"] = int((z < ZMIN).sum())
    hide = ~(z >= ZMIN) | (z == 0)
    for v in ["Zh", "Zdr", "Kdp", "Rhohv"]:
        if v in out:
            out[v] = np.where(hide, np.nan, out[v])
    if "Zh_att" in out:
        out["Zh_att"] = np.where(~(out["Zh_att"] >= ZMIN) | (z == 0), np.nan, out["Zh_att"])
    out["Zmax"] = np.nanmax(np.where(np.isfinite(out["Zh"]), out["Zh"], -999.), axis=0)
    out["core"] = out["Zmax"] >= ZCORE
    return out


def profile(run, var, mask2d, q=(25, 50, 75)):
    """Percentiles of var at each height bin, inside the columns of mask2d."""
    x = run[var]
    H = run["H"]
    m = np.broadcast_to(mask2d[None], x.shape) & np.isfinite(x)
    edges = np.arange(0, HMAX + DZ, DZ)
    res = np.full((len(edges) - 1, len(q)), np.nan)
    n = np.zeros(len(edges) - 1, int)
    for k in range(len(edges) - 1):
        sel = m & (H >= edges[k]) & (H < edges[k + 1])
        vals = x[sel]
        n[k] = vals.size
        if vals.size >= 20:
            res[k] = np.percentile(vals, q)
    return 0.5 * (edges[1:] + edges[:-1]), res, n


def fig_profiles(runs, labels, name, own_mask=True):
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.2), sharey=True)
    for ax, var in zip(axes, ["Zh", "Zdr", "Kdp"]):
        for i, (r, lab) in enumerate(zip(runs, labels)):
            mask = r["core"] if own_mask else runs[0]["core"]
            h, p, _ = profile(r, var, mask)
            c = RUN_COLORS[i % len(RUN_COLORS)]
            ax.fill_betweenx(h, p[:, 0], p[:, 2], color=c, alpha=0.15, lw=0)
            ax.plot(p[:, 1], h, color=c, lw=2, label=f"{lab} Q50")
            ax.plot(p[:, 0], h, color=c, lw=0.6, ls="--")
            ax.plot(p[:, 2], h, color=c, lw=0.6, ls="--")
        ax.set_xlim(*XLIM[var])
        ax.set_xlabel(f"{var} ({UNITS[var]})")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Height above lowest level (km)")
    axes[0].set_ylim(0, HMAX)
    axes[0].legend(fontsize=8, loc="upper right")
    ncols = [int(r["core"].sum()) for r in runs]
    fig.suptitle(f"Convective cores (max_z Zh ≥ {ZCORE:.0f} dBZ): median (line) and Q25-Q75 (shading). "
                 f"Core columns: " + ", ".join(f"{l} {n}" for l, n in zip(labels, ncols)), fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=130)
    plt.close(fig)


def fig_contents(run, name):
    ds = run["ds"]
    hydros = [str(h) for h in ds["hydrometeor"].values]
    M = ds["Contents"].values.astype(float)  # g/m3 (operadar output)
    fig, ax = plt.subplots(figsize=(5, 5.2))
    for k, h in enumerate(hydros):
        x = M[k]
        run[f"M_{h}"] = x
        hh, p, _ = profile(run, f"M_{h}", run["core"])
        ax.fill_betweenx(hh, p[:, 0], p[:, 2], color=HYDRO_COLORS[h], alpha=0.15, lw=0)
        ax.plot(p[:, 1], hh, color=HYDRO_COLORS[h], lw=2, label=HYDRO_NAMES[h])
    ax.set_xlabel("Content (g m⁻³)")
    ax.set_ylabel("Height above lowest level (km)")
    ax.set_ylim(0, HMAX)
    ax.set_xlim(left=0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Hydrometeor contents in the convective cores", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=130)
    plt.close(fig)


def fig_freq_cfad(runs, labels, var, name):
    """Frequency CFADs (% per bin, normalised at each height) in the cores, and the difference to run 0."""
    edges_h = np.arange(0, HMAX + DZ, DZ)
    hists = []
    for r in runs:
        x, H = r[var], r["H"]
        m = np.broadcast_to(r["core"][None], x.shape) & np.isfinite(x)
        h2, _, _ = np.histogram2d(H[m], x[m], bins=[edges_h, BINS[var]])
        tot = h2.sum(axis=1, keepdims=True)
        hists.append(np.where(tot > 0, 100. * h2 / np.maximum(tot, 1), np.nan))
    n = len(runs)
    fig, axes = plt.subplots(1, 2 * n - 1, figsize=(4.2 * (2 * n - 1), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    vmax = np.nanpercentile(np.concatenate([h.ravel() for h in hists]), 99.5)
    for i in range(n):
        pc = axes[i].pcolormesh(BINS[var], edges_h, hists[i], cmap="viridis", vmin=0, vmax=vmax)
        axes[i].set_title(labels[i], fontsize=10)
        axes[i].set_xlabel(f"{var} ({UNITS[var]})")
    fig.colorbar(pc, ax=axes[n - 1], label="% at each height")
    for j in range(1, n):
        d = hists[j] - hists[0]
        lim = np.nanpercentile(np.abs(d), 99) or 1
        pc = axes[n - 1 + j].pcolormesh(BINS[var], edges_h, d, cmap="RdBu_r", vmin=-lim, vmax=lim)
        axes[n - 1 + j].set_title(f"{labels[j]} − {labels[0]}", fontsize=10)
        axes[n - 1 + j].set_xlabel(f"{var} ({UNITS[var]})")
        fig.colorbar(pc, ax=axes[n - 1 + j], label="Δ %")
    for ax in axes:
        ax.set_xlim(*XLIM[var])
    axes[0].set_ylabel("Height above lowest level (km)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=120)
    plt.close(fig)


def fig_scatter(ref, run, lref, lrun, name):
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.3))
    for ax, var in zip(axes, ["Zh", "Zdr", "Kdp", "Ah"]):
        a, b = ref[var], run[var]
        m = np.isfinite(a) & np.isfinite(b)
        bins = BINS[var]
        h, _, _ = np.histogram2d(a[m], b[m], bins=[bins, bins])
        pc = ax.pcolormesh(bins, bins, h.T, norm=LogNorm(vmin=1), cmap="magma_r")
        ax.plot(bins[[0, -1]], bins[[0, -1]], "k--", lw=0.8)
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(bins[0], bins[-1])
        ax.set_xlabel(f"{lref} {var} ({UNITS[var]})")
        ax.set_ylabel(f"{lrun} {var}")
        fig.colorbar(pc, ax=ax, label="points")
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=120)
    plt.close(fig)


def fig_maps(runs, labels, name):
    lon, lat = runs[0]["ds"]["lon"].values, runs[0]["ds"]["lat"].values
    n = len(runs)
    fig, axes = plt.subplots(2, n, figsize=(4.6 * n, 8), squeeze=False)
    for i, (r, lab) in enumerate(zip(runs, labels)):
        z = np.where(r["Zmax"] > -900, r["Zmax"], np.nan)
        pc = axes[0, i].pcolormesh(lon, lat, z, cmap="turbo", vmin=0, vmax=60)
        axes[0, i].contour(lon, lat, r["core"].astype(float), [0.5], colors="k", linewidths=0.6)
        axes[0, i].set_title(f"{lab}: max_z Zh", fontsize=10)
        fig.colorbar(pc, ax=axes[0, i], label="dBZ")
        if i == 0:
            zdr = r["Zdr"][-1]
            pc = axes[1, 0].pcolormesh(lon, lat, zdr, cmap="turbo", vmin=0, vmax=3)
            axes[1, 0].set_title(f"{lab}: Zdr lowest level", fontsize=10)
            fig.colorbar(pc, ax=axes[1, 0], label="dB")
        else:
            z0 = np.where(runs[0]["Zmax"] > -900, runs[0]["Zmax"], np.nan)
            pc = axes[1, i].pcolormesh(lon, lat, z - z0, cmap="RdBu_r", vmin=-5, vmax=5)
            axes[1, i].set_title(f"{lab} − {labels[0]}: max_z Zh", fontsize=10)
            fig.colorbar(pc, ax=axes[1, i], label="dB")
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=110)
    plt.close(fig)


def fig_section(runs, labels, name):
    ref = runs[0]
    j, i = np.unravel_index(np.argmax(ref["Zmax"]), ref["Zmax"].shape)
    x = ref["ds"]["lon"].values[j]
    n = len(runs)
    fig, axes = plt.subplots(3, n, figsize=(4.6 * n, 9), squeeze=False, sharex=True, sharey=True)
    for c, (r, lab) in enumerate(zip(runs, labels)):
        H = r["H"][:, j, :]
        X = np.broadcast_to(x[None], H.shape)
        for k, (var, lim, cm) in enumerate([("Zh", (0, 60), "turbo"), ("Zdr", (-0.5, 4), "turbo"),
                                            ("Kdp", (0, 3), "turbo")]):
            pc = axes[k, c].pcolormesh(X, H, r[var][:, j, :], cmap=cm, vmin=lim[0], vmax=lim[1], shading="nearest")
            axes[k, c].contour(X, H, r["T"][:, j, :], [0], colors="w", linewidths=0.8)
            axes[k, c].set_title(f"{lab}: {var}", fontsize=10)
            fig.colorbar(pc, ax=axes[k, c], label=UNITS[var])
    for ax in axes[:, 0]:
        ax.set_ylabel("Height (km)")
        ax.set_ylim(0, HMAX)
    for ax in axes[-1]:
        ax.set_xlabel("Longitude (°)")
    fig.suptitle(f"West-east section through the reflectivity maximum (lat {ref['ds']['lat'].values[j, i]:.2f}°), "
                 "white line: 0 °C", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=110)
    plt.close(fig)


def regions(ref):
    ds = ref["ds"]
    hydros = [str(h) for h in ds["hydrometeor"].values]
    M = ds["Contents"].values
    wet = np.zeros(ref["T"].shape, bool)
    for h in ("wg", "wh"):
        if h in hydros:
            wet |= M[hydros.index(h)] > 0
    T = ref["T"]
    return {"all": np.ones(T.shape, bool), "liquid (T>0, no wet ice)": (T > 0) & ~wet,
            "melting layer (wet ice)": wet, "ice (T<0, no wet ice)": (T <= 0) & ~wet}


def stats(ref, run):
    res = {}
    regs = regions(ref)
    for var in VARS:
        if var not in ref or var not in run:
            continue
        a, b = ref[var], run[var]
        res[var] = {"nan_ref": int(np.isnan(a).sum()), "nan_run": int(np.isnan(b).sum()),
                    "min_run": float(np.nanmin(b)), "max_run": float(np.nanmax(b))}
        for rn, rm in regs.items():
            # differences where both are defined (detectable echoes, see load)
            m = rm & np.isfinite(a) & np.isfinite(b)
            if m.sum() == 0:
                continue
            d = b[m] - a[m]
            res[var][rn] = {"n": int(m.sum()), "bias": float(d.mean()), "mad": float(np.abs(d).mean()),
                            "p01": float(np.percentile(d, 1)), "p99": float(np.percentile(d, 99)),
                            "maxabs": float(np.abs(d).max())}
    return res


def main(argv):
    pairs = [a.split("=", 1) for a in argv]
    labels = [p[0] for p in pairs]
    runs = [load(p[1]) for p in pairs]
    tag = os.environ.get("FIGTAG", "")
    fig_profiles(runs, labels, f"profiles_cores{tag}.png")
    fig_profiles(runs, labels, f"profiles_cores_refmask{tag}.png", own_mask=False)
    for var in ["Zh", "Zdr", "Kdp"]:
        fig_freq_cfad(runs, labels, var, f"cfad_{var}{tag}.png")
    fig_maps(runs, labels, f"maps_Zmax{tag}.png")
    fig_section(runs, labels, f"section{tag}.png")
    fig_contents(runs[0], f"contents_cores{tag}.png")
    allstats = {}
    for lab, r in zip(labels[1:], runs[1:]):
        fig_scatter(runs[0], r, labels[0], lab, f"scatter_{lab}{tag}.png")
        allstats[lab] = stats(runs[0], r)
    # sample sizes and median profiles for the report
    allstats["_cores"] = {lab: int(r["core"].sum()) for lab, r in zip(labels, runs)}
    allstats["_zero_dbz"] = {lab: r["n_zero_dbz"] for lab, r in zip(labels, runs)}
    allstats["_below_zmin"] = {lab: r["n_below_zmin"] for lab, r in zip(labels, runs)}
    allstats["_valid_zh"] = {lab: int(np.isfinite(r["Zh"]).sum()) for lab, r in zip(labels, runs)}
    allstats["_profiles"] = {}
    for lab, r in zip(labels, runs):
        allstats["_profiles"][lab] = {}
        for var in ["Zh", "Zdr", "Kdp"]:
            h, p, n = profile(r, var, r["core"])
            allstats["_profiles"][lab][var] = {f"{hh:.2f}": [None if np.isnan(v) else round(float(v), 3) for v in pp]
                                               for hh, pp in zip(h, p)}
    with open(f"{OUT}/stats{tag}.json", "w") as f:
        json.dump(allstats, f, indent=1)
    print(json.dumps({k: v for k, v in allstats.items() if k != "_profiles"}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
