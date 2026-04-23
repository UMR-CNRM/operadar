#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:59:47 2026

@author: augrosc, borderiesm

Radar attenuation module.
=========================
Provides:
  - compute_gaz_attenuation       : gas attenuation via pyMPM
  - compute_relative_humidity     : RH from T, rv, P
  - compute_extinction            : total (gas + hydrometeor) extinction
  - compute_attenuated_zh_3D      : vectorized 3-D two-way attenuation
  - W-band empirical relationships (94 GHz)
  - Visualization utilities (activated only under __main__)
"""

# ======================================================================
# Standard imports
# ======================================================================
import numpy as np

# Visualization imports are optional: only loaded when needed
# so that the module can be imported in non-graphical environments
# (e.g. HPC batch jobs) without requiring matplotlib.
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as mticker
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    _MATPLOTLIB_AVAILABLE = True
except ImportError:
    _MATPLOTLIB_AVAILABLE = False

# pyMPM is now included in operadar
from operadar.radar.pyMPM.MPM import MPM


# ======================================================================
# Plot styling constants
# (defined at module level so every plot function can use them;
#  they are plain Python scalars/strings → safe even without matplotlib)
# ======================================================================

# -- Font sizes --
FONTSIZE_TITLE  = 15
FONTSIZE_LABEL  = 13
FONTSIZE_TICK   = 11
FONTSIZE_LEGEND = 10
FONTSIZE_CBAR   = 12

# -- Overlay colours --
RADAR_COLOR = "deepskyblue"
CLOUD_COLOR = "cyan"
STRAT_COLOR = "#00FF7F"    # spring green
CONV_COLOR  = "#FF69B4"    # hot pink
BB_COLOR    = "lightyellow"

# -- Background colours (dark theme) --
BG_DARK    = "#1a1a2e"    # figure background
BG_PANEL   = "#16213e"    # axes background
BG_PROFILE = "#0f0f23"    # profile panel background


# ======================================================================
# Gas attenuation
# ======================================================================

def compute_gaz_attenuation(temperatureC, pressure, qv, LAM):
    """
    Computes gaz attenuation using pyMPM

    Parameters
    ----------
    temperatureC : 3D array
        Temperature [C]
    pressure : 3D array
        Pressure [Pa]
    qv : 3D array
        Specific humidity [kg/kg]
    LAM: float
        radar wavelength (meters)

    Returns
    -------
    kext_gaz : 3D array
        Gas extinction coefficient (m⁻¹), shape identical to temperature.
    """
    # Convert specific humidity -> mixing ratio rv
    rv = qv / (1 - qv)

    # Compute relative humidity
    rh = compute_relative_humidity(temperatureC, rv, pressure)

    # Frequency in Hz (c / λ)
    freq = 3e8 / LAM  # Hz

    # Call pyMPM (requires frequency in GHz, pressure in hPa, temperature in °C)
    # MPM returns shape (Nfreq, ...) even for a single frequency → squeeze axis 0
    kext_gaz_raw = (
        1e-3 / 4.343
        * MPM(
            freq * 1e-9,
            pressure * 1e-2,
            temperatureC,
            rh,
            wa='None',
            wae='None',
            R='None',
            output_type='att'
        )
    )

    # MPM prepends a frequency axis when called with a scalar frequency;
    # remove it so kext_gaz has the same shape as temperature.
    kext_gaz = np.squeeze(kext_gaz_raw, axis=0) if kext_gaz_raw.ndim == temperatureC.ndim + 1 else kext_gaz_raw

    return kext_gaz


def compute_relative_humidity(temperatureC, rv, pressure_Pa):
    """
    Computes relative humidity (%) from:
        - temperature in Celcius
        - water vapor mixing ratio rv (kg/kg)
        - pressure in Pascals

    Formula:
        Psat : Saturation vapor pressure using Tetens formula (hPa)
        Pe   : Partial vapor pressure (hPa)
        RH   : 100 * Pe / Psat

    Parameters
    ----------
    temperatureC : array-like or float
        Air temperature [C]
    rv : array-like or float
        Water vapor mixing ratio [kg/kg]
    pressure_Pa : array-like or float
        Atmospheric pressure [Pa]

    Returns
    -------
    rh : array-like or float
        Relative humidity in %
    """

    # Saturation vapor pressure (hPa)
    Psat = 6.107 * np.power(10., (7.5 * temperatureC) / (237.3 + temperatureC))

    # Partial vapor pressure (hPa)
    Pe = (rv * pressure_Pa / (rv + 0.622)) / 100.0

    # Relative humidity in %
    rh = 100.0 * (Pe / Psat)

    return rh


def compute_extinction(temperatureC: np.ndarray,
                       pressure: np.ndarray,
                       qv: np.ndarray,
                       radar_lam: float,
                       ah: np.ndarray,
                       ) -> np.ndarray:
    """Computes gaz and hydrometeor extinction cross-section

    Args:
        temperature (np.ndarray): 3D array (deg Celcius)
        pressure (np.ndarray): 3D array (Pa)
        qv (np.ndarray): 3D array specific humidity (kg/kg)
        radar_lam: radar wavelength (m)
        ah : 3D array (dB/km)

    Returns:
        kext (np.ndarray): 3D array of total extinction coefficient (m⁻¹)
    """
    # Gas extinction cross-section
    kext_gaz = compute_gaz_attenuation(temperatureC, pressure, qv, radar_lam)

    # Hydrometeor extinction cross-section
    kext_hydro = 1e-3 * ah / 4.343

    # Defensive check: squeeze any spurious leading dimensions
    # (e.g. frequency axis returned by MPM for a single frequency)
    target_shape = ah.shape
    if kext_gaz.shape != target_shape:
        kext_gaz = kext_gaz.reshape(target_shape)

    # Total extinction cross-section
    kext_total = kext_hydro + kext_gaz

    # Final shape check before returning
    if kext_total.shape != target_shape:
        raise ValueError(
            f"compute_extinction: shape mismatch after combining terms. "
            f"kext_hydro={kext_hydro.shape}, kext_gaz={kext_gaz.shape}, "
            f"expected={target_shape}"
        )

    return kext_total


# ======================================================================
# Core attenuation function
# ======================================================================

def compute_attenuated_zh_3D(
        kext:np.ndarray,
        model_altitude:np.ndarray,
        zh:np.ndarray,
        radar_altitude:float,
        ) -> np.ndarray:
    """
    Compute attenuated radar reflectivity (dBZ) using vertical integration.
    Optimized for speed using vectorized cumulative sums.

    Parameters
    ----------
    kext : ndarray (Nz, Ny, Nx)
        Extinction coefficient (m⁻¹).
    model_altitude : ndarray (Nz, Ny, Nx)
        Altitude of each grid point (m).
    zh : ndarray (Nz, Ny, Nx)
        Unattenuated reflectivity (dBZ).
    radar_altitude : float
        Radar altitude (m).

    Returns
    -------
    zh_att : ndarray (Nz, Ny, Nx)
        Attenuated reflectivity (dBZ), in the same vertical order as the input.
    """
    # ------------------------------------------------------------------
    # 1. Setup and type conversion
    # ------------------------------------------------------------------
    kext = np.asarray(kext, dtype=float)
    zh   = np.asarray(zh,   dtype=float)
    model_altitude  = np.asarray(model_altitude,  dtype=float)

    if not (kext.shape == zh.shape == model_altitude.shape):
        raise ValueError(
            f"Shape mismatch: kext={kext.shape}, zh={zh.shape}, model_altitude={model_altitude.shape}"
        )

    Nz, Ny, Nx = zh.shape

    # ------------------------------------------------------------------
    # 2. Convert zh to linear units
    # ------------------------------------------------------------------
    zh_lin = 10.0 ** (zh / 10.0)

    # ------------------------------------------------------------------
    # 3. Sort every column so altitude increases with index k
    # ------------------------------------------------------------------
    sort_idx = np.argsort(model_altitude, axis=0)          # (Nz, Ny, Nx)

    # Build broadcast-safe helper indices
    y_idx = np.arange(Ny)[None, :, None]        # (1, Ny, 1)
    x_idx = np.arange(Nx)[None, None, :]        # (1, 1, Nx)

    # Apply the sort to every 3-D field
    model_altitude_sorted    = model_altitude[sort_idx,    y_idx, x_idx]   # (Nz, Ny, Nx)
    kext_sorted   = kext[sort_idx,   y_idx, x_idx]
    zh_lin_sorted = zh_lin[sort_idx, y_idx, x_idx]

    # ------------------------------------------------------------------
    # 4. Vectorized integration of extinction (optical depth)
    # ------------------------------------------------------------------
    dz        = np.diff(model_altitude_sorted, axis=0)                            # (Nz-1, Ny, Nx)
    kext_mean = 0.5 * (kext_sorted[:-1] + kext_sorted[1:])            # (Nz-1, Ny, Nx)
    d_tau     = kext_mean * dz                                         # (Nz-1, Ny, Nx)

    cum_tau          = np.zeros((Nz, Ny, Nx))
    #cum_tau[1:, ...] = np.cumsum(d_tau, axis=0)                        # (Nz, Ny, Nx)
    # Where d_tau is NaN, use 0, otherwise use d_tau
    cum_tau[1:, ...] = np.cumsum(np.where(np.isnan(d_tau), 0, d_tau), axis=0)

    # ------------------------------------------------------------------
    # 5. Find the level closest to the radar altitude in each column
    # ------------------------------------------------------------------
    idx_radar     = np.abs(model_altitude_sorted - radar_altitude).argmin(axis=0) # (Ny, Nx)
    radar_cum_val = cum_tau[idx_radar, y_idx[0], x_idx[0]]             # (Ny, Nx)

    # ------------------------------------------------------------------
    # 6. Two-way attenuation
    # ------------------------------------------------------------------
    att        = np.abs(cum_tau - radar_cum_val[np.newaxis, ...])      # (Nz, Ny, Nx)
    zh_att_lin = zh_lin_sorted * np.exp(-2.0 * att)

    # ------------------------------------------------------------------
    # 7. Convert back to dBZ
    # ------------------------------------------------------------------
    zh_att_sorted = np.full((Nz, Ny, Nx), np.nan)
    valid = zh_att_lin > 1e-15
    zh_att_sorted[valid] = 10.0 * np.log10(zh_att_lin[valid])

    # ------------------------------------------------------------------
    # 8. Invert the sort → restore the original vertical ordering
    # ------------------------------------------------------------------
    unsort_idx = np.argsort(sort_idx, axis=0)
    zh_att     = zh_att_sorted[unsort_idx, y_idx, x_idx]

    return zh_att


# ======================================================================
# W-band empirical relationships (94 GHz)
# Based on Lhermitte (1990), Matrosov (2007), Kollias et al. (2007)
# ======================================================================

def rain_zh_wband(RR):
    """
    Rain reflectivity at W-band (dBZ).
    Z = 22 * R^1.4  (Lhermitte 1990, mm^6/m^3 vs mm/h)
    """
    return 10.0 * np.log10(22.0 * RR ** 1.4 + 1e-15)


def rain_kext_wband(RR):
    """
    Rain extinction at W-band (m⁻¹).
    k = 0.0238 * R^0.824  (dB/km) → convert to m⁻¹
    k [dB/km] → k [m⁻¹] : multiply by ln(10)/10 / 1000
    Reference: Matrosov (2007)
    """
    k_db_km = 0.0238 * RR ** 0.824   # dB/km  -- one-way
    return k_db_km * (np.log(10.0) / 10.0) / 1000.0   # m⁻¹


def cloud_kext_wband(LWC):
    """
    Cloud liquid water extinction at W-band (m⁻¹).
    k ≈ 0.0096 * LWC  [dB/km per g/m³]  (Lhermitte 1990)
    LWC in g/m³.
    """
    k_db_km = 0.0096 * LWC
    return k_db_km * (np.log(10.0) / 10.0) / 1000.0   # m⁻¹


def cloud_zh_wband(LWC):
    """
    Cloud droplet reflectivity at W-band (dBZ).
    Z ≈ 0.031 * LWC^1.56  (Atlas 1954 adapted for W-band, mm^6/m^3)
    """
    return 10.0 * np.log10(0.031 * LWC ** 1.56 + 1e-15)


def wvapor_kext_wband(height_m):
    """
    Water vapor continuum extinction at W-band (m⁻¹).
    Approximated as exponential decay with scale height ~2 km.
    Surface value ~0.003 dB/km (mid-latitude standard atmosphere).
    """
    k_db_km = 0.003 * np.exp(-height_m / 2000.0)
    return k_db_km * (np.log(10.0) / 10.0) / 1000.0   # m⁻¹


# ======================================================================
# Scene builder
# ======================================================================

def build_wband_scene(Nz=80, Ny=1, Nx=120):
    """
    Construct a 2-D (z, x) W-band scene with:
      - Water vapor background (all heights)
      - Cloud layer (2–4 km)
      - Stratiform rain column (0–4 km, x = 20–60 km)
      - Convective rain column (0–6 km, x = 80–100 km)
      - Melting layer bright-band (4–4.5 km)

    Returns
    -------
    model_altitude, zh, kext : each (Nz, Ny, Nx)
    z1d, x1d      : 1-D coordinate arrays
    scene_info    : dict with layer boundaries for annotation
    """
    rng = np.random.default_rng(0)

    z1d = np.linspace(0, 12000, Nz)
    x1d = np.linspace(0, 200e3, Nx)

    model_altitude  = np.zeros((Nz, Ny, Nx))
    zh   = np.full((Nz, Ny, Nx), -30.0)
    kext = np.zeros((Nz, Ny, Nx))

    # Broadcast z and x to full (Nz, Ny, Nx) shape from the start
    z3 = np.broadcast_to(z1d[:, None, None], (Nz, Ny, Nx)).copy()
    x3 = np.broadcast_to(x1d[None, None, :], (Nz, Ny, Nx)).copy()

    # ------------------------------------------------------------------
    # 1. Water vapor background
    # ------------------------------------------------------------------
    kext += wvapor_kext_wband(z3)

    # ------------------------------------------------------------------
    # 2. Cloud layer (2–4 km), LWC = 0.1–0.3 g/m³
    # ------------------------------------------------------------------
    cloud_mask = (z3 >= 2000) & (z3 <= 4000)
    LWC_cloud  = 0.2 + 0.1 * np.sin(np.pi * (z3 - 2000) / 2000)

    zh_cloud   = cloud_zh_wband(LWC_cloud)
    kext_cloud = cloud_kext_wband(LWC_cloud)

    zh[cloud_mask]   = np.maximum(
        zh[cloud_mask],
        zh_cloud[cloud_mask] + rng.normal(0, 0.5, cloud_mask.sum())
    )
    kext[cloud_mask] += kext_cloud[cloud_mask]

    # ------------------------------------------------------------------
    # 3. Stratiform rain (x = 20–60 km, z = 0–4 km)
    # ------------------------------------------------------------------
    strat_mask = (x3 >= 20e3) & (x3 <= 60e3) & (z3 <= 4000)
    RR_strat   = np.clip(
        4.0 * (1.0 - z3 / 4000.0) + rng.normal(0, 0.3, (Nz, Ny, Nx)),
        0.1, None
    )

    zh_strat   = rain_zh_wband(RR_strat)
    kext_strat = rain_kext_wband(RR_strat)

    zh[strat_mask]    = np.maximum(zh[strat_mask], zh_strat[strat_mask])
    kext[strat_mask] += kext_strat[strat_mask]

    # ------------------------------------------------------------------
    # 4. Convective rain (x = 80–100 km, z = 0–6 km)
    # ------------------------------------------------------------------
    conv_mask = (x3 >= 80e3) & (x3 <= 100e3) & (z3 <= 6000)
    RR_conv   = np.clip(
        50.0 * np.exp(-z3 / 3000.0)
        * (1.0 - np.abs(x3 - 90e3) / 10e3)
        + rng.normal(0, 2.0, (Nz, Ny, Nx)),
        0.5, None
    )

    zh_conv   = rain_zh_wband(RR_conv)
    kext_conv = rain_kext_wband(RR_conv)

    zh[conv_mask]    = np.maximum(zh[conv_mask], zh_conv[conv_mask])
    kext[conv_mask] += kext_conv[conv_mask]

    # ------------------------------------------------------------------
    # 5. Melting layer bright-band (4.0–4.5 km)
    # ------------------------------------------------------------------
    bb_mask   = (z3 >= 4000) & (z3 <= 4500)
    precip_x  = (x3 >= 20e3) & (x3 <= 60e3)
    bb_active = bb_mask & precip_x

    kext_bb = rain_kext_wband(np.full((Nz, Ny, Nx), 5.0))

    zh[bb_active]    += 8.0
    kext[bb_active]  += kext_bb[bb_active]

    # ------------------------------------------------------------------
    # 6. Altitude array
    # ------------------------------------------------------------------
    model_altitude[:] = z3

    scene_info = {
        "cloud_z"  : (2000, 4000),
        "strat_x"  : (20,   60),
        "conv_x"   : (80,  100),
        "bb_z"     : (4.0,  4.5),
        "radar_alt": 0.0,
    }
    return model_altitude, zh, kext, z1d, x1d, scene_info


# ======================================================================
# Analytical reference
# ======================================================================

def analytical_att_dbz(zh0, kext_const, z_levels, z_radar):
    """Two-way attenuation for constant kext (dBZ)."""
    att_db = (20.0 * kext_const / np.log(10.0)) * np.abs(z_levels - z_radar)
    return zh0 - att_db


# ======================================================================
# Plot helpers  (private, require matplotlib)
# ======================================================================

def _check_mpl():
    """Raise a clean error if matplotlib was not imported."""
    if not _MATPLOTLIB_AVAILABLE:
        raise RuntimeError(
            "matplotlib is required for plotting but could not be imported. "
            "Install it with:  pip install matplotlib"
        )


def _style_ax(ax, xlabel=None, ylabel=None, title=None,
              xlim=None, ylim=None):
    """Apply consistent dark-theme styling to an axis."""
    _check_mpl()
    if xlabel: ax.set_xlabel(xlabel, fontsize=FONTSIZE_LABEL,  color="white")
    if ylabel: ax.set_ylabel(ylabel, fontsize=FONTSIZE_LABEL,  color="white")
    if title:  ax.set_title(title,   fontsize=FONTSIZE_TITLE,
                            fontweight="bold", pad=8, color="white")
    if xlim:   ax.set_xlim(xlim)
    if ylim:   ax.set_ylim(ylim)
    ax.tick_params(labelsize=FONTSIZE_TICK, colors="white")
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ax.spines.values():
        spine.set_edgecolor("grey")
    return ax


def _add_colorbar(fig, ax, cf, label, ticks=None):
    """Add a well-sized colorbar with label."""
    _check_mpl()
    cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=25, shrink=0.95)
    cb.set_label(label, fontsize=FONTSIZE_CBAR, color="white")
    cb.ax.tick_params(labelsize=FONTSIZE_TICK, colors="white")
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
    if ticks is not None:
        cb.set_ticks(ticks)
    return cb


def _add_scene_overlays(ax, info, x_km, z_km, show_legend=False):
    """
    Draw layer boundaries and region markers as clean overlays.
    Returns legend handles for optional external legend.
    """
    _check_mpl()
    z_max = z_km[-1]

    # Melting layer shading
    ax.axhspan(info["bb_z"][0], info["bb_z"][1],
               color=BB_COLOR, alpha=0.25, zorder=2)

    # Cloud top / base dashed lines
    for z_val in info["cloud_z"]:
        ax.axhline(z_val / 1e3, color=CLOUD_COLOR,
                   lw=1.8, ls="--", alpha=0.9, zorder=3)

    # Stratiform region
    ax.axvspan(info["strat_x"][0], info["strat_x"][1],
               color=STRAT_COLOR, alpha=0.08, zorder=1)
    ax.axvline(info["strat_x"][0], color=STRAT_COLOR,
               lw=1.5, ls=":", alpha=0.9, zorder=3)
    ax.axvline(info["strat_x"][1], color=STRAT_COLOR,
               lw=1.5, ls=":", alpha=0.9, zorder=3)

    # Convective region
    ax.axvspan(info["conv_x"][0], info["conv_x"][1],
               color=CONV_COLOR, alpha=0.08, zorder=1)
    ax.axvline(info["conv_x"][0], color=CONV_COLOR,
               lw=1.5, ls=":", alpha=0.9, zorder=3)
    ax.axvline(info["conv_x"][1], color=CONV_COLOR,
               lw=1.5, ls=":", alpha=0.9, zorder=3)

    # Radar altitude line
    ax.axhline(info["radar_alt"] / 1e3, color=RADAR_COLOR,
               lw=2.0, ls="-.", alpha=0.95, zorder=4)

    # Region text labels
    y_label = z_max * 0.93
    for x_ctr, label, color in [
        (sum(info["strat_x"]) / 2, "Stratiform",  STRAT_COLOR),
        (sum(info["conv_x"])  / 2, "Convective",  CONV_COLOR),
    ]:
        ax.text(
            x_ctr, y_label, label,
            color=color, fontsize=9, fontweight="bold",
            ha="center", va="top",
            bbox=dict(fc="black", alpha=0.35, boxstyle="round,pad=0.2"),
            zorder=5
        )

    # Legend handles (returned for optional use)
    handles = [
        Line2D([0], [0], color=CLOUD_COLOR, lw=1.8, ls="--",
               label="Cloud top / base"),
        Line2D([0], [0], color=RADAR_COLOR, lw=2.0, ls="-.",
               label=f"Radar alt. ({info['radar_alt']:.0f} m)"),
        Rectangle((0, 0), 1, 1, fc=BB_COLOR,    alpha=0.5,
                  label="Melting layer"),
        Rectangle((0, 0), 1, 1, fc=STRAT_COLOR, alpha=0.3,
                  label="Stratiform region"),
        Rectangle((0, 0), 1, 1, fc=CONV_COLOR,  alpha=0.3,
                  label="Convective region"),
    ]
    if show_legend:
        leg = ax.legend(handles=handles, fontsize=FONTSIZE_LEGEND,
                        loc="upper right", framealpha=0.75,
                        edgecolor="grey", ncol=2,
                        facecolor=BG_PANEL)
        for text in leg.get_texts():
            text.set_color("white")
    return handles


def _dark_contour_panel(fig, ax, x_km, z_km, data,
                        levels, cmap, norm,
                        cbar_label, title, cbar_ticks=None):
    """
    Filled-contour panel with dark theme, white grid, colorbar.
    Returns the QuadContourSet for further use if needed.
    """
    _check_mpl()
    ax.set_facecolor(BG_PANEL)
    cf = ax.contourf(x_km, z_km, data,
                     levels=levels, cmap=cmap, norm=norm,
                     extend="both", zorder=0)
    ax.contour(x_km, z_km, data,
               levels=levels[::2], colors="white",
               linewidths=0.4, alpha=0.25, zorder=1)
    _add_colorbar(fig, ax, cf, cbar_label, ticks=cbar_ticks)
    _style_ax(ax, ylabel="Altitude (km)", title=title,
              xlim=(x_km[0], x_km[-1]), ylim=(0, 12))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax.xaxis.set_major_locator(mticker.MultipleLocator(25))
    ax.grid(which="major", color="white", alpha=0.12, lw=0.6, zorder=0)
    ax.grid(which="minor", color="white", alpha=0.05, lw=0.3, zorder=0)
    return cf


# ======================================================================
# Public plot functions
# ======================================================================

def plot_wband_scene():
    """
    Full W-band cross-section figure (dark theme, 5 panels):
      Row 0 : Unattenuated Zh          (full width)
      Row 1 : Attenuated  Zh           (full width)
      Row 2 : Two-way attenuation (dB) (full width)
      Row 3L: Extinction coefficient   (half width)
      Row 3R: Vertical profiles        (half width)
    """
    _check_mpl()

    # ---- build scene --------------------------------------------------
    alt, zh, kext, z1d, x1d, info = build_wband_scene()
    radar_alt = info["radar_alt"]
    zh_att    = compute_attenuated_zh_3D(kext, model_altitude, zh, radar_alt)

    x_km   = x1d / 1e3
    z_km   = z1d / 1e3
    zh2d   = zh[:, 0, :]
    zha2d  = zh_att[:, 0, :]
    att2d  = zh2d - zha2d
    kext2d = kext[:, 0, :] * 1e4    # ×10⁻⁴ m⁻¹

    # ---- colour maps & norms ------------------------------------------
    zh_levels  = np.arange(-20, 51, 5)
    zh_cmap    = plt.cm.turbo
    zh_norm    = mcolors.BoundaryNorm(zh_levels, zh_cmap.N)

    att_levels = np.arange(0, 42, 2)
    att_cmap   = plt.cm.afmhot_r
    att_norm   = mcolors.BoundaryNorm(att_levels, att_cmap.N)

    kext_levels = np.linspace(0, 6, 25)
    kext_cmap   = plt.cm.plasma

    # ---- figure layout ------------------------------------------------
    fig = plt.figure(figsize=(20, 22), facecolor=BG_DARK)
    gs  = gridspec.GridSpec(
        4, 2, figure=fig,
        height_ratios=[1, 1, 1, 1.15],
        hspace=0.45, wspace=0.22,
        left=0.07, right=0.97,
        top=0.93, bottom=0.05
    )

    ax_zh   = fig.add_subplot(gs[0, :])
    ax_zha  = fig.add_subplot(gs[1, :])
    ax_att  = fig.add_subplot(gs[2, :])
    ax_kext = fig.add_subplot(gs[3, 0])
    ax_prof = fig.add_subplot(gs[3, 1])

    # ==================================================================
    # Panel 0 – Unattenuated Zh
    # ==================================================================
    _dark_contour_panel(
        fig, ax_zh, x_km, z_km, zh2d,
        zh_levels, zh_cmap, zh_norm,
        cbar_label="dBZ",
        title="⬆  Unattenuated reflectivity   $Z_H$",
        cbar_ticks=zh_levels[::2]
    )
    _add_scene_overlays(ax_zh, info, x_km, z_km, show_legend=True)

    # ==================================================================
    # Panel 1 – Attenuated Zh
    # ==================================================================
    _dark_contour_panel(
        fig, ax_zha, x_km, z_km, zha2d,
        zh_levels, zh_cmap, zh_norm,
        cbar_label="dBZ",
        title="📡  Attenuated reflectivity   $Z_H^{att}$  (ground radar)",
        cbar_ticks=zh_levels[::2]
    )
    _add_scene_overlays(ax_zha, info, x_km, z_km)

    # ==================================================================
    # Panel 2 – Two-way attenuation
    # ==================================================================
    _dark_contour_panel(
        fig, ax_att, x_km, z_km, att2d,
        att_levels, att_cmap, att_norm,
        cbar_label="dB",
        title="🔥  Two-way attenuation   $Z_H - Z_H^{att}$",
        cbar_ticks=att_levels[::4]
    )
    _add_scene_overlays(ax_att, info, x_km, z_km)
    ax_att.set_xlabel("Range (km)", fontsize=FONTSIZE_LABEL, color="white")

    # ==================================================================
    # Panel 3L – Extinction coefficient
    # ==================================================================
    ax_kext.set_facecolor(BG_PANEL)
    cf_k = ax_kext.contourf(x_km, z_km, kext2d,
                             levels=kext_levels, cmap=kext_cmap,
                             extend="max", zorder=0)
    ax_kext.contour(x_km, z_km, kext2d,
                    levels=kext_levels[::3], colors="white",
                    linewidths=0.4, alpha=0.2, zorder=1)
    _add_colorbar(fig, ax_kext, cf_k,
                  label=r"$k_{ext}$  (×10⁻⁴ m⁻¹)")
    _add_scene_overlays(ax_kext, info, x_km, z_km)
    _style_ax(ax_kext,
              xlabel="Range (km)", ylabel="Altitude (km)",
              title="Extinction coefficient  $k_{ext}$",
              xlim=(x_km[0], x_km[-1]), ylim=(0, 12))
    ax_kext.yaxis.set_major_locator(mticker.MultipleLocator(2))
    ax_kext.xaxis.set_major_locator(mticker.MultipleLocator(25))
    ax_kext.grid(color="white", alpha=0.1, lw=0.5)

    # ==================================================================
    # Panel 3R – Vertical profiles
    # ==================================================================
    ax_prof.set_facecolor(BG_PROFILE)

    ix_conv  = np.argmin(np.abs(x1d - 90e3))
    ix_strat = np.argmin(np.abs(x1d - 40e3))
    ix_clear = np.argmin(np.abs(x1d -  5e3))

    profile_specs = [
        (ix_conv,  "Convective", CONV_COLOR,  "-",  "--", 2.5),
        (ix_strat, "Stratiform", STRAT_COLOR, "-",  "--", 2.0),
        (ix_clear, "Clear air",  "#87CEEB",   "-",  "--", 1.5),
    ]

    for ix, label, color, ls_u, ls_a, lw in profile_specs:
        ax_prof.plot(zh2d[:, ix],  z_km, color=color, lw=lw,
                     ls=ls_u, label=f"{label} (unatten.)", alpha=0.95)
        ax_prof.plot(zha2d[:, ix], z_km, color=color, lw=lw,
                     ls=ls_a, label=f"{label} (atten.)",  alpha=0.65)

    # Shading between unattenuated and attenuated for convective column
    ax_prof.fill_betweenx(
        z_km,
        zh2d[:, ix_conv], zha2d[:, ix_conv],
        color=CONV_COLOR, alpha=0.15,
        label="Atten. deficit (conv.)"
    )

    # Layer markers
    for z_val in info["cloud_z"]:
        ax_prof.axhline(z_val / 1e3, color=CLOUD_COLOR,
                        lw=1.4, ls="--", alpha=0.7)
    ax_prof.axhspan(*info["bb_z"], color=BB_COLOR,
                    alpha=0.2, label="Melting layer")
    ax_prof.axhline(info["radar_alt"] / 1e3, color=RADAR_COLOR,
                    lw=1.8, ls="-.", alpha=0.9, label="Radar alt.")
    ax_prof.axvline(0, color="white", lw=0.8, alpha=0.3)

    _style_ax(ax_prof,
              xlabel="Reflectivity (dBZ)", ylabel="Altitude (km)",
              title="Vertical profiles",
              xlim=(-30, 55), ylim=(0, 12))
    ax_prof.xaxis.set_major_locator(mticker.MultipleLocator(10))
    ax_prof.yaxis.set_major_locator(mticker.MultipleLocator(2))
    ax_prof.yaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax_prof.grid(which="major", color="white", alpha=0.12, lw=0.6)
    ax_prof.grid(which="minor", color="white", alpha=0.05, lw=0.3)

    leg = ax_prof.legend(fontsize=8.5, loc="upper right",
                         framealpha=0.6, edgecolor="grey",
                         facecolor=BG_PANEL, ncol=1)
    for text in leg.get_texts():
        text.set_color("white")

    # ==================================================================
    # Super title
    # ==================================================================
    fig.text(
        0.5, 0.965,
        "W-band (94 GHz) radar  —  Simulated attenuation cross-section",
        ha="center", va="center",
        fontsize=18, fontweight="bold", color="white",
        fontfamily="monospace"
    )
    fig.text(
        0.5, 0.950,
        "Ground-based geometry  |  Convective + Stratiform + Cloud + Melting layer",
        ha="center", va="center",
        fontsize=12, color="#aaaacc", style="italic"
    )

    plt.savefig("wband_cross_section.png", dpi=160,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show()
    print("Figure saved → wband_cross_section.png")


def plot_wband_analytical(z1d, zh_ref, zh_num, zh0, kext_val, z_radar):
    """
    Side-by-side validation figure:
      Left  – reflectivity profile (analytical vs numerical)
      Right – two-way attenuation profile
    Dark theme, large fonts, shaded attenuation deficit.
    """
    _check_mpl()

    att_ref = zh0 - zh_ref
    att_num = zh0 - zh_num
    z_km    = z1d / 1e3

    fig, axes = plt.subplots(
        1, 2, figsize=(14, 8), sharey=True,
        facecolor=BG_DARK
    )
    fig.subplots_adjust(wspace=0.12, left=0.09, right=0.97,
                        top=0.88, bottom=0.10)

    for ax in axes:
        ax.set_facecolor(BG_PANEL)
        ax.tick_params(labelsize=FONTSIZE_TICK, colors="white")
        ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
        ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))
        ax.grid(which="major", color="white", alpha=0.12, lw=0.6)
        ax.grid(which="minor", color="white", alpha=0.05, lw=0.3)
        ax.spines[["top", "right"]].set_visible(False)
        for spine in ax.spines.values():
            spine.set_edgecolor("grey")

    # ------------------------------------------------------------------
    # Left panel – reflectivity
    # ------------------------------------------------------------------
    ax = axes[0]
    ax.plot(zh_ref, z_km, color="white",   lw=3,  ls="--",
            label="Analytical", zorder=3)
    ax.plot(zh_num, z_km, color="#FF6B6B", lw=0,  marker="o",
            ms=5, label="Numerical", zorder=4, alpha=0.85)
    ax.fill_betweenx(z_km, zh_ref, zh_num,
                     color="#FF6B6B", alpha=0.15,
                     label="Residual error")
    ax.axhline(z_radar / 1e3, color=RADAR_COLOR, lw=2, ls="-.",
               label=f"Radar ({z_radar:.0f} m)", zorder=5)

    # Annotation: attenuation at top of column
    ax.annotate(
        f"Δ = {zh0 - zh_num[-1]:.1f} dB\nat 10 km",
        xy=(zh_num[-1], z_km[-1]),
        xytext=(zh_num[-1] + 4, z_km[-1] - 1.0),
        color="white", fontsize=10,
        arrowprops=dict(arrowstyle="->", color="white", lw=1.2),
        bbox=dict(fc=BG_PANEL, ec="grey", boxstyle="round,pad=0.3")
    )

    _style_ax(ax,
              xlabel="Reflectivity (dBZ)",
              ylabel="Altitude (km)",
              title=f"Reflectivity profile\n"
                    f"$Z_{{H,0}}$ = {zh0:.1f} dBZ  |  RR = 50 mm/h")
    ax.set_ylim(0, 11)
    leg = ax.legend(fontsize=FONTSIZE_LEGEND, loc="lower left",
                    framealpha=0.6, edgecolor="grey", facecolor=BG_PANEL)
    for t in leg.get_texts():
        t.set_color("white")

    # ------------------------------------------------------------------
    # Right panel – attenuation
    # ------------------------------------------------------------------
    ax = axes[1]
    ax.fill_betweenx(z_km, att_ref,
                     color="#FF8C00", alpha=0.25,
                     label="Attenuation (analytical)")
    ax.plot(att_ref, z_km, color="white",   lw=3,  ls="--",
            label="Analytical", zorder=3)
    ax.plot(att_num, z_km, color="#FF6B6B", lw=0,  marker="o",
            ms=5, label="Numerical", zorder=4, alpha=0.85)
    ax.axhline(z_radar / 1e3, color=RADAR_COLOR, lw=2, ls="-.",
               label=f"Radar ({z_radar:.0f} m)", zorder=5)
    ax.axvline(0, color="white", lw=0.8, alpha=0.3)

    # Info box
    k_dbkm = kext_val / (np.log(10.0) / 10.0) * 1000.0
    ax.text(
        0.97, 0.05,
        f"$k_{{ext}}$ = {kext_val*1e4:.2f}×10⁻⁴ m⁻¹\n"
        f"= {k_dbkm:.2f} dB km⁻¹ (one-way)\n"
        f"= {2*k_dbkm:.2f} dB km⁻¹ (two-way)",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=10, color="#aaaacc",
        bbox=dict(fc=BG_PANEL, ec="grey",
                  boxstyle="round,pad=0.4", alpha=0.8)
    )

    _style_ax(ax,
              xlabel="Two-way attenuation (dB)",
              title="Two-way attenuation\nAnalytical vs Numerical")
    leg = ax.legend(fontsize=FONTSIZE_LEGEND, loc="upper left",
                    framealpha=0.6, edgecolor="grey", facecolor=BG_PANEL)
    for t in leg.get_texts():
        t.set_color("white")

    # Super title
    fig.text(
        0.5, 0.94,
        "W-band (94 GHz)  —  Analytical validation  |  "
        "Heavy rain column (RR = 50 mm/h)",
        ha="center", fontsize=15, fontweight="bold", color="white",
        fontfamily="monospace"
    )

    plt.savefig("wband_analytical_check.png", dpi=160,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show()
    print("Figure saved → wband_analytical_check.png")


# ======================================================================
# Tests
# ======================================================================

def run_wband_tests():
    """
    Three quantitative tests for compute_attenuated_zh_3D at W-band.

    Returns
    -------
    Tuple passed directly to plot_wband_analytical.
    """
    print("=" * 65)
    print("W-band (94 GHz) attenuation tests")
    print("=" * 65)

    # ------------------------------------------------------------------
    # Test W1 – Heavy rain column, constant kext, analytical comparison
    # ------------------------------------------------------------------
    RR_heavy  = 50.0
    kext_rain = rain_kext_wband(RR_heavy)
    zh0_rain  = rain_zh_wband(RR_heavy)

    Nz, Ny, Nx = 100, 3, 3
    z1d     = np.linspace(0, 10000, Nz)
    model_altitude     = np.broadcast_to(z1d[:, None, None], (Nz, Ny, Nx)).copy()
    zh_in   = np.full((Nz, Ny, Nx), zh0_rain)
    kext_in = np.full((Nz, Ny, Nx), kext_rain)
    z_radar = 0.0

    zh_out  = compute_attenuated_zh_3D(kext_in, alt, zh_in, z_radar)
    zh_ref  = analytical_att_dbz(zh0_rain, kext_rain, z1d, z_radar)

    max_err  = np.nanmax(np.abs(zh_out[:, 1, 1] - zh_ref))
    att_10km = zh0_rain - zh_out[-1, 1, 1]

    status = "PASS" if max_err < 0.05 else "FAIL"
    print(f"\n{status}  W1 – Heavy rain (RR={RR_heavy} mm/h)")
    print(f"       kext      = {kext_rain*1e4:.3f} × 10⁻⁴ m⁻¹"
          f"  ({kext_rain*(np.log(10.0)/10.0)**-1*1000:.2f} dB/km one-way)")
    print(f"       zh_in     = {zh0_rain:.1f} dBZ")
    print(f"       Attenuation over 10 km = {att_10km:.1f} dB  (two-way)")
    print(f"       Max error vs analytical = {max_err:.4f} dBZ")

    # ------------------------------------------------------------------
    # Test W2 – Layered scene: cloud + rain, monotone attenuation check
    # ------------------------------------------------------------------
    Nz2, Ny2, Nx2 = 60, 1, 1
    z1d2 = np.linspace(0, 6000, Nz2)
    model_altitude2 = np.broadcast_to(z1d2[:, None, None], (Nz2, Ny2, Nx2)).copy()

    LWC_cloud  = 0.3
    kext_cloud = cloud_kext_wband(LWC_cloud)
    RR_light   = 2.0
    kext_light = rain_kext_wband(RR_light)

    kext2 = np.zeros((Nz2, Ny2, Nx2))
    kext2[z1d2 < 2000, 0, 0] = kext_light
    kext2[(z1d2 >= 2000) & (z1d2 < 4000), 0, 0] = kext_cloud

    zh2     = np.full((Nz2, Ny2, Nx2), 20.0)
    zh_out2 = compute_attenuated_zh_3D(kext2, model_altitude2, zh2, radar_altitude=0.0)

    att2    = zh2[:, 0, 0] - zh_out2[:, 0, 0]
    is_mono = np.all(np.diff(att2) >= -1e-10)
    status  = "PASS" if is_mono else "FAIL"
    print(f"\n{status}  W2 – Layered cloud+rain, attenuation monotone = {is_mono}")
    print(f"       Attenuation at 6 km = {att2[-1]:.2f} dB (two-way)")
    print(f"       kext rain  = {kext_light*1e4:.3f} × 10⁻⁴ m⁻¹")
    print(f"       kext cloud = {kext_cloud*1e4:.3f} × 10⁻⁴ m⁻¹")

    # ------------------------------------------------------------------
    # Test W3 – Airborne radar at 8 km
    # ------------------------------------------------------------------
    z_radar_aloft = 8000.0
    zh_out3 = compute_attenuated_zh_3D(kext_in, model_altitude, zh_in, z_radar_aloft)
    att3    = zh_in[:, 1, 1] - zh_out3[:, 1, 1]

    below        = z1d < z_radar_aloft
    above        = z1d > z_radar_aloft
    is_mono_down = np.all(np.diff(att3[below]) <= 1e-10)
    is_mono_up   = np.all(np.diff(att3[above]) >= -1e-10)

    status = "PASS" if (is_mono_down and is_mono_up) else "FAIL"
    print(f"\n{status}  W3 – Airborne radar at {z_radar_aloft/1e3:.0f} km")
    print(f"       Attenuation monotone below radar = {is_mono_down}")
    print(f"       Attenuation monotone above radar = {is_mono_up}")
    print(f"       Max two-way att (surface) = {att3[0]:.1f} dB")

    print("\n" + "=" * 65)
    return z1d, zh_ref, zh_out[:, 1, 1], zh0_rain, kext_rain, z_radar


# ======================================================================
# Entry point
# ======================================================================

if __name__ == "__main__":
    results = run_wband_tests()
    plot_wband_analytical(*results)
    plot_wband_scene()