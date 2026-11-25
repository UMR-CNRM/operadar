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
  - compute_attenuated_zh_3D      : vectorized 3-D two-way attenuation along the vertical, for a
                                    vertically pointing radar at radar_altitude (ground or airborne)
Tests: tests/test_attenuation.py
"""

# ======================================================================
# Standard imports
# ======================================================================
import numpy as np

# pyMPM is now included in operadar
from operadar.radar.pyMPM.MPM import MPM


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

    # Hydrometeor extinction cross-section. Ah is NaN where there is no hydrometeor: no hydrometeor
    # extinction there, but the gas extinction must be kept (with NaN, the gas attenuation in clear air
    # and the half-layers at the cloud edges were ignored)
    kext_hydro = 1e-3 * np.nan_to_num(ah, nan=0.) / 4.343

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
