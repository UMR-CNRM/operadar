#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:59:47 2026

@author: augros
"""

import numpy as np
import pyMPM

def compute_gaz_attenuation(temperature, pressure, qv, LAM):
    """
    Computes gaz attenuation using pyMPM

    Parameters
    ----------
    temperature : 3D array
        Temperature [K]
    pressure : 3D array
        Pressure [Pa]
    qv : 3D array
        Specific humidity [kg/kg]
    LAM: float
        radar wavelength (meters)

    Returns
    -------
    kext_gaz : 3D array
        Gas extinction coefficient
i    """

    # Convert specific humidity -> mixing ratio rv
    rv = qv / (1 - qv)

    # Compute relative humidity
    rh = compute_relative_humidity(temperature, rv, pressure)

    # Frequency in Hz (c / λ)
    freq = 3e8 / LAM  # Hz
    #freq_array= np.array([freq,])

    # Call pyMPM (requires frequency in GHz, pressure in hPa, temperature in °C)
    kext_gaz = (
        1e-3 / 4.343
        * pyMPM.MPM(
            freq * 1e-9,
            pressure * 1e-2,
            temperature - 273.15,
            rh,
            wa='None',
            wae='None',
            R='None',
            output_type='att'
        )
    )

    return kext_gaz


def compute_relative_humidity(temperature_K, rv, pressure_Pa):
    """
    Computes relative humidity (%) from:
        - temperature in Kelvin
        - water vapor mixing ratio rv (kg/kg)
        - pressure in Pascals

    Formula:
        Psat : Saturation vapor pressure using Tetens formula (hPa)
        Pe   : Partial vapor pressure (hPa)
        RH   : 100 * Pe / Psat

    Parameters
    ----------
    temperature_K : array-like or float
        Air temperature [K]
    rv : array-like or float
        Water vapor mixing ratio [kg/kg]
    pressure_Pa : array-like or float
        Atmospheric pressure [Pa]

    Returns
    -------
    rh : array-like or float
        Relative humidity in %
    """

    # Convert temperature to Celsius
    temperature_C = temperature_K - 273.15

    # Saturation vapor pressure (hPa)
    Psat = 6.107 * np.power(10., (7.5 * temperature_C) / (237.3 + temperature_C))

    # Partial vapor pressure (hPa)
    Pe = (rv * pressure_Pa / (rv + 0.622)) / 100.0

    # Relative humidity in %
    rh = 100.0 * (Pe / Psat)

    return rh


def compute_extinction(temperature: np.ndarray,
                pressure: np.ndarray,
                qv: np.ndarray,
                radar_lam:float,
                ah: np.ndarray,
                )-> np.ndarray :
    """Computes gaz and hydrometeor extinction cross-section 

    Args:
        temperature (np.ndarray): 3D array (°C)
        pressure (np.ndarray): 3D array (Pa)
        qv (np.ndarray): 3D array specific humidity (kg/kg)
        radar_lam: radar wavelength (m)
        ah : 3D array (dB/km)

    Returns:
        kext np.ndarray 3D array of extinction cross section
    """
    
    # Gaz extinction cross-section
    kext_gaz=compute_gaz_attenuation(temperature, pressure, qv, radar_lam)
    
    # Hydrometeor extinction cross-section
    kext_hydro=1e-3 * ah / 4.343

    # Total extinction cross-section
    kext_total = kext_hydro + kext_gaz
    
    return kext_total


def compute_attenuated_zh_3D(kext, zh, model_altitude, radar_altitude):
    """
    Compute attenuated radar reflectivity (dBZ) using vertical integration
    of extinction, for 3D fields with per-column altitude.

    Parameters
    ----------
    kext : ndarray (Nz, Ny, Nx)
        Extinction coefficient (m⁻¹), non-negative.
    zh : ndarray (Nz, Ny, Nx)
        Unattenuated reflectivity (dBZ).
    model_altitude : ndarray (Nz, Ny, Nx)
        Altitude of each grid point (m). Must be monotonic along the z-axis
        for each (y, x) column (either increasing or decreasing).
    radar_altitude : float
        Radar altitude (m).

    Returns
    -------
    zh_att : ndarray (Nz, Ny, Nx)
        Attenuated reflectivity (dBZ), including two-way attenuation.

    Notes
    -----
    - Integration is performed independently for each (y, x) column.
    - Attenuation is zero at the level closest to the radar altitude.
    - Only vertical attenuation is considered.
    """

    kext = np.asarray(kext)
    zh = np.asarray(zh)
    z = np.asarray(model_altitude)

    Nz, Ny, Nx = zh.shape

    # Convert reflectivity to linear scale
    zh_lin = 10 ** (zh / 10.0)

    # Ensure altitude increases with z for each column
    flip_mask = z[0, ...] > z[-1, ...]

    if np.any(flip_mask):
        z = z.copy()
        kext = kext.copy()
        zh_lin = zh_lin.copy()

        z[:, flip_mask] = z[::-1, flip_mask]
        kext[:, flip_mask] = kext[::-1, flip_mask]
        zh_lin[:, flip_mask] = zh_lin[::-1, flip_mask]

    # Compute dz (layer thickness)
    dz = np.diff(z, axis=0)  # (Nz-1, Ny, Nx)

    # Mean extinction between levels
    kext_mean = 0.5 * (kext[1:, ...] + kext[:-1, ...])

    # Find radar index per column (closest level)
    idx = np.abs(z - radar_altitude).argmin(axis=0)  # (Ny, Nx)

    # Initialize attenuation
    att = np.zeros_like(zh_lin)

    # Loop over vertical levels (vectorized over y,x)
    for k in range(Nz):
        # mask where this level is below radar
        below = k < idx
        above = k > idx

        if k > 0:
            seg = kext_mean[k-1, ...] * dz[k-1, ...]

            # accumulate downward
            att[k, ...][below] = att[k+1, ...][below] + seg[below] if k < Nz-1 else seg[below]

        if k < Nz - 1:
            seg = kext_mean[k, ...] * dz[k, ...]

            # accumulate upward
            att[k, ...][above] = att[k-1, ...][above] + seg[above] if k > 0 else seg[above]

    # Enforce zero attenuation at radar level
    for j in range(Ny):
        for i in range(Nx):
            att[idx[j, i], j, i] = 0.0

    # Apply two-way attenuation
    zh_att_lin = zh_lin * np.exp(-2.0 * att)

    # Convert back to dBZ
    zh_att = np.full_like(zh, np.nan)
    mask = zh_att_lin > 0
    zh_att[mask] = 10.0 * np.log10(zh_att_lin[mask])

    # Flip back columns if needed
    if np.any(flip_mask):
        zh_att[:, flip_mask] = zh_att[::-1, flip_mask]

    return zh_att    
