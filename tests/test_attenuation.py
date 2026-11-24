#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Attenuation of a vertically pointing radar (operadar/radar/attenuation.py): two-way integration of the
extinction along the vertical from the radar, conversion of Ah, and gas attenuation (pyMPM).
These tests replace the W-band checks that were run as __main__ of attenuation.py."""
import numpy as np
import pytest

from operadar.radar.attenuation import (compute_attenuated_zh_3D, compute_extinction, compute_gaz_attenuation,
                                        compute_relative_humidity)

DB = 10. / np.log(10.)         # 4.343: one-way power extinction (m-1) -> dB/m
LAM = {"C": 0.0532, "X": 0.0319, "Ka": 0.0084, "W": 0.00316}   # wavelengths (m) of the tables


def column(values, ny=2, nx=3):
    """(Nz,) profile -> (Nz, ny, nx) field."""
    return np.broadcast_to(np.asarray(values, dtype=float)[:, None, None], (len(values), ny, nx)).copy()


def kext_from_db_km(a_db_km):
    """One-way attenuation (dB/km) -> extinction coefficient (m-1)."""
    return a_db_km / DB / 1000.


# ----------------------------------------------------------------------------------------------------------
# compute_attenuated_zh_3D
# ----------------------------------------------------------------------------------------------------------
def test_constant_extinction_matches_analytical_two_way_attenuation():
    # ground radar, 1 dB/km one-way everywhere: Zh_att = Zh - 2 dB/km x height
    z = np.linspace(0., 10000., 101)
    zh_att = compute_attenuated_zh_3D(column(np.full(z.size, kext_from_db_km(1.))), column(z),
                                      column(np.full(z.size, 30.)), radar_altitude=0.)
    np.testing.assert_allclose(zh_att, column(30. - 2. * z / 1000.), atol=1e-9)


def test_no_extinction_leaves_zh_unchanged():
    z = np.linspace(0., 5000., 51)
    zh = column(np.linspace(-20., 40., z.size))
    np.testing.assert_allclose(compute_attenuated_zh_3D(column(np.zeros(z.size)), column(z), zh, 0.), zh,
                               atol=1e-12)


def test_attenuating_layer_then_clear_air():
    # rain layer 0-2 km at 5 dB/km one-way, nothing above: 20 dB two-way above the layer (plateau),
    # increasing through the layer
    z = np.arange(0., 6001., 250.)
    kext = column(np.where(z <= 2000., kext_from_db_km(5.), 0.))
    att = 20. - compute_attenuated_zh_3D(kext, column(z), column(np.full(z.size, 20.)), 0.)[:, 0, 0]
    assert np.all(np.diff(att) >= -1e-12)
    np.testing.assert_allclose(att[z <= 2000.], 2. * 5. * z[z <= 2000.] / 1000., atol=1e-9)
    # the trapezoid between 2000 m (5 dB/km) and 2250 m (0) adds half a layer: 20 + 2 x 2.5 x 0.25
    np.testing.assert_allclose(att[z >= 2250.], 21.25, atol=1e-9)


def test_airborne_radar_attenuation_grows_away_from_the_radar_on_both_sides():
    z = np.linspace(0., 10000., 101)
    zh_att = compute_attenuated_zh_3D(column(np.full(z.size, kext_from_db_km(1.))), column(z),
                                      column(np.full(z.size, 30.)), radar_altitude=8000.)
    np.testing.assert_allclose(zh_att[:, 0, 0], 30. - 2. * np.abs(z - 8000.) / 1000., atol=1e-9)


def test_vertical_order_of_the_levels_does_not_matter():
    # AROME levels go from the top to the ground: the result must not depend on the order
    z = np.linspace(0., 8000., 41)
    rng = np.random.default_rng(1)
    kext = column(kext_from_db_km(rng.uniform(0., 3., z.size)))
    zh = column(rng.uniform(-10., 40., z.size))
    up = compute_attenuated_zh_3D(kext, column(z), zh, 0.)
    down = compute_attenuated_zh_3D(kext[::-1], column(z[::-1]), zh[::-1], 0.)
    np.testing.assert_allclose(down[::-1], up, atol=1e-10)


def test_columns_are_independent_with_terrain():
    # two columns with different ground altitudes and extinctions: each one is integrated on its own
    z0 = np.linspace(0., 4000., 21)
    alt = np.stack([z0, z0 + 1500.], axis=1)[:, :, None]                  # (Nz, 2, 1)
    kext = np.stack([np.full(z0.size, kext_from_db_km(1.)), np.full(z0.size, kext_from_db_km(4.))],
                    axis=1)[:, :, None]
    zh = np.full(alt.shape, 25.)
    zh_att = compute_attenuated_zh_3D(kext, alt, zh, radar_altitude=0.)
    # the radar is at the level closest to its altitude: the lowest level of each column
    np.testing.assert_allclose(zh_att[:, 0, 0], 25. - 2. * 1. * z0 / 1000., atol=1e-9)
    np.testing.assert_allclose(zh_att[:, 1, 0], 25. - 2. * 4. * z0 / 1000., atol=1e-9)


def test_nan_reflectivity_stays_nan_without_breaking_the_column():
    z = np.linspace(0., 3000., 31)
    zh = column(np.full(z.size, 20.))
    zh[10:15] = np.nan                                       # no hydrometeor there
    zh_att = compute_attenuated_zh_3D(column(np.full(z.size, kext_from_db_km(1.))), column(z), zh, 0.)
    assert np.all(np.isnan(zh_att[10:15]))
    np.testing.assert_allclose(zh_att[20:, 0, 0], 20. - 2. * z[20:] / 1000., atol=1e-9)


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        compute_attenuated_zh_3D(np.zeros((3, 2, 2)), np.zeros((3, 2, 2)), np.zeros((4, 2, 2)), 0.)


# ----------------------------------------------------------------------------------------------------------
# compute_extinction (Ah + gas)
# ----------------------------------------------------------------------------------------------------------
def atmosphere(nz=11, rh_surface=0.6):
    """Standard atmosphere 0-10 km: T (C), p (Pa), qv (kg/kg) with a given surface relative humidity."""
    z = np.linspace(0., 10000., nz)
    t_k = 288.15 - 0.0065 * z
    p = 101325. * (t_k / 288.15) ** 5.2559
    tc = t_k - 273.15
    e = rh_surface * np.exp(-z / 2000.) * 610.7 * 10. ** (7.5 * tc / (237.3 + tc))    # Pa
    qv = 0.622 * e / (p - 0.378 * e)
    return z, column(tc), column(p), column(qv)


@pytest.mark.parametrize("band", ["C", "W"])
def test_ah_is_converted_from_one_way_db_per_km(band):
    _, tc, p, qv = atmosphere()
    ah = np.full(tc.shape, 2.)
    gas = compute_extinction(tc, p, qv, LAM[band], np.zeros(tc.shape))
    total = compute_extinction(tc, p, qv, LAM[band], ah)
    np.testing.assert_allclose(total - gas, kext_from_db_km(2.), rtol=1e-4)   # code uses 4.343


def test_gas_extinction_is_kept_where_there_is_no_hydrometeor():
    # Ah is NaN outside the hydrometeors: the gas extinction must remain (it was lost before Nov. 2026)
    _, tc, p, qv = atmosphere()
    ah = np.full(tc.shape, np.nan)
    ah[:3] = 1.
    kext = compute_extinction(tc, p, qv, LAM["W"], ah)
    gas = compute_extinction(tc, p, qv, LAM["W"], np.zeros(tc.shape))
    assert np.all(np.isfinite(kext))
    np.testing.assert_allclose(kext[3:], gas[3:], rtol=1e-12)


# ----------------------------------------------------------------------------------------------------------
# Gas attenuation (pyMPM) and relative humidity
# ----------------------------------------------------------------------------------------------------------
def standard_surface(rho_v=7.5e-3, nz=2):
    """ITU-R P.676 reference: 1013.25 hPa, 15 C, water vapour density rho_v (kg/m3)."""
    tc, p = 15., 101325.
    rho_dry = p / (287.05 * (tc + 273.15))
    return column(np.full(nz, tc), 1, 1), column(np.full(nz, p), 1, 1), \
        column(np.full(nz, rho_v / (rho_dry + rho_v)), 1, 1)


@pytest.mark.parametrize("band, lo, hi", [("C", 0.006, 0.012),    # one-way dB/km, ITU-R P.676 orders
                                          ("Ka", 0.08, 0.16),     # of magnitude at 5.6, 35.7 and
                                          ("W", 0.35, 0.60)])     # 94.9 GHz
def test_gas_attenuation_of_the_standard_atmosphere(band, lo, hi):
    tc, p, qv = standard_surface()
    a_db_km = compute_gaz_attenuation(tc, p, qv, LAM[band]) * DB * 1000.
    assert a_db_km.shape == tc.shape
    assert np.all((a_db_km > lo) & (a_db_km < hi)), a_db_km.ravel()


def test_gas_attenuation_at_w_band_comes_mostly_from_water_vapour():
    tc, p, qv = standard_surface()
    moist = compute_gaz_attenuation(tc, p, qv, LAM["W"])
    dry = compute_gaz_attenuation(tc, p, np.zeros(qv.shape), LAM["W"])
    assert np.all(dry > 0.) and np.all(moist > 5. * dry)


def test_gas_attenuation_decreases_with_height():
    _, tc, p, qv = atmosphere()
    k = compute_gaz_attenuation(tc, p, qv, LAM["W"])[:, 0, 0]
    assert np.all(np.diff(k) < 0.)


def test_relative_humidity():
    tc = np.array([0., 15., 30.])
    p = np.full(3, 100000.)
    esat = 610.7 * 10. ** (7.5 * tc / (237.3 + tc))           # Pa (Tetens, as in the code)
    rv_sat = 0.622 * esat / (p - esat)
    np.testing.assert_allclose(compute_relative_humidity(tc, rv_sat, p), 100., rtol=1e-6)
    np.testing.assert_allclose(compute_relative_humidity(tc, 0.5 * rv_sat, p),
                               100. * 0.5 * rv_sat / (0.5 * rv_sat + 0.622) * p / esat, rtol=1e-6)
