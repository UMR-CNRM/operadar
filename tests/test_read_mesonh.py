#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reading of MesoNH netCDF files (operadar/read/mesonh.py, operadar/read/with_netCDF4.py) on a small synthetic
file. Before Nov. 2026 MesoNH files could not be read at all (invalid annotation in with_netCDF4.py, functions
called without their module), the last row and column of the domain were dropped, and qv was not read."""
import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")
mesonh = pytest.importorskip("operadar.read.mesonh")

NZ, NJ, NI = 4, 6, 5   # levels, rows (latitude), columns (longitude)
MOMENTS = {'cc': 1, 'rr': 1, 'ss': 1, 'gg': 1, 'ii': 1, 'wg': 1}


@pytest.fixture(scope="module")
def mnh_file(tmp_path_factory):
    path = tmp_path_factory.mktemp("mnh") / "fake.MXX.nc"
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("time", 1)
        d.createDimension("level", NZ)
        d.createDimension("nj", NJ)
        d.createDimension("ni", NI)
        d.createVariable("XHAT", "f8", ("ni",))[:] = np.arange(NI) * 2500.
        d.createVariable("YHAT", "f8", ("nj",))[:] = np.arange(NJ) * 2500.
        d.createVariable("ZHAT", "f8", ("level",))[:] = np.arange(NZ) * 1000.
        d.createVariable("ZS", "f8", ("nj", "ni"))[:] = 10.
        lon, lat = np.meshgrid(50. + np.arange(NI), -20. + np.arange(NJ))
        d.createVariable("longitude", "f8", ("nj", "ni"))[:] = lon
        d.createVariable("latitude", "f8", ("nj", "ni"))[:] = lat
        d.createVariable("RHOREFZ", "f8", ("level",))[:] = np.linspace(1.2, 0.8, NZ)
        shape = (1, NZ, NJ, NI)
        d.createVariable("PABST", "f8", ("time", "level", "nj", "ni"))[:] = np.full(shape, 90000.)
        d.createVariable("THT", "f8", ("time", "level", "nj", "ni"))[:] = np.full(shape, 300.)
        # water vapour: a distinct value at each point, to check that it is read
        d.createVariable("RVT", "f8", ("time", "level", "nj", "ni"))[:] = \
            1e-3 * (1. + np.arange(NZ * NJ * NI).reshape(shape) / 1000.)
        for name in ("RCT", "RRT", "RIT", "RST", "RGT", "RHT"):
            d.createVariable(name, "f8", ("time", "level", "nj", "ni"))[:] = np.full(shape, 1e-4)
        cit = d.createVariable("CIT", "f8", ("time", "level", "nj", "ni"))
        cit[:] = np.full(shape, 1e5)
        cit.units = "m-3"
        # LIMA prognostic concentrations (kg-1)
        for name, value in (("CRAIN", 2e4), ("CICE", 3e5), ("CCLOUD", 1e8)):
            var = d.createVariable(name, "f8", ("time", "level", "nj", "ni"))
            var[:] = np.full(shape, value)
            var.units = "kg-1"
    return path


def read(path, subDomain=None, micro="ICE3", moments=MOMENTS):
    return mesonh.read_mesonh(filePath=path, micro=micro, subDomain=subDomain, hydrometeorMoments=moments,
                              real_case=True, verbose=False)


def test_whole_domain_keeps_every_point(mnh_file):
    X, Y, Z, lon, lat, M, Nc, Tc, p, qv = read(mnh_file)
    assert Tc.shape == (NZ, NJ, NI)
    assert lon.shape == (NJ, NI) and X.shape == (NI,) and Y.shape == (NJ,)
    assert M['rr'].shape == (NZ, NJ, NI)


def test_water_vapour_is_read(mnh_file):
    *_, qv = read(mnh_file)
    with netCDF4.Dataset(mnh_file) as d:
        expected = d.variables["RVT"][0]
    np.testing.assert_allclose(qv, expected)


def test_subdomain_keeps_its_last_row_and_column(mnh_file):
    # points strictly inside: lon 51, 52, 53 and lat -19, -18
    X, Y, Z, lon, lat, M, Nc, Tc, p, qv = read(mnh_file, subDomain=[50.5, 53.5, -19.5, -17.5])
    assert lon.shape == (2, 3)
    np.testing.assert_allclose(lon[0], [51., 52., 53.])
    np.testing.assert_allclose(lat[:, 0], [-19., -18.])
    assert Tc.shape == (NZ, 2, 3) and qv.shape == (NZ, 2, 3)


# ----------------------------------------------------------------------------------------------------------
# Number concentrations: chosen from the number of moments of each species, not from the name of the scheme
# ----------------------------------------------------------------------------------------------------------
LIMA = {'cc': 2, 'rr': 2, 'ss': 1, 'gg': 1, 'ii': 2, 'wg': 1}
LIMC = {'cc': 2, 'rr': 2, 'ss': 1, 'gg': 1, 'ii': 1, 'wg': 1}    # "LIMA chaud": 2-moment warm species only


def air_density(Nc):
    """Dry air density used by the reader, from a concentration in kg-1 (CRAIN = 2e4 kg-1)."""
    return Nc['rr'] / 2e4


@pytest.mark.parametrize("micro", ["ICE3", "ICJW"])
def test_one_moment_ice_concentration_is_cit_in_m3(mnh_file, micro):
    # CIT is in m-3: it must not be multiplied by the air density (it was before Nov. 2026)
    Nc = read(mnh_file, micro=micro)[6]
    np.testing.assert_allclose(Nc['ii'], 1e5)
    assert np.all(Nc['rr'] == 0.)


def test_lima_reads_the_prognostic_concentrations_in_kg(mnh_file):
    Nc = read(mnh_file, micro="LIMA", moments=LIMA)[6]
    rho = air_density(Nc)
    assert np.all((rho > 0.5) & (rho < 1.5))
    np.testing.assert_allclose(Nc['ii'], 3e5 * rho)
    np.testing.assert_allclose(Nc['cc'], 1e8 * rho)


def test_limc_two_moment_warm_species_one_moment_ice(mnh_file):
    Nc = read(mnh_file, micro="LIMC", moments=LIMC)[6]
    rho = air_density(Nc)
    assert np.all((rho > 0.5) & (rho < 1.5))
    np.testing.assert_allclose(Nc['cc'], 1e8 * rho)
    np.testing.assert_allclose(Nc['ii'], 1e5)          # CIT, not CICE


def test_missing_two_moment_concentration_raises(tmp_path, mnh_file):
    import shutil
    path = tmp_path / "no_crain.nc"
    shutil.copy(mnh_file, path)
    with netCDF4.Dataset(path, "a") as d:
        d.renameVariable("CRAIN", "OTHER")
    with pytest.raises(KeyError):
        read(path, micro="LIMC", moments=LIMC)


def test_one_moment_ice_without_cit_uses_800(tmp_path, mnh_file):
    import shutil
    path = tmp_path / "no_cit.nc"
    shutil.copy(mnh_file, path)
    with netCDF4.Dataset(path, "a") as d:
        d.renameVariable("CIT", "OTHER")
    np.testing.assert_allclose(read(path, micro="LIMC", moments=LIMC)[6]['ii'], 800.)
