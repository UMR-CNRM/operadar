"""Conversion of the summed linear reflectivities into polarimetric variables."""
import numpy as np

from operadar.radar.dualpol_variables import compute_dpol_var


def test_zero_linear_reflectivity_is_not_zero_dbz():
    # point 0: no signal from the tables (Zhhlin = 0), point 1: 100 mm6/m3 with Zdr = 3 dB, point 2: masked
    zhh = np.array([0., 100., np.nan])
    zvv = np.array([0., 100. / 10 ** 0.3, np.nan])
    out = compute_dpol_var({"Zhhlin": zhh, "Zvvlin": zvv, "Kdp": np.array([0., 0.1, np.nan])},
                           ["Zh", "Zdr", "Kdp"])
    assert np.isnan(out["Zh"][0]), "a zero linear reflectivity must not be written as 0 dBZ"
    assert np.isnan(out["Zdr"][0])
    np.testing.assert_allclose(out["Zh"][1], 20.)
    np.testing.assert_allclose(out["Zdr"][1], 3.)
    assert np.isnan(out["Zh"][2]) and np.isnan(out["Zdr"][2])


def rhohv_of_sum(species):
    """Sum the scattering terms of several species as operadar does, then compute rhohv."""
    keys = ["REdeltaco", "IMdeltaco", "sighh", "sigvv"]
    summed = {k: np.array([sum(s[k] for s in species)]) for k in keys}
    return compute_dpol_var(summed, ["Rhohv"])["Rhohv"][0]


def species(sighh, sigvv, rho):
    """Scattering terms of one species with a given rhohv (real covariance)."""
    return {"sighh": sighh, "sigvv": sigvv, "REdeltaco": rho * np.sqrt(sighh * sigvv) / (4 * np.pi), "IMdeltaco": 0.}


def test_rhohv_single_species():
    np.testing.assert_allclose(rhohv_of_sum([species(2., 1., 0.98)]), 0.98)


def test_rhohv_identical_species_unchanged():
    s = species(2., 1., 0.98)
    np.testing.assert_allclose(rhohv_of_sum([s, s]), 0.98)


def test_rhohv_mixture_decorrelation():
    # two species of equal Zh, Zdr 3 dB and 0 dB, each with rhohv = 1: the mixture is decorrelated
    a, b = species(1., 10 ** -0.3, 1.), species(1., 1., 1.)
    expected = (np.sqrt(10 ** -0.3) + 1.) / np.sqrt(2. * (1. + 10 ** -0.3))
    np.testing.assert_allclose(rhohv_of_sum([a, b]), expected)
    assert rhohv_of_sum([a, b]) < 0.99


def test_rhohv_capped_at_one():
    assert rhohv_of_sum([species(1., 1., 1.0008)]) == 1.
