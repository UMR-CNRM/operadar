"""Tests of the mixed-phase parametrizations (operadar/microphysics/mixed_phase.py).

Run from the repository root: pytest tests/test_mixed_phase.py
"""
import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pytest

try:
    from operadar.microphysics.mixed_phase import (apply_mixed_phase_parametrization,
                                                   compute_liquid_water_fraction)
except ImportError:
    # Environment without operadar's heavy dependencies (epygram, pandas...): load mixed_phase.py alone.
    # Only the two functions tested here are used; they need nothing from the stubbed helper modules.
    for name in ("operadar", "operadar.utils", "operadar.utils.masking", "operadar.utils.make_links"):
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["operadar.utils.masking"].mask_bright_band = None
    sys.modules["operadar.utils.make_links"].link_keys_with_available_hydrometeors = None
    _path = Path(__file__).resolve().parents[1] / "operadar" / "microphysics" / "mixed_phase.py"
    _spec = importlib.util.spec_from_file_location("mixed_phase", _path)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    apply_mixed_phase_parametrization = _mod.apply_mixed_phase_parametrization
    compute_liquid_water_fraction = _mod.compute_liquid_water_fraction

# 4 grid points: [melting layer with graupel and hail, melting layer with graupel only,
#                 hail without rain (not in the melting layer), graupel without rain]
RR = [1e-3, 2e-3, 0.0, 0.0]
GG = [1e-3, 1e-3, 0.0, 3e-3]
HH = [2e-3, 0.0, 4e-3, 0.0]
BB = np.array([True, True, False, False])


def run(parametrization, wet_species):
    contents = {"rr": np.array(RR), "gg": np.array(GG), "hh": np.array(HH)}
    Fw = compute_liquid_water_fraction(contents, BB, wet_species)
    for w in wet_species:  # as compute_mixed_phase does
        contents[w] = np.copy(contents[w[1] * 2])
    return apply_mixed_phase_parametrization(contents, Fw, BB, wet_species, parametrization), Fw


def total(c):
    return sum(c[k] for k in ("rr", "gg", "hh", "wg", "wh") if k in c)


@pytest.mark.parametrize("parametrization", ["Fw_pos", "Fw_posg"])
@pytest.mark.parametrize("wet_species", [["wg"], ["wg", "wh"]], ids=["wg", "wg+wh"])
def test_total_content_is_conserved(parametrization, wet_species):
    """No hydrometeor is counted twice: the total content (rain + graupel + hail, dry and wet) is unchanged."""
    c, _ = run(parametrization, wet_species)
    np.testing.assert_allclose(total(c), np.array(RR) + np.array(GG) + np.array(HH), rtol=1e-12)


@pytest.mark.parametrize("parametrization", ["Fw_pos", "Fw_posg"])
def test_wet_hail_only_in_melting_layer(parametrization):
    """With wet hail, hail is wet in the melting layer and dry elsewhere, never both."""
    c, _ = run(parametrization, ["wg", "wh"])
    assert np.all(c["hh"][BB] == 0)
    assert np.all(c["wh"][BB] >= np.array(HH)[BB])  # all the hail (plus, for Fw_pos, its share of rain)
    assert np.all(c["wh"][~BB] == 0)
    np.testing.assert_array_equal(c["hh"][~BB], np.array(HH)[~BB])
    np.testing.assert_array_equal(c["gg"][~BB], np.array(GG)[~BB])
    assert np.all(c["wg"][~BB] == 0)


def test_fw_pos_shares_rain_between_wet_graupel_and_wet_hail():
    """Fw_pos: the rain of the melting layer is shared between wet graupel and wet hail in proportion to the
    graupel and hail contents, so that each gets the liquid fraction Fw = rr / (rr + gg + hh)."""
    c, Fw = run("Fw_pos", ["wg", "wh"])
    rr, gg, hh = RR[0], GG[0], HH[0]
    assert c["wg"][0] == pytest.approx(gg + rr * gg / (gg + hh))
    assert c["wh"][0] == pytest.approx(hh + rr * hh / (gg + hh))
    for wet, dry in (("wg", gg), ("wh", hh)):
        assert 1 - dry / c[wet][0] == pytest.approx(Fw[0])  # liquid fraction of the wet particle
    assert c["rr"][0] == 0 and c["gg"][0] == 0 and c["hh"][0] == 0
    # without hail, all the rain goes into wet graupel (unchanged behaviour)
    assert c["wg"][1] == pytest.approx(RR[1] + GG[1])


def test_fw_pos_without_wet_hail_unchanged():
    """Fw_pos without wet hail: wet graupel = graupel + rain in the melting layer, dry hail is left as it is."""
    c, _ = run("Fw_pos", ["wg"])
    np.testing.assert_allclose(c["wg"], [RR[0] + GG[0], RR[1] + GG[1], 0, 0])
    np.testing.assert_array_equal(c["hh"], HH)
    assert np.all(c["rr"][BB] == 0) and np.all(c["gg"][BB] == 0)
