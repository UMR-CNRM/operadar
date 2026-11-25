#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests of the NetCDF lookup tables (operadar/read/table_io.py) and of the radar sensitivity threshold
(operadar/radar/sensitivity.py). The modules are loaded directly, without the operadar package (whose
import needs epygram)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


table_io = _load("table_io", "operadar/read/table_io.py")
sensitivity = _load("sensitivity", "operadar/radar/sensitivity.py")

INT_HEADER = ("micro;canting;SIGBETA;DIEL;ARfunc;ARcnst;DSTYfunc;LAM;Frim;ELEVmin;ELEVstep;ELEVmax;Tcmin;Tcstep;"
              "Tcmax;expMmin;expMstep;expMmax;expCCmin;expCCstep;expCCmax;Fwmin;Fwstep;Fwmax")
INT_COLS = ("Tc;ELEV;P3;M;RR;sighh;sigvv;REdeltaco;IMdeltaco;zhh;zdr;kdp;rhohv;Ah;Av;sighhR;sigvvR;REdeltacoR;"
            "IMdeltacoR;zhhR;zdrR;kdpR;rhohvR;AhR;AvR;Mfrac")


def write_int_table(path, two_moment=False, rng=None):
    """Small TmatCoefInt text table: Tc (-10, 0), ELEV (0, 10), Fw (0, 0.5) or N (1, 1e9 m-3), M (3 values)."""
    rng = rng or np.random.default_rng(0)
    p3_values = [1.0, 1e9] if two_moment else [0.0, 0.5]
    lines = [INT_HEADER,
             "ICE3;orient;40.0;LBwetgr;RYwg;0.80;RHOX; 53.20;  1.0;  0.00; 10.00; 10.00;-10.00; 10.00;  0.00;"
             " -3.00;  0.50; -2.00;  0.00;  9.00;  9.00;  0.00;  0.50;  " + ("0.00" if two_moment else "0.50"),
             INT_COLS]
    names = INT_COLS.split(";")[4:]
    rows = []
    for tc in (-10.0, 0.0):
        for el in (0.0, 10.0):
            for p3 in p3_values:
                for m in (1e-3, 10 ** -2.5, 1e-2):
                    vals = [f"{x:.4E}" for x in rng.uniform(1e-3, 10, len(names))]
                    vals[names.index("IMdeltaco")] = "0.1000E-41"  # single-precision subnormal
                    if (tc, el, p3) == (0.0, 10.0, p3_values[0]):
                        vals[names.index("rhohv")] = "  Infinity"  # as written by TmatInt (underflow)
                    rows.append(f"{tc:5.1f};{el:4.1f};{p3:.3E};{m:.4E};" + ";".join(vals))
    Path(path).write_text("\n".join(lines + rows) + "\n")


def test_int_table_netcdf_round_trip(tmp_path):
    """A TmatCoefInt text table converted to NetCDF gives the same values (float32) in the row order of the
    text table, the header parameters as attributes, 'Infinity' read as NaN, and dims (Tc, ELEV, Fw, M)."""
    import xarray as xr
    txt = tmp_path / "TmatCoefInt_ICE3_1M_Cwg"
    write_int_table(txt)
    ref = table_io.read_table_columns(str(txt))          # text (before conversion)
    nc = table_io.text_to_netcdf(str(txt))
    got = table_io.read_table_columns(str(txt))          # now reads <table>.nc
    assert table_io.table_path(str(txt)) == nc
    for c in ref:
        atol = 1e-44 if c == "IMdeltaco" else 0  # subnormal: float32 precision is lower there
        assert np.allclose(got[c], ref[c], rtol=1e-6, atol=atol, equal_nan=True), c
    assert np.isnan(got["rhohv"]).sum() == np.isnan(ref["rhohv"]).sum() > 0
    params = table_io.read_table_params(str(txt))
    assert params["Tcmin"] == -10.0 and params["expMstep"] == 0.5 and params["canting"] == "orient"
    with xr.open_dataset(nc) as ds:
        assert list(ds.zhh.dims) == ["Tc", "ELEV", "Fw", "M"]
        assert ds.sizes == {"Tc": 2, "ELEV": 2, "Fw": 2, "M": 3}
        assert ds.IMdeltaco.dtype == np.float32 and float(ds.IMdeltaco.min()) > 0  # subnormal kept


def test_two_moment_table_has_n_dimension(tmp_path):
    """2-moment tables (and pristine ice) have N instead of Fw as third dimension."""
    import xarray as xr
    txt = tmp_path / "TmatCoefInt_LIMA_2M_Crr"
    write_int_table(txt, two_moment=True)
    nc = table_io.text_to_netcdf(str(txt))
    with xr.open_dataset(nc) as ds:
        assert list(ds.sighh.dims) == ["Tc", "ELEV", "N", "M"]


def test_incomplete_table_is_rejected(tmp_path):
    """A table with missing rows (e.g. a T-matrix that stopped) is not converted."""
    txt = tmp_path / "TmatCoefInt_ICE3_1M_Cwg"
    write_int_table(txt)
    lines = txt.read_text().splitlines()
    txt.write_text("\n".join(lines[:-2]) + "\n")
    with pytest.raises(ValueError):
        table_io.text_to_netcdf(str(txt))


def test_sensitivity_threshold_range_dependent():
    """Zh, Zdr, Kdp and Rhohv are removed where Zh < Zmin_1km + 20 log10(r / 1 km); Ah is kept."""
    zh = np.array([-35.0, -25.0, -15.0, -5.0, np.nan])
    r = np.array([1e3, 1e3, 10e3, 100e3, 1e3])          # Zmin: -30, -30, -10, +10, -30
    d = {"Zh": zh.copy(), "Zdr": np.ones(5), "Kdp": np.ones(5), "Rhohv": np.ones(5), "Ah": np.ones(5)}
    out = sensitivity.apply_sensitivity_threshold(d, band="C", slant_range=r, verbose=False)
    kept = [False, True, False, False, False]
    for var in ("Zh", "Zdr", "Kdp", "Rhohv"):
        assert list(np.isfinite(out[var])) == kept, var
    assert np.all(out["Ah"] == 1)


def test_sensitivity_threshold_constant_and_attenuated():
    """Without range the threshold is constant (band default, W: -50 dBZ); with attenuation the
    attenuated Zh decides."""
    d = {"Zh": np.array([-45.0, -55.0]), "Zh_att": np.array([-60.0, -55.0]), "Rhohv": np.ones(2)}
    out = sensitivity.apply_sensitivity_threshold(d, band="W", verbose=False)
    assert np.all(np.isnan(out["Rhohv"])) and np.all(np.isnan(out["Zh"]))
    d = {"Zh": np.array([-45.0, -55.0]), "Rhohv": np.ones(2)}
    out = sensitivity.apply_sensitivity_threshold(d, band="W", zmin_1km=-60.0, verbose=False)
    assert np.all(np.isfinite(out["Rhohv"]))


def test_vertical_range():
    """Vertical pointing radar: range = |altitude - radar altitude| (zenith and nadir)."""
    r = sensitivity.vertical_range(np.array([0.0, 2000.0, 8000.0]), 3000.0, (3,))
    assert list(r) == [3000.0, 1000.0, 5000.0]
