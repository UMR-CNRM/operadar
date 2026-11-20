"""Tests of the lookup-table generator (tables_generator/ + generate_tables.sh).

Run from the repository root with:
    pytest tests/tables_generator -v
Fortran tests need gfortran and a LAPACK/refBLAS build: set LAPACK_DIR to the folder
holding liblapack.a and librefblas.a (default: tables_generator/src/libFortran/lapack).
They are skipped otherwise. Parameter-file and shell-script tests only need Python.

Tests marked xfail(strict=True) document bugs found in the code review: they fail
today and will turn into an XPASS error once the bug is fixed, as a reminder to
remove the marker.
"""
from __future__ import annotations

import math
import re
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from harness import (LAM_MM, MICRO_KEYS, PARAM, PARAM_LAYOUT, REPO, SRC, Sandbox, build, default_param,
                     format_param, lapack_dir, read_coefdiff, read_coefdiff_header, read_coefint)

RHOI = 916.0  # solid ice density used in Tmatrix.f (kg/m3)

needs_fortran = pytest.mark.skipif(shutil.which("gfortran") is None or lapack_dir() is None,
                                   reason="gfortran or LAPACK (LAPACK_DIR) not available")


# ---------------------------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def binaries(tmp_path_factory):
    if shutil.which("gfortran") is None or lapack_dir() is None:
        pytest.skip("gfortran or LAPACK (LAPACK_DIR) not available")
    return build(tmp_path_factory.mktemp("build"))


@pytest.fixture
def sandbox(tmp_path, binaries):
    return Sandbox(tmp_path / "gen", *binaries)


def tmat_rows(sandbox, **param):
    return read_coefdiff(sandbox.run_tmat(default_param(**param)))


def resu_rows(sandbox, typeh, band):
    """D, RHOP ... from the TmatResu side file written by Tmat."""
    lines = (sandbox.root / "tables" / typeh / f"TmatResu_{band}{typeh}").read_text().splitlines()[1:]
    return [dict(D=float(l.split()[0]), RHOP=float(l.split()[4])) for l in lines]


# ---------------------------------------------------------------------------------------------
# 1. Parameter files (pure Python)
# ---------------------------------------------------------------------------------------------
PARAM_FILES = sorted(p for p in PARAM.glob("TmatParam_*") if not p.name.endswith(".ref"))

# Parameter files with known defects (see review); strict xfail so a fix is noticed.
KNOWN_BAD: dict[str, str] = {}
PARAM_CASES = [pytest.param(p, id=p.name, marks=pytest.mark.xfail(strict=True, reason=KNOWN_BAD[p.name]))
               if p.name in KNOWN_BAD else pytest.param(p, id=p.name) for p in PARAM_FILES]


def parse_param(path: Path) -> dict:
    lines = path.read_text().splitlines()
    out = {}
    for (key, fmt), line in zip(PARAM_LAYOUT, lines):
        k, _, v = line.partition(":")
        assert k == key, f"{path.name}: expected key '{key}' got '{k}' (Tmatrix.f reads lines by position)"
        width = int(re.match(r"[AF](\d+)", fmt).group(1))
        assert len(v) <= width, f"{path.name}: '{line}' does not fit the Fortran edit descriptor {fmt}"
        out[k] = v if fmt.startswith("A") else float(v)
    for key, line in zip(MICRO_KEYS, lines[len(PARAM_LAYOUT):]):
        k, _, v = line.partition(":")
        assert k == key and v in ("1M", "2M"), f"{path.name}: bad microphysics line '{line}'"
        out[k] = v
    return out


@pytest.mark.parametrize("path", PARAM_CASES)
def test_param_file_layout(path):
    p = parse_param(path)
    for name in ("ELEV", "Tc", "expD", "Fw"):
        lo, hi, step = p[name + "min"], p[name + "max"], p[name + "step"]
        assert hi >= lo, f"{name}max < {name}min"
        assert step > 0
        n = (hi - lo) / step
        assert abs(n - round(n)) < 1e-6, f"{name} range is not a whole number of steps"
    m = re.match(r"TmatParam_(L|S|C|X|Ku|Ka|K|W)(\w\w)(_|$)", path.name)
    assert m, "unexpected file name"
    assert p["band"] == m.group(1)
    assert p["type"] == m.group(2)
    assert abs(p["LAM"] - LAM_MM[p["band"]]) < 0.06, f"LAM={p['LAM']} inconsistent with band {p['band']}"
    assert p["canting"] in ("nocant", "analyt")
    assert p["DIEL"] in ("Liebe91", "RY19dry", "LBwetgr", "MGwMA08")
    assert p["ARfunc"] in ("BR02", "CNST", "AUds", "RYdg", "RYwg")


VERTICAL_FILES = [p for p in PARAM_FILES if p.name.endswith("_vertical")]


@pytest.mark.parametrize("path", VERTICAL_FILES, ids=lambda p: p.name)
def test_vertical_configs_include_90deg(path):
    p = parse_param(path)
    assert p["ELEVmax"] == 90.0


def test_mass_diameter_constants_consistent_between_schemes():
    """Tmat (scattering per diameter) always reads a,b from ICE3_constants.txt, whereas TmatInt
    integrates with <MICRO>_constants.txt. A per-diameter table is only reusable across schemes
    if a,b are identical for each species."""
    def ab(f):
        rows = [l.split() for l in (PARAM / f).read_text().splitlines()[1:] if l.strip()]
        return {r[0]: (float(r[1]), float(r[2])) for r in rows}
    ref = ab("ICE3_constants.txt")
    for f in ("LIMA_constants.txt", "ICJW_constants.txt", "LIMC_constants.txt"):
        for sp, val in ab(f).items():
            if sp in ref:
                assert val == ref[sp], f"{f}: a,b of '{sp}' differ from ICE3 ({val} vs {ref[sp]})"


# ---------------------------------------------------------------------------------------------
# 2. Tmat: scattering per diameter
# ---------------------------------------------------------------------------------------------
@needs_fortran
def test_tmat_header_roundtrip(sandbox):
    path = sandbox.run_tmat(default_param(Tcmin=-5.0, Tcmax=5.0, Tcstep=5.0, ELEVmax=10.0, ELEVstep=10.0))
    h = read_coefdiff_header(path)
    assert float(h["Tcmin"]) == -5.0 and float(h["Tcmax"]) == 5.0 and float(h["ELEVmax"]) == 10.0
    rows = read_coefdiff(path)
    assert len(rows) == 3 * 2 * 20  # nTc * nELEV * nD (nFw = 1)
    # loop order is Tc, ELEV, Fw, D (TmatInt relies on it to compute kTmat)
    assert [r["Tc"] for r in rows[::40]] == [-5.0, 0.0, 5.0]
    assert [r["ELEV"] for r in rows[:40:20]] == [0.0, 10.0]


@needs_fortran
@pytest.mark.parametrize("band", ["S", "C", "X"])
def test_tmatrix_matches_rayleigh_for_small_drops(sandbox, band):
    for r in tmat_rows(sandbox, band=band, expDmin=-1.0, expDmax=-0.3, expDstep=0.1):
        assert abs(r["zhh"] - r["zhhR"]) < 0.05, r
        assert abs(r["zdr"] - r["zdrR"]) < 0.01, r


@needs_fortran
def test_spheres_have_no_polarimetric_signature(sandbox):
    for r in tmat_rows(sandbox, band="C", ARfunc="CNST", ARcnst=1.0, expDmin=-1.0, expDmax=0.8):
        assert abs(r["zdr"]) < 1e-3 and abs(r["rhohv"] - 1) < 1e-6 and abs(r["kdp"]) < 1e-8, r


@needs_fortran
def test_rain_low_elevation_signature(sandbox):
    rows = tmat_rows(sandbox, band="S", expDmin=0.0, expDmax=0.8, expDstep=0.1)
    zdr = [r["zdr"] for r in rows]
    assert all(z > 0 for z in zdr) and zdr == sorted(zdr), zdr
    assert all(r["kdp"] > 0 for r in rows)
    assert all(r["Ah"] > 0 and r["Ah"] >= r["Av"] for r in rows)


@needs_fortran
def test_rain_vertical_incidence_is_isotropic(sandbox):
    rows = tmat_rows(sandbox, band="C", ELEVmin=90.0, ELEVmax=90.0, expDmin=0.0, expDmax=0.8)
    for r in rows:
        assert abs(r["zdr"]) < 1e-3 and abs(r["kdp"]) < 1e-6, r


@needs_fortran
def test_analytical_canting_with_zero_width_equals_no_canting(sandbox):
    common = dict(band="C", expDmin=0.0, expDmax=0.8)
    a = tmat_rows(sandbox, canting="analyt", sigbeta=0.0, **common)
    b = tmat_rows(sandbox, canting="nocant", **common)
    for ra, rb in zip(a, b):
        for k in ("sighh", "sigvv", "kdp", "Ah", "Av", "REdeltaco"):
            assert ra[k] == pytest.approx(rb[k], rel=1e-4, abs=1e-30), k


@needs_fortran
def test_canting_reduces_zdr_and_kdp(sandbox):
    common = dict(band="C", expDmin=0.3, expDmax=0.8)
    a = tmat_rows(sandbox, canting="analyt", sigbeta=20.0, **common)
    b = tmat_rows(sandbox, canting="nocant", **common)
    for ra, rb in zip(a, b):
        assert 0 < ra["zdr"] < rb["zdr"] and 0 < ra["kdp"] < rb["kdp"]


@needs_fortran
def test_rayleigh_sphere_reflectivity_is_d6(sandbox):
    """For a small water sphere Z = D**6 whatever the band: the radar constant must use |K|^2."""
    for r in tmat_rows(sandbox, band="W", Tcmin=0.0, Tcmax=0.0, ARfunc="CNST", ARcnst=1.0,
                       expDmin=-2.0, expDmax=-1.5, expDstep=0.1):
        assert r["zhhR"] == pytest.approx(60 * math.log10(r["D"]), abs=0.02)


@needs_fortran
def test_wet_graupel_canting_does_not_depend_on_loop_history(sandbox):
    common = dict(type="wg", band="C", canting="analyt", sigbeta=40.0, DIEL="LBwetgr", ARfunc="RYwg",
                  ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, Fwmin=0.5, Fwmax=0.5, Fwstep=0.05, expDstep=0.1)
    single = tmat_rows(sandbox, expDmin=0.5, expDmax=0.5, **common)[0]
    many = tmat_rows(sandbox, expDmin=0.0, expDmax=0.5, **common)[-1]
    assert single["D"] == many["D"]
    assert single["zdr"] == pytest.approx(many["zdr"], rel=1e-4)


@needs_fortran
def test_wet_graupel_without_water_equals_dry_graupel(sandbox):
    common = dict(band="C", ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, expDmin=-1.0, expDmax=1.0, expDstep=0.1)
    wet = tmat_rows(sandbox, type="wg", DIEL="LBwetgr", ARfunc="RYwg", **common)
    dry = tmat_rows(sandbox, type="gg", DIEL="RY19dry", ARfunc="RYdg", **common)
    for rw, rd in zip(wet, dry):
        assert rw["zhh"] == pytest.approx(rd["zhh"], abs=0.01)
        assert rw["zdr"] == pytest.approx(rd["zdr"], abs=0.001)


@needs_fortran
@pytest.mark.parametrize("typeh,arfunc,arcnst,dmin", [("ss", "AUds", 0.75, -1.0), ("ii", "CNST", 0.2, -2.0)])
def test_particle_density_not_above_solid_ice(sandbox, typeh, arfunc, arcnst, dmin):
    sandbox.run_tmat(default_param(type=typeh, ARfunc=arfunc, ARcnst=arcnst, DIEL="RY19dry", Tcmin=-10.0,
                                   Tcmax=-10.0, expDmin=dmin, expDmax=0.0, expDstep=0.5))
    for r in resu_rows(sandbox, typeh, "S"):
        assert r["RHOP"] <= RHOI, r


@needs_fortran
def test_tmat_fails_without_type_subfolder(tmp_path, binaries):
    """Tmat writes to ../tables/<type>/ and crashes if that folder is missing; generate_tables.sh
    must therefore create it (see test_default_mode_generates_all_tables)."""
    sb = Sandbox(tmp_path / "gen", *binaries)
    (sb.root / "param" / "tmp_config").write_text(format_param(default_param()))
    r = subprocess.run(["./Tmat"], cwd=sb.root / "src", capture_output=True, text=True, timeout=60)
    assert not (sb.root / "tables" / "rr").exists()
    assert "Cannot open file" in r.stderr or r.returncode != 0


# ---------------------------------------------------------------------------------------------
# 3. TmatInt: integration over the PSD
# ---------------------------------------------------------------------------------------------
ICE3_RAIN = dict(a=524.0, b=3.0, C=8e6, X=-1.0, c=842.0, d=0.8)


def ice3_rain_moments(M):
    """Analytic Z (mm6 m-3) and rain rate (mm/h) for the ICE3 exponential rain PSD."""
    p = ICE3_RAIN
    lam = (M / (p["a"] * p["C"] * math.gamma(4))) ** (1 / (p["X"] - p["b"]))
    ntot = p["C"] * lam ** p["X"]
    z = ntot * lam * math.gamma(7) / lam ** 7 * 1e18
    rr = 3600 * p["a"] * p["c"] * ntot * lam * math.gamma(4 + p["d"]) / lam ** (4 + p["d"])
    return 10 * math.log10(z), rr


@pytest.fixture(scope="module")
def sphere_rain_table(tmp_path_factory, binaries):
    sb = Sandbox(tmp_path_factory.mktemp("int") / "gen", *binaries)
    diff = sb.run_tmat(default_param(band="S", ARfunc="CNST", ARcnst=1.0, expDmin=-2.0, expDmax=0.9,
                                     expDstep=0.02))
    return read_coefint(sb.run_tmatint(diff, "t", "rr", "S"))


@needs_fortran
def test_coefint_header(sphere_rain_table):
    header, rows = sphere_rain_table
    assert header["micro"] == "ICE3" and float(header["LAM"]) == 106.2
    assert len(rows) == round((-1.5 + 7.0) / 0.01) + 1
    assert rows[0]["M"] == pytest.approx(1e-7, rel=1e-3) and rows[-1]["M"] == pytest.approx(10 ** -1.5, rel=1e-3)


@needs_fortran
def test_integrated_rayleigh_reflectivity_matches_analytic(sphere_rain_table):
    _, rows = sphere_rain_table
    for r in rows:
        if r["M"] <= 1e-3:  # above, the 8 mm cut-off of the diameter range truncates the D^6 moment
            z, _ = ice3_rain_moments(r["M"])
            assert r["zhhR"] == pytest.approx(z, abs=0.05), r["M"]


@needs_fortran
def test_integrated_rain_rate_matches_analytic(sphere_rain_table):
    _, rows = sphere_rain_table
    for r in rows:
        if 1e-6 <= r["M"] <= 1e-3:
            _, rr = ice3_rain_moments(r["M"])
            assert r["RR"] == pytest.approx(rr, rel=0.03), r["M"]


@needs_fortran
def test_integrated_rain_variables_increase_with_content(sandbox):
    diff = sandbox.run_tmat(default_param(band="C", expDmin=-2.0, expDmax=0.9, expDstep=0.02))
    _, rows = read_coefint(sandbox.run_tmatint(diff, "t", "rr", "C"))
    rows = [r for r in rows if r["M"] <= 1e-2]
    for k in ("zhh", "zdr", "kdp", "Ah"):
        v = [r[k] for r in rows]
        assert all(b >= a for a, b in zip(v, v[1:])), k


@needs_fortran
def test_two_moment_ice_reflectivity_trends(sandbox):
    diff = sandbox.run_tmat(default_param(type="ii", band="C", DIEL="RY19dry", ARfunc="CNST", ARcnst=1.0,
                                          Tcmin=-20.0, Tcmax=-20.0, expDmin=-2.0, expDmax=0.5, expDstep=0.05))
    _, rows = read_coefint(sandbox.run_tmatint(diff, "t", "ii", "C"))
    by_n = {}
    for r in rows:
        by_n.setdefault(r["P3"], []).append(r)
    assert len(by_n) == 91  # expCC 0..9 step 0.1
    n_lo, n_hi = sorted(by_n)[40], sorted(by_n)[60]  # 1e4 and 1e6 m-3
    for lo, hi in zip(by_n[n_lo], by_n[n_hi]):
        assert lo["M"] == hi["M"]
        if 1e-6 <= lo["M"] <= 1e-4:
            assert lo["zhhR"] > hi["zhhR"]  # same mass, fewer (bigger) crystals -> higher Z


@pytest.fixture(scope="module")
def graupel_tables(tmp_path_factory, binaries):
    sb = Sandbox(tmp_path_factory.mktemp("gr") / "gen", *binaries)
    common = dict(band="C", ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, expDmin=-1.0, expDmax=1.4, expDstep=0.02)
    wg = sb.run_tmat(default_param(type="wg", DIEL="LBwetgr", ARfunc="RYwg", Fwmin=0.0, Fwmax=0.5,
                                   Fwstep=0.5, **common))
    gg = sb.run_tmat(default_param(type="gg", DIEL="RY19dry", ARfunc="RYdg", **common))
    _, rwg = read_coefint(sb.run_tmatint(wg, "t", "wg", "C"))
    _, rgg = read_coefint(sb.run_tmatint(gg, "t", "gg", "C"))
    return rwg, rgg


@needs_fortran
def test_wet_graupel_integration_without_water_equals_dry_graupel(graupel_tables):
    rwg, rgg = graupel_tables
    dry_part = [r for r in rwg if r["P3"] == 0.0]
    wet_part = [r for r in rwg if r["P3"] == 0.5]
    assert len(dry_part) == len(rgg) == len(wet_part)
    for a, b in zip(dry_part, rgg):
        if a["M"] >= 1e-4:
            assert a["zhh"] == pytest.approx(b["zhh"], abs=0.02), a["M"]
    for a, b in zip(wet_part, rgg):
        if 1e-5 <= a["M"] <= 1e-3:
            assert a["zhh"] > b["zhh"]  # melting graupel is brighter than dry graupel


@needs_fortran
@pytest.mark.xfail(strict=True, reason="review: the M/Mint renormalisation cancels for non-melting species "
                                       "(PSD keeps the truncated mass) but not for wg/ws (PSD scaled to M)")
def test_wet_graupel_and_graupel_same_normalisation_at_low_content(graupel_tables):
    rwg, rgg = graupel_tables
    for a, b in zip([r for r in rwg if r["P3"] == 0.0], rgg):
        assert a["zhh"] == pytest.approx(b["zhh"], abs=0.02), a["M"]


# ---------------------------------------------------------------------------------------------
# 4. PSD subroutine of TmatrixIntegration.f90, compiled with a small driver
# ---------------------------------------------------------------------------------------------
PSD_DRIVER = """
program t
  implicit none
  real*8 :: lamb, n
  ! call 1 : another content, to leave a different value in lamb
  call PSD(1e-6, 1e-3, 273.15, 'ICE3', 0.0, 'rr', 524., 3., 1., 1., 8e6, -1., 15., 1, lamb, n)
  ! call 2 : the one we check
  call PSD(%(M)s, %(D)s, 273.15, 'ICE3', 0.0, 'rr', 524., 3., 1., 1., 8e6, -1., 15., 1, lamb, n)
  write(*,'(2ES24.15)') lamb, n
end program
"""


def run_psd(tmp_path, M, D):
    src = (SRC / "TmatrixIntegration.f90").read_text()
    psd = src[src.index("SUBROUTINE PSD("):]
    (tmp_path / "psd.f90").write_text(PSD_DRIVER % dict(M=M, D=D) + "\n" + psd)
    subprocess.run(["gfortran", "-O0", "-w", "psd.f90", "-o", "psd"], cwd=tmp_path, check=True,
                   capture_output=True)
    out = subprocess.run(["./psd"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout
    return map(float, out.split())


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
@pytest.mark.xfail(strict=True, reason="review: in PSD (ICE3, 1-moment) No=c*lamb**x is computed before lamb")
def test_psd_ice3_rain_marshall_palmer(tmp_path):
    lamb, n = run_psd(tmp_path, "1e-3", "1e-3")
    p = ICE3_RAIN
    lam_ref = (1e-3 / (p["a"] * p["C"] * 6)) ** (1 / (p["X"] - p["b"]))
    assert lamb == pytest.approx(lam_ref, rel=1e-5)
    assert n == pytest.approx(p["C"] * math.exp(-lam_ref * 1e-3), rel=1e-4)  # N(D) = N0 exp(-lam D)


# ---------------------------------------------------------------------------------------------
# 5. generate_tables.sh, run with fake Tmat / TmatInt that mimic the file contracts
# ---------------------------------------------------------------------------------------------
FAKE_TMAT = """#!/bin/bash
# mimics Tmat: reads ../param/tmp_config, writes ../tables/<type>/TmatCoefDiff_<band><type>
echo "$PWD" >> ../calls_tmat
t=$(sed -n 's/^type://p' ../param/tmp_config); b=$(sed -n 's/^band://p' ../param/tmp_config)
[ -d ../tables/$t ] || { echo "Cannot open file '../tables/$t/TmatCoefDiff_$b$t'" >&2; exit 2; }
echo fake > ../tables/$t/TmatCoefDiff_$b$t
"""
FAKE_TMATINT = """#!/bin/bash
# mimics TmatInt <exec_dir> <repo> <type> <band> <micro> <moment>
echo "$@" >> $1/../calls_tmatint
[ -f $1/../tables/$2/TmatCoefDiff_$4$3 ] || { echo "missing CoefDiff" >&2; exit 2; }
echo fake > $1/../tables/$2/TmatCoefInt_$5_$6_$4$3
"""


@pytest.fixture
def script_root(tmp_path):
    root = tmp_path / "repo"
    gen = root / "tables_generator"
    (gen / "src").mkdir(parents=True)
    (gen / "tables").mkdir()
    shutil.copytree(PARAM, gen / "param")
    shutil.copy(REPO / "generate_tables.sh", root)
    for name, body in (("Tmat", FAKE_TMAT), ("TmatInt", FAKE_TMATINT)):
        f = gen / "src" / name
        f.write_text(body)
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return root


def run_script(root, *args):
    return subprocess.run(["bash", "generate_tables.sh", *args], cwd=root, capture_output=True, text=True,
                          timeout=60)


HYDRO_LIST = ["cs", "cl", "rr", "ii", "ss", "gg", "hh", "wg"]


def test_script_rejects_unknown_band(script_root):
    r = run_script(script_root, "--default", "--band", "Q")
    assert r.returncode != 0 and "invalid band" in r.stdout


@pytest.mark.xfail(strict=True, reason="review: generate_tables.sh never creates tables/<type>/, so Tmat "
                                       "cannot write its output")
def test_default_mode_generates_all_tables(script_root):
    run_script(script_root, "--default", "--band", "C")
    out = script_root / "tables_generator" / "tables" / "default"
    missing = [h for h in HYDRO_LIST if not (out / f"TmatCoefInt_ICE3_1M_C{h}").exists()]
    assert not missing, f"no table produced for {missing}"


def test_default_mode_works_when_type_folders_exist(script_root):
    for h in HYDRO_LIST:
        (script_root / "tables_generator" / "tables" / h).mkdir()
    run_script(script_root, "--default", "--band", "C")
    out = script_root / "tables_generator" / "tables" / "default"
    for h in HYDRO_LIST:
        assert (out / f"TmatCoefDiff_C{h}").exists() and (out / f"TmatCoefInt_ICE3_1M_C{h}").exists(), h


def test_default_mode_reuses_existing_tables(script_root):
    out = script_root / "tables_generator" / "tables" / "default"
    out.mkdir()
    for h in HYDRO_LIST:
        (out / f"TmatCoefDiff_C{h}").write_text("x")
        (out / f"TmatCoefInt_ICE3_1M_C{h}").write_text("x")
    run_script(script_root, "--default", "--band", "C")
    assert not (script_root / "tables_generator" / "calls_tmat").exists()
    assert not (script_root / "tables_generator" / "calls_tmatint").exists()


@pytest.mark.xfail(strict=True, reason="review: MODE is set to 'newConf' but the final case tests "
                                       "'--newConf', so --newConf falls into edit mode")
def test_newconf_mode_runs_generator(script_root):
    for h in HYDRO_LIST:
        (script_root / "tables_generator" / "tables" / h).mkdir()
    r = run_script(script_root, "--newConf", "myconf", "--band", "C")
    assert "toto_C.txt" not in r.stdout
    assert (script_root / "tables_generator" / "calls_tmat").exists(), r.stdout


@pytest.mark.xfail(strict=True, reason="review: valid_arf uses ARF_LIST but the array is ARfunc_LIST")
def test_arf_option_accepts_valid_value(script_root):
    r = run_script(script_root, "--hydro", "rr", "--band", "C", "--arf", "CNST")
    assert "invalid axis ratio" not in r.stdout
