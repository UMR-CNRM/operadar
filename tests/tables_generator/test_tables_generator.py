"""Tests of the lookup-table generator (tables_generator/ + generate_tables.sh).

Run from the repository root with:
    LAPACK_DIR=<folder with liblapack.a and librefblas.a> pytest tests/tables_generator
then open tests/tables_generator/report/test_report.html (see README.md in this folder).

Fortran tests need gfortran and a LAPACK/refBLAS build (LAPACK_DIR, default
tables_generator/src/libFortran/lapack); they are skipped otherwise.

Tests marked xfail(strict=True) document known bugs: they fail until the bug is fixed,
then are reported as "BUG FIXED?" as a reminder to remove the marker. (None at the moment.)
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
# Fixtures and helpers
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
    """D, RHOP from the TmatResu side file written by Tmat."""
    lines = (sandbox.root / "tables" / typeh / f"TmatResu_{band}{typeh}").read_text().splitlines()[1:]
    return [dict(D=float(l.split()[0]), RHOP=float(l.split()[4])) for l in lines]


def col(rows, key):
    return [r[key] for r in rows]


# =============================================================================================
PARAM_FILES = sorted(p for p in PARAM.glob("TmatParam_*") if not p.name.endswith(".ref"))
KNOWN_BAD: dict[str, str] = {}  # file name -> reason (strict xfail) for known-defective parameter files
PARAM_CASES = [pytest.param(p, id=p.name, marks=pytest.mark.xfail(strict=True, reason=KNOWN_BAD[p.name]))
               if p.name in KNOWN_BAD else pytest.param(p, id=p.name) for p in PARAM_FILES]
VERTICAL_FILES = [p for p in PARAM_FILES if p.name.endswith("_vertical")]


class TestParameterFiles:
    """1. Parameter files (param/TmatParam_*, *_constants.txt)

    Tmat reads its parameter file line by line with fixed Fortran formats, so a missing
    line or a value that does not fit silently shifts or truncates everything after it.
    These checks only need Python.
    """

    @staticmethod
    def parse(path: Path) -> dict:
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
    def test_param_file_layout(self, path):
        """Each TmatParam file has the 26 lines in the order Tmat reads them, every value fits its Fortran
        format, min/max/step ranges are consistent, and type/band/wavelength match the file name.

        Wavelengths expected per band (mm): L 190.4, S 106.2, C 53.2, X 31.9, Ku 20.0, K 15.0, Ka 8.4, W 3.15."""
        p = self.parse(path)
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

    @pytest.mark.parametrize("path", VERTICAL_FILES, ids=lambda p: p.name)
    def test_vertical_configs_include_90deg(self, path):
        """The '*_vertical' configurations (for vertically pointing radars) contain the 90 degree elevation."""
        assert self.parse(path)["ELEVmax"] == 90.0

    @pytest.mark.parametrize("path", [p for p in PARAM_FILES if re.match(r"TmatParam_[A-Za-z]+ss", p.name)],
                             ids=lambda p: p.name)
    def test_snow_tables_start_at_10_microns(self, path):
        """Snow parameter files start at D = 0.01 mm (expDmin = -2): the ICE3 snow distribution at low content
        is made of very small particles, which a 0.1 mm lower limit misses (Zh of -290 dBZ at 1 mg/m3)."""
        assert self.parse(path)["expDmin"] <= -2.0

    def test_mass_diameter_constants_consistent_between_schemes(self, report):
        """The mass-diameter coefficients a, b of each species are identical in all *_constants.txt files.

        Tmat (scattering per diameter) always reads a, b from ICE3_constants.txt whereas TmatInt integrates
        with the constants of the chosen scheme, so a per-diameter table is only valid for every scheme
        if a and b are the same."""
        def ab(f):
            rows = [l.split() for l in (PARAM / f).read_text().splitlines()[1:] if l.strip()]
            return {r[0]: (float(r[1]), float(r[2])) for r in rows}
        files = ["ICE3_constants.txt", "LIMA_constants.txt", "ICJW_constants.txt", "LIMC_constants.txt"]
        tabs = {f: ab(f) for f in files}
        species = list(tabs["ICE3_constants.txt"])
        report.table("a [kg m^-b], b per species and scheme", ["species"] + [f.split("_")[0] for f in files],
                     [[sp] + ["%g, %g" % tabs[f][sp] if sp in tabs[f] else "-" for f in files] for sp in species])
        ref = tabs["ICE3_constants.txt"]
        for f in files[1:]:
            for sp, val in tabs[f].items():
                if sp in ref:
                    assert val == ref[sp], f"{f}: a,b of '{sp}' differ from ICE3 ({val} vs {ref[sp]})"


# =============================================================================================
@needs_fortran
class TestTmat:
    """2. Tmat: scattering coefficients for one particle (TmatCoefDiff tables)

    Tmat is run on small configurations in a temporary folder. Its results are compared with
    physical limits that are known exactly: Rayleigh scattering for small particles, spheres,
    vertical incidence, zero canting.
    """

    def test_tmat_header_roundtrip(self, sandbox):
        """The TmatCoefDiff file has the expected header and rows ordered Tc > ELEV > Fw > D,
        which is the order TmatInt relies on to find a row."""
        path = sandbox.run_tmat(default_param(Tcmin=-5.0, Tcmax=5.0, Tcstep=5.0, ELEVmax=10.0, ELEVstep=10.0))
        h = read_coefdiff_header(path)
        assert float(h["Tcmin"]) == -5.0 and float(h["Tcmax"]) == 5.0 and float(h["ELEVmax"]) == 10.0
        rows = read_coefdiff(path)
        assert len(rows) == 3 * 2 * 20  # nTc * nELEV * nD (nFw = 1)
        assert [r["Tc"] for r in rows[::40]] == [-5.0, 0.0, 5.0]
        assert [r["ELEV"] for r in rows[:40:20]] == [0.0, 10.0]

    @pytest.mark.parametrize("band", ["S", "C", "X"])
    def test_tmatrix_matches_rayleigh_for_small_drops(self, sandbox, report, band):
        """For raindrops much smaller than the wavelength (D < 0.5 mm), T-matrix and Rayleigh give the same
        Zh (within 0.05 dB) and Zdr (within 0.01 dB). For bigger drops they diverge (Mie effects), mostly at X band."""
        rows = tmat_rows(sandbox, band=band, expDmin=-1.0, expDmax=0.9, expDstep=0.05)
        d = col(rows, "D")
        with report.figure(f"Rain, {band} band, 10 C, elevation 0: T-matrix vs Rayleigh (tested for D < 0.5 mm)",
                           ncols=2) as ax:
            ax[0].semilogx(d, col(rows, "zhh"), label="T-matrix")
            ax[0].semilogx(d, col(rows, "zhhR"), "--", label="Rayleigh")
            ax[0].set(xlabel="D (mm)", ylabel="Zh of one drop per m3 (dBZ)")
            ax[1].semilogx(d, col(rows, "zdr"), label="T-matrix")
            ax[1].semilogx(d, col(rows, "zdrR"), "--", label="Rayleigh")
            for a in ax:
                a.axvspan(min(d), 0.5, color="green", alpha=0.08, label="tested range")
            ax[1].set(xlabel="D (mm)", ylabel="Zdr (dB)")
        for r in rows:
            if r["D"] < 0.5:
                assert abs(r["zhh"] - r["zhhR"]) < 0.05, r
                assert abs(r["zdr"] - r["zdrR"]) < 0.01, r

    def test_spheres_have_no_polarimetric_signature(self, sandbox, report):
        """Spherical drops (axis ratio 1) give Zdr = 0, rhohv = 1 and Kdp = 0 at all sizes."""
        rows = tmat_rows(sandbox, band="C", ARfunc="CNST", ARcnst=1.0, expDmin=-1.0, expDmax=0.8)
        with report.figure("Spheres at C band", ncols=2) as ax:
            ax[0].semilogx(col(rows, "D"), col(rows, "zdr"), "o-")
            ax[0].set(xlabel="D (mm)", ylabel="Zdr (dB)", ylim=(-0.1, 0.1))
            ax[1].semilogx(col(rows, "D"), col(rows, "rhohv"), "o-")
            ax[1].set(xlabel="D (mm)", ylabel="rhohv", ylim=(0.99, 1.01))
        for r in rows:
            assert abs(r["zdr"]) < 1e-3 and abs(r["rhohv"] - 1) < 1e-6 and abs(r["kdp"]) < 1e-8, r

    def test_rain_signature_low_and_vertical_elevation(self, sandbox, report):
        """S band: oblate raindrops seen horizontally give Zdr > 0 increasing with D, Kdp > 0 and Ah >= Av.
        Seen from below (elevation 90 deg) they look circular: Zdr = 0 and Kdp = 0."""
        low = tmat_rows(sandbox, band="S", expDmin=-0.5, expDmax=0.85, expDstep=0.05)
        vert = tmat_rows(sandbox, band="S", ELEVmin=90.0, ELEVmax=90.0, expDmin=-0.5, expDmax=0.85, expDstep=0.05)
        with report.figure("Rain at S band, elevation 0 vs 90 deg", ncols=2) as ax:
            for rows, lab in ((low, "elev 0"), (vert, "elev 90")):
                ax[0].plot(col(rows, "D"), col(rows, "zdr"), label=lab)
                ax[1].plot(col(rows, "D"), col(rows, "kdp"), label=lab)
            ax[0].set(xlabel="D (mm)", ylabel="Zdr (dB)")
            ax[1].set(xlabel="D (mm)", ylabel="Kdp of one drop per m3 (deg/km)")
        big = [r for r in low if r["D"] >= 1]
        zdr = col(big, "zdr")
        assert all(z > 0 for z in zdr) and zdr == sorted(zdr), zdr
        assert all(r["kdp"] > 0 and r["Ah"] >= r["Av"] for r in big)
        for r in vert:
            assert abs(r["zdr"]) < 1e-3 and abs(r["kdp"]) < 1e-6, r

    def test_canting(self, sandbox, report):
        """Canting option 'analyt': with a zero canting width it equals 'nocant' exactly, and a wider
        canting distribution (20 deg) reduces Zdr and Kdp."""
        common = dict(band="C", expDmin=0.0, expDmax=0.8, expDstep=0.05)
        nocant = tmat_rows(sandbox, canting="nocant", **common)
        zero = tmat_rows(sandbox, canting="analyt", sigbeta=0.0, **common)
        wide = tmat_rows(sandbox, canting="analyt", sigbeta=20.0, **common)
        with report.figure("Rain at C band: effect of the canting angle distribution width", ncols=2) as ax:
            for rows, lab, st in ((nocant, "nocant", "-"), (zero, "analyt, 0 deg", "--"), (wide, "analyt, 20 deg", "-")):
                ax[0].plot(col(rows, "D"), col(rows, "zdr"), st, label=lab)
                ax[1].plot(col(rows, "D"), col(rows, "kdp"), st, label=lab)
            ax[0].set(xlabel="D (mm)", ylabel="Zdr (dB)")
            ax[1].set(xlabel="D (mm)", ylabel="Kdp (deg/km)")
        for ra, rb in zip(zero, nocant):
            for k in ("sighh", "sigvv", "kdp", "Ah", "Av", "REdeltaco"):
                assert ra[k] == pytest.approx(rb[k], rel=1e-4, abs=1e-30), k
        for rw, rb in zip(wide, nocant):
            assert 0 < rw["zdr"] < rb["zdr"] and 0 < rw["kdp"] < rb["kdp"]

    def test_rayleigh_sphere_reflectivity_is_d6(self, sandbox, report):
        """For a small water sphere, Z = D^6 at every band. This checks the radar constant uses |K|^2
        (the old code used Re(K^2): +0.4 dB at W band)."""
        res = {}
        for band in ("S", "C", "X", "Ka", "W"):
            rows = tmat_rows(sandbox, band=band, Tcmin=0.0, Tcmax=0.0, ARfunc="CNST", ARcnst=1.0,
                             expDmin=-2.0, expDmax=-1.5, expDstep=0.1)
            res[band] = [r["zhhR"] - 60 * math.log10(r["D"]) for r in rows]
        with report.figure("Water sphere D = 0.01 mm at 0 C: Zh - 10 log10(D^6) per band (should be 0)") as ax:
            ax.bar(list(res), [v[0] for v in res.values()])
            ax.axhline(0, color="k", lw=0.8)
            ax.set(ylabel="error (dB)", ylim=(-0.5, 0.5))
        for band, errs in res.items():
            assert max(abs(x) for x in errs) < 0.02, band

    def test_wet_graupel_canting_does_not_depend_on_loop_history(self, sandbox, report):
        """Wet graupel with analytic canting: computing a diameter alone or inside a whole table gives
        the same result (the old code overwrote the canting width at each diameter)."""
        common = dict(type="wg", band="C", canting="analyt", sigbeta=40.0, DIEL="LBwetgr", ARfunc="RYwg",
                      ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, Fwmin=0.5, Fwmax=0.5, Fwstep=0.05, expDstep=0.1)
        table = tmat_rows(sandbox, expDmin=0.0, expDmax=0.5, **common)
        alone = [tmat_rows(sandbox, expDmin=e, expDmax=e, **common)[0] for e in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)]
        with report.figure("Wet graupel, Fw = 0.5, canting 40 deg, C band") as ax:
            ax.plot(col(table, "D"), col(table, "zdr"), "o-", label="whole table in one run")
            ax.plot(col(alone, "D"), col(alone, "zdr"), "x", ms=10, label="each diameter run alone")
            ax.set(xlabel="D (mm)", ylabel="Zdr (dB)")
        for a, b in zip(alone, table):
            assert a["D"] == b["D"]
            assert a["zdr"] == pytest.approx(b["zdr"], rel=1e-4)

    def test_wet_graupel_without_water_equals_dry_graupel(self, sandbox, report):
        """Wet graupel with liquid fraction Fw = 0 gives the same scattering as dry graupel."""
        common = dict(band="C", ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, expDmin=-1.0, expDmax=1.0, expDstep=0.1)
        wet = tmat_rows(sandbox, type="wg", DIEL="LBwetgr", ARfunc="RYwg", **common)
        dry = tmat_rows(sandbox, type="gg", DIEL="RY19dry", ARfunc="RYdg", **common)
        with report.figure("Graupel at C band, 0 C", ncols=2) as ax:
            for rows, lab, st in ((dry, "dry graupel (gg)", "-"), (wet, "wet graupel (wg), Fw = 0", "--")):
                ax[0].semilogx(col(rows, "D"), col(rows, "zhh"), st, label=lab)
                ax[1].semilogx(col(rows, "D"), col(rows, "zdr"), st, label=lab)
            ax[0].set(xlabel="D (mm)", ylabel="Zh (dBZ)")
            ax[1].set(xlabel="D (mm)", ylabel="Zdr (dB)")
        for rw, rd in zip(wet, dry):
            assert rw["zhh"] == pytest.approx(rd["zhh"], abs=0.01)
            assert rw["zdr"] == pytest.approx(rd["zdr"], abs=0.001)

    def test_mgwma08_dielectric_option(self, sandbox, report):
        """Wet graupel with the MGwMA08 dielectric option: with Fw = 0 it equals dry graupel, and a fully melted
        particle (Fw = 1) is brighter than the dry one (the water volume fraction used to be left undefined).

        Zh is not monotonic in Fw for small particles (the particle shrinks while melting), so only the two ends
        are tested."""
        common = dict(band="C", ARcnst=0.8, Tcmin=0.0, Tcmax=0.0, expDmin=0.0, expDmax=1.0, expDstep=0.1)
        wet = tmat_rows(sandbox, type="wg", DIEL="MGwMA08", ARfunc="RYwg", Fwmin=0.0, Fwmax=1.0, Fwstep=0.25,
                        **common)
        dry = tmat_rows(sandbox, type="gg", DIEL="RY19dry", ARfunc="RYdg", **common)
        by_fw = {}
        for r in wet:
            by_fw.setdefault(r["Fw"], []).append(r)
        with report.figure("Wet graupel with MGwMA08, C band, 0 C") as ax:
            ax.plot(col(dry, "D"), col(dry, "zhh"), "k--", label="dry graupel")
            for fw, rows in by_fw.items():
                ax.plot(col(rows, "D"), col(rows, "zhh"), label=f"Fw = {fw}")
            ax.set(xlabel="D (mm)", ylabel="Zh (dBZ)")
        for rw, rd in zip(by_fw[0.0], dry):
            assert rw["zhh"] == pytest.approx(rd["zhh"], abs=0.01)
        for r0, r1 in zip(by_fw[0.0], by_fw[1.0]):
            assert r1["zhh"] > r0["zhh"], r0["D"]

    @pytest.mark.parametrize("typeh,arfunc,arcnst,dmin", [("ss", "AUds", 0.75, -1.0), ("ii", "CNST", 0.2, -2.0)],
                             ids=["snow", "pristine ice"])
    def test_particle_density_not_above_solid_ice(self, sandbox, report, typeh, arfunc, arcnst, dmin):
        """The particle density derived from the mass-diameter law never exceeds solid ice (916 kg/m3).
        Without a cap, the ICE3 laws give up to 1900 kg/m3 for the smallest crystals."""
        sandbox.run_tmat(default_param(type=typeh, ARfunc=arfunc, ARcnst=arcnst, DIEL="RY19dry", Tcmin=-10.0,
                                       Tcmax=-10.0, expDmin=dmin, expDmax=1.0, expDstep=0.1))
        rows = resu_rows(sandbox, typeh, "S")
        with report.figure(f"Density of {typeh} ({arfunc}, axis ratio {arcnst})") as ax:
            ax.loglog(col(rows, "D"), col(rows, "RHOP"), "o-", label="particle density")
            ax.axhline(RHOI, color="r", ls="--", label="solid ice")
            ax.set(xlabel="equivalent diameter D (mm)", ylabel="density (kg/m3)")
        for r in rows:
            assert r["RHOP"] <= RHOI, r

    def test_tmat_exit_status_when_tmatrix_does_not_converge(self, sandbox):
        """When the T-matrix does not converge (here fully melted 100 mm hail at W band, beyond the method's
        limits), Tmat stops with a non-zero exit status, so that a truncated table cannot go unnoticed."""
        p = default_param(type="wh", band="W", DIEL="LBwetgr", ARfunc="RYwg", ARcnst=0.8, Tcmin=0.0, Tcmax=0.0,
                          Fwmin=1.0, Fwmax=1.0, expDmin=2.0, expDmax=2.0)
        (sandbox.root / "param" / "tmp_config").write_text(format_param(p))
        (sandbox.root / "tables" / "wh").mkdir()
        r = subprocess.run(["./Tmat"], cwd=sandbox.root / "src", capture_output=True, text=True, timeout=600)
        assert r.returncode != 0 and "CONVERGENCE IS NOT OBTAINED" in r.stderr

    def test_large_hail_at_w_band_converges(self, sandbox):
        """Dry hail of 100 mm at W band (the upper limit of TmatParam_Whh_default) is computed: it needs larger
        T-matrix arrays (NPN1 = 160 in ampld.par.f; with 100 the table stopped near 70 mm)."""
        rows = tmat_rows(sandbox, type="hh", band="W", DIEL="RY19dry", ARfunc="RYdg", ARcnst=0.8, Tcmin=0.0,
                         Tcmax=0.0, expDmin=2.0, expDmax=2.0)
        assert len(rows) == 1 and rows[0]["D"] == 100.0 and -50 < rows[0]["zhh"] < 100

    def test_tmat_fails_without_type_subfolder(self, tmp_path, binaries):
        """Tmat writes into tables/<type>/ and stops if that folder is missing, which is why
        generate_tables.sh must create it (see section 5)."""
        sb = Sandbox(tmp_path / "gen", *binaries)
        (sb.root / "param" / "tmp_config").write_text(format_param(default_param()))
        r = subprocess.run(["./Tmat"], cwd=sb.root / "src", capture_output=True, text=True, timeout=60)
        assert not (sb.root / "tables" / "rr").exists()
        assert "Cannot open file" in r.stderr or r.returncode != 0


# =============================================================================================
ICE3_RAIN = dict(a=524.0, b=3.0, C=8e6, X=-1.0, c=842.0, d=0.8)


def ice3_rain_moments(M):
    """Analytic Z (dBZ) and rain rate (mm/h) for the ICE3 exponential rain PSD and content M (kg/m3)."""
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
class TestTmatInt:
    """3. TmatInt: integration over the particle size distribution (TmatCoefInt tables)

    TmatInt sums the per-diameter coefficients over the ICE3/LIMA size distribution for each
    content M. For Rayleigh spheres with an exponential distribution, Z and the rain rate have
    exact formulas, which gives a direct check of the integration.
    """

    def test_coefint_header(self, sphere_rain_table):
        """The TmatCoefInt file has the expected header and one row per content M
        (10^-7 to 10^-1.5 kg/m3, step 0.01 in log10)."""
        header, rows = sphere_rain_table
        assert header["micro"] == "ICE3" and float(header["LAM"]) == 106.2
        assert len(rows) == round((-1.5 + 7.0) / 0.01) + 1
        assert rows[0]["M"] == pytest.approx(1e-7, rel=1e-3)
        assert rows[-1]["M"] == pytest.approx(10 ** -1.5, rel=1e-3)

    def test_integrated_reflectivity_and_rain_rate_match_analytic(self, sphere_rain_table, report):
        """For spherical drops at S band (Rayleigh) and the ICE3 Marshall-Palmer distribution, the table gives
        the exact Z within 0.05 dB (M <= 1 g/m3) and the exact rain rate within 3 % (1 mg/m3 <= M <= 1 g/m3).

        Above 1 g/m3 the 8 mm upper limit of the diameter range cuts the largest drops, so Z is underestimated:
        this is visible on the figure but not tested."""
        _, rows = sphere_rain_table
        m = col(rows, "M")
        za, ra = zip(*(ice3_rain_moments(x) for x in m))
        with report.figure("Rain spheres, S band, ICE3: table vs exact formula", ncols=2) as ax:
            g = [x * 1e3 for x in m]
            ax[0].semilogx(g, col(rows, "zhhR"), label="TmatInt (Rayleigh)")
            ax[0].semilogx(g, za, "--", label="exact")
            ax[0].set(xlabel="rain content (g/m3)", ylabel="Zh (dBZ)")
            ax[1].semilogx(g, [a - b for a, b in zip(col(rows, "zhhR"), za)], label="Zh error (dB)")
            ax[1].semilogx(g, [10 * math.log10(max(a, 1e-12) / b) for a, b in zip(col(rows, "RR"), ra)],
                           label="rain rate error (dB)")
            ax[1].axvspan(1e-4, 1, color="green", alpha=0.08, label="tested range")
            ax[1].set(xlabel="rain content (g/m3)", ylabel="table - exact (dB)", ylim=(-1, 0.3))
        for r, z, rr in zip(rows, za, ra):
            if r["M"] <= 1e-3:
                assert r["zhhR"] == pytest.approx(z, abs=0.05), r["M"]
            if 1e-6 <= r["M"] <= 1e-3:
                assert r["RR"] == pytest.approx(rr, rel=0.03), r["M"]

    def test_integrated_rain_variables_increase_with_content(self, sandbox, report):
        """Oblate rain at C band: Zh, Zdr, Kdp and Ah increase with the rain content (up to 10 g/m3)."""
        diff = sandbox.run_tmat(default_param(band="C", expDmin=-2.0, expDmax=0.9, expDstep=0.02))
        _, rows = read_coefint(sandbox.run_tmatint(diff, "t", "rr", "C"))
        rows = [r for r in rows if r["M"] <= 1e-2]
        g = [r["M"] * 1e3 for r in rows]
        with report.figure("Rain at C band, 10 C, ICE3", nrows=2, ncols=2, height=3.0) as ax:
            for a, (k, lab) in zip(ax, (("zhh", "Zh (dBZ)"), ("zdr", "Zdr (dB)"), ("kdp", "Kdp (deg/km)"),
                                        ("Ah", "Ah (dB/km)"))):
                a.semilogx(g, col(rows, k))
                a.set(xlabel="rain content (g/m3)", ylabel=lab)
        for k in ("zhh", "zdr", "kdp", "Ah"):
            v = col(rows, k)
            assert all(b >= a for a, b in zip(v, v[1:])), k

    def test_two_moment_ice_reflectivity_trends(self, sandbox, report):
        """Two-moment pristine ice: the table covers 91 concentrations (10^0 to 10^9 m-3), and for the same
        mass, fewer (hence bigger) crystals give a higher Zh."""
        diff = sandbox.run_tmat(default_param(type="ii", band="C", DIEL="RY19dry", ARfunc="CNST", ARcnst=1.0,
                                              Tcmin=-20.0, Tcmax=-20.0, expDmin=-2.0, expDmax=0.5, expDstep=0.05))
        _, rows = read_coefint(sandbox.run_tmatint(diff, "t", "ii", "C"))
        by_n = {}
        for r in rows:
            by_n.setdefault(r["P3"], []).append(r)
        with report.figure("Pristine ice at C band, -20 C (Rayleigh column)") as ax:
            for n in sorted(by_n)[20::20]:
                rr = [r for r in by_n[n] if r["zhhR"] > -100]
                ax.semilogx([r["M"] * 1e3 for r in rr], col(rr, "zhhR"), label=f"N = {n:.0e} m-3")
            ax.set(xlabel="ice content (g/m3)", ylabel="Zh (dBZ)")
        assert len(by_n) == 91
        n_lo, n_hi = sorted(by_n)[40], sorted(by_n)[60]
        for lo, hi in zip(by_n[n_lo], by_n[n_hi]):
            assert lo["M"] == hi["M"]
            if 1e-6 <= lo["M"] <= 1e-4:
                assert lo["zhhR"] > hi["zhhR"]

    def test_two_moment_table_with_several_liquid_fractions(self, sandbox, report):
        """A 2-moment species (LIMA rain) integrated from a T-matrix table that also contains other liquid
        fractions gives the same result as from a table with Fw = 0 only (rows used to be read at the wrong place)."""
        # two temperatures: with one temperature and one elevation the wrong row index happened to be right
        common = dict(band="C", Tcmin=0.0, Tcmax=10.0, Tcstep=10.0, expDmin=-2.0, expDmax=0.9, expDstep=0.05)
        one = sandbox.run_tmat(default_param(**common))
        _, r1 = read_coefint(sandbox.run_tmatint(one, "a", "rr", "C", micro="LIMA", moment="2M"))
        two = sandbox.run_tmat(default_param(Fwmin=0.0, Fwmax=0.5, Fwstep=0.5, **common))
        _, r2 = read_coefint(sandbox.run_tmatint(two, "b", "rr", "C", micro="LIMA", moment="2M"))
        n = sorted({r["P3"] for r in r1})[50]  # 1e5 m-3
        with report.figure("LIMA rain (2 moments), N = 1e5 m-3, C band, 10 C") as ax:
            for rows, lab, st in ((r1, "T-matrix table with Fw = 0 only", "-"), (r2, "table with Fw = 0 and 0.5", "--")):
                rows = [r for r in rows if r["P3"] == n and r["Tc"] == 10.0]
                ax.semilogx([r["M"] * 1e3 for r in rows], col(rows, "zhh"), st, label=lab)
            ax.set(xlabel="rain content (g/m3)", ylabel="Zh (dBZ)")
        assert len(r1) == len(r2)
        for a, b in zip(r1, r2):
            assert a["zhh"] == b["zhh"] and a["kdp"] == b["kdp"], (a["M"], a["P3"])

    def test_snow_lower_diameter_limit(self, sandbox, report):
        """ICE3 snow at C band, -10 C: with the diameter range starting at 0.01 mm, Zh increases with the content
        from 0.01 to 2 g/m3, and is unchanged above 0.3 g/m3 compared with a 0.1 mm start.

        The ICE3 snow distribution has a mean diameter of 0.02 mm at 0.01 g/m3 and 1.4 micron at 0.001 g/m3,
        so below about 0.01 g/m3 even a 0.01 mm start misses most of it (Zh too low). Above about 2.4 g/m3 the
        largest snowflakes exceed the 50 mm upper limit and Zh stops increasing. Both limits are visible on the
        figure but not tested."""
        common = dict(type="ss", band="C", DIEL="RY19dry", ARfunc="AUds", ARcnst=0.75, Tcmin=-10.0, Tcmax=-10.0,
                      expDmax=1.7, expDstep=0.02)
        res = {}
        for dmin in (-1.0, -2.0):
            diff = sandbox.run_tmat(default_param(expDmin=dmin, **common))
            res[dmin] = read_coefint(sandbox.run_tmatint(diff, f"d{int(-dmin)}", "ss", "C"))[1]
        with report.figure("ICE3 snow, C band, -10 C: effect of the smallest diameter of the table") as ax:
            for dmin, rows in res.items():
                rows = [r for r in rows if r["zhh"] > -60]
                ax.semilogx([r["M"] * 1e3 for r in rows], col(rows, "zhh"), label=f"Dmin = {10 ** dmin:g} mm")
            ax.axvspan(1e-2, 2, color="green", alpha=0.08, label="tested range")
            ax.set(xlabel="snow content (g/m3)", ylabel="Zh (dBZ)")
        rows = [r for r in res[-2.0] if 1e-5 <= r["M"] <= 2e-3]
        z = col(rows, "zhh")
        assert all(b > a for a, b in zip(z, z[1:]))
        old = {r["M"]: r["zhh"] for r in res[-1.0]}
        assert all(abs(r["zhh"] - old[r["M"]]) < 0.1 for r in res[-2.0] if 3e-4 <= r["M"] <= 3e-3)

    def test_unknown_type_gives_clear_error(self, sandbox):
        """TmatInt stops with a clear message (instead of a Fortran end-of-file crash) when the type is not in
        the constants file."""
        r = subprocess.run([str(sandbox.root / "src" / "TmatInt"), str(sandbox.root / "src"), "t", "zz", "C",
                            "ICE3", "1M"], capture_output=True, text=True, timeout=30)
        assert r.returncode != 0 and "not found" in r.stderr

    def test_wet_graupel_integration_without_water_equals_dry_graupel(self, graupel_tables, report):
        """Integrated wet graupel with Fw = 0 equals dry graupel (within 0.02 dB for M >= 0.1 g/m3),
        and melting graupel (Fw = 0.5) is brighter than dry graupel."""
        rwg, rgg = graupel_tables
        dry_part = [r for r in rwg if r["P3"] == 0.0]
        wet_part = [r for r in rwg if r["P3"] == 0.5]
        with report.figure("Graupel at C band, 0 C, ICE3") as ax:
            for rows, lab, st in ((rgg, "dry graupel (gg)", "-"), (dry_part, "wet graupel, Fw = 0", "--"),
                                  (wet_part, "wet graupel, Fw = 0.5", "-")):
                rows = [r for r in rows if r["zhh"] > -100]
                ax.semilogx([r["M"] * 1e3 for r in rows], col(rows, "zhh"), st, label=lab)
            ax.set(xlabel="graupel content (g/m3)", ylabel="Zh (dBZ)")
        assert len(dry_part) == len(rgg) == len(wet_part)
        for a, b in zip(dry_part, rgg):
            if a["M"] >= 1e-4:
                assert a["zhh"] == pytest.approx(b["zhh"], abs=0.02), a["M"]
        for a, b in zip(wet_part, rgg):
            if 1e-5 <= a["M"] <= 1e-3:
                assert a["zhh"] > b["zhh"]

    def test_wet_graupel_and_graupel_same_normalisation_at_low_content(self, graupel_tables, report):
        """Same comparison as above, but at all contents: wet graupel with Fw = 0 and dry graupel agree
        within 0.02 dB down to 10^-7 kg/m3 (wet species used to be rescaled to the full content M,
        dry ones not: all species now keep the mass inside the diameter range of the table)."""
        rwg, rgg = graupel_tables
        dry_part = [r for r in rwg if r["P3"] == 0.0]
        with report.figure("Zh(wet graupel, Fw = 0) - Zh(dry graupel)") as ax:
            pts = [(a["M"] * 1e3, a["zhh"] - b["zhh"]) for a, b in zip(dry_part, rgg) if a["zhh"] > -100]
            ax.semilogx(*zip(*pts))
            ax.set(xlabel="graupel content (g/m3)", ylabel="difference (dB)")
        for a, b in zip(dry_part, rgg):
            assert a["zhh"] == pytest.approx(b["zhh"], abs=0.02), a["M"]


# =============================================================================================
PSD_DRIVER = """
program t
  implicit none
  real*8 :: lamb, n
  real :: d
  integer :: i
  do i = 0, 40
    d = 1e-4 * 10**(i*0.05)
    ! first call with another content, to leave a different value in lamb (as in TmatInt loops)
    call PSD(1e-6, d, 273.15, 'ICE3', 0.0, 'rr', 524., 3., 1., 1., 8e6, -1., 15., 1, lamb, n)
    call PSD(%(M)s, d, 273.15, 'ICE3', 0.0, 'rr', 524., 3., 1., 1., 8e6, -1., 15., 1, lamb, n)
    write(*,'(3ES24.15)') d, lamb, n
  end do
end program
"""


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
class TestPSD:
    """4. PSD subroutine of TmatrixIntegration.f90

    The PSD subroutine is compiled alone with a small driver program and compared with the
    analytic Marshall-Palmer distribution used by ICE3 for rain.
    """

    def test_psd_ice3_rain_marshall_palmer(self, tmp_path, report):
        """For ICE3 rain at 1 g/m3, PSD returns N(D) = 8e6 exp(-lambda D) (m-4), whatever was computed
        in the previous call (No used to be computed from the lambda of the previous call)."""
        src = (SRC / "TmatrixIntegration.f90").read_text()
        (tmp_path / "psd.f90").write_text(PSD_DRIVER % dict(M="1e-3") + "\n" + src[src.index("SUBROUTINE PSD("):])
        subprocess.run(["gfortran", "-O0", "-w", "psd.f90", "-o", "psd"], cwd=tmp_path, check=True,
                       capture_output=True)
        out = subprocess.run(["./psd"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout
        d, lamb, n = zip(*(map(float, l.split()) for l in out.splitlines()))
        p = ICE3_RAIN
        lam_ref = (1e-3 / (p["a"] * p["C"] * 6)) ** (1 / (p["X"] - p["b"]))
        n_ref = [p["C"] * math.exp(-lam_ref * x) for x in d]
        with report.figure("ICE3 rain size distribution at 1 g/m3") as ax:
            ax.semilogy([x * 1e3 for x in d], n, label="PSD subroutine")
            ax.semilogy([x * 1e3 for x in d], n_ref, "--", label="Marshall-Palmer (exact)")
            ax.set(xlabel="D (mm)", ylabel="N(D) (m-4)")
        assert lamb[0] == pytest.approx(lam_ref, rel=1e-5)
        for a, b in zip(n, n_ref):
            assert a == pytest.approx(b, rel=1e-4)


# =============================================================================================
FAKE_TMAT = """#!/bin/bash
# mimics Tmat <param> <out coef> <out resu> <first Tc index> <last Tc index>:
# 3 header lines, then one line per temperature index
echo "$@" >> ../calls_tmat
t=$(sed -n 's/^type://p' "$1"); [ -n "$t" ] || { echo "cannot read $1" >&2; exit 2; }
{ echo h1; echo h2; echo h3; for ((i=$4; i<=$5; i++)); do echo "row $i"; done; } > "$2"
{ echo r1; for ((i=$4; i<=$5; i++)); do echo "resu $i"; done; } > "$3"
"""
FAKE_TMATINT = """#!/bin/bash
# mimics TmatInt <exec_dir> <repo> <type> <band> <micro> <moment>
echo "$@" >> $1/../calls_tmatint
[ -f $1/../tables/$2/TmatCoefDiff_$4$3 ] || { echo "missing CoefDiff" >&2; exit 2; }
echo fake > $1/../tables/$2/TmatCoefInt_$5_$6_$4$3
"""
HYDRO_LIST = ["cs", "cl", "rr", "ii", "ss", "gg", "hh", "wg"]


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


class TestGenerateTablesScript:
    """5. generate_tables.sh

    The script is run in a temporary copy of the repository where Tmat and TmatInt are replaced by
    small fake programs that only create the expected files (and fail like the real ones when an
    input is missing). This checks the options and the file handling of the script in a second,
    without computing any table.
    """

    def test_parallel_blocks_are_concatenated_in_order(self, script_root, report):
        """With -j 4, Tmat is run as 4 processes on consecutive blocks of temperatures and their outputs are joined
        in temperature order, with the header written once."""
        run_script(script_root, "--default", "--band", "C", "-j", "4")
        gen = script_root / "tables_generator"
        calls = (gen / "calls_tmat").read_text().splitlines()
        lines = (gen / "tables" / "default" / "TmatCoefDiff_Crr").read_text().splitlines()
        ntc = 61  # TmatParam_Crr_default: -20..40 C
        report.table("Tmat processes started for rain (arguments: param, outputs, first and last Tc index)",
                     ["call"], [[c.split("/")[-1] if False else " ".join(a.split("/")[-1] for a in c.split())]
                                for c in calls if "Crr" in c])
        assert sum("Crr" in c for c in calls) == 4
        assert lines[:3] == ["h1", "h2", "h3"]
        assert lines[3:] == [f"row {i}" for i in range(ntc)]

    def test_script_rejects_unknown_band(self, script_root):
        """An unknown band (--band Q) is rejected with an error message."""
        r = run_script(script_root, "--default", "--band", "Q")
        assert r.returncode != 0 and "invalid band" in r.stdout

    def test_default_mode_generates_all_tables(self, script_root, report):
        """On a fresh clone, 'generate_tables.sh --default --band C' produces a TmatCoefInt table for all
        8 hydrometeor types in tables/default/."""
        r = run_script(script_root, "--default", "--band", "C")
        out = script_root / "tables_generator" / "tables" / "default"
        report.table("Tables produced", ["type", "TmatCoefDiff", "TmatCoefInt"],
                     [[h, "yes" if (out / f"TmatCoefDiff_C{h}").exists() else "no",
                       "yes" if (out / f"TmatCoefInt_ICE3_1M_C{h}").exists() else "no"] for h in HYDRO_LIST])
        missing = [h for h in HYDRO_LIST if not (out / f"TmatCoefInt_ICE3_1M_C{h}").exists()]
        assert not missing, f"no table produced for {missing}"

    def test_default_mode_works_when_type_folders_exist(self, script_root):
        """If the tables/<type>/ folders already exist, --default produces both tables for all 8 types."""
        for h in HYDRO_LIST:
            (script_root / "tables_generator" / "tables" / h).mkdir()
        run_script(script_root, "--default", "--band", "C")
        out = script_root / "tables_generator" / "tables" / "default"
        for h in HYDRO_LIST:
            assert (out / f"TmatCoefDiff_C{h}").exists() and (out / f"TmatCoefInt_ICE3_1M_C{h}").exists(), h

    def test_default_mode_reuses_existing_tables(self, script_root):
        """Tables already present in tables/default/ are reused: Tmat and TmatInt are not run again."""
        out = script_root / "tables_generator" / "tables" / "default"
        out.mkdir()
        for h in HYDRO_LIST:
            (out / f"TmatCoefDiff_C{h}").write_text("x")
            (out / f"TmatCoefInt_ICE3_1M_C{h}").write_text("x")
        run_script(script_root, "--default", "--band", "C")
        assert not (script_root / "tables_generator" / "calls_tmat").exists()
        assert not (script_root / "tables_generator" / "calls_tmatint").exists()

    def test_newconf_mode_runs_generator(self, script_root, report):
        """'--newConf myconf --band C' copies each base file TmatParam_C<type> to TmatParam_C<type>_myconf and
        produces both tables for the 8 types in tables/myconf/."""
        r = run_script(script_root, "--newConf", "myconf", "--band", "C")
        gen = script_root / "tables_generator"
        out = gen / "tables" / "myconf"
        report.table("Files produced", ["type", "param copy", "TmatCoefDiff", "TmatCoefInt"],
                     [[h] + ["yes" if f.exists() else "no" for f in (gen / "param" / f"TmatParam_C{h}_myconf",
                      out / f"TmatCoefDiff_C{h}", out / f"TmatCoefInt_ICE3_1M_C{h}")] for h in HYDRO_LIST])
        assert r.returncode == 0, r.stdout
        for h in HYDRO_LIST:
            assert (gen / "param" / f"TmatParam_C{h}_myconf").exists(), h
            assert (out / f"TmatCoefInt_ICE3_1M_C{h}").exists(), h

    def test_script_can_be_launched_from_another_directory(self, script_root, tmp_path):
        """The script finds its folders from its own location, so it also works when launched from elsewhere."""
        r = subprocess.run(["bash", str(script_root / "generate_tables.sh"), "--default", "--band", "C"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, r.stdout
        assert (script_root / "tables_generator" / "tables" / "default" / "TmatCoefInt_ICE3_1M_Crr").exists()

    def test_script_exit_status_reports_failures(self, script_root):
        """If one step fails (here a missing parameter file), the other types are still processed but the
        script ends with an error status and a message, so a failed generation does not go unnoticed."""
        (script_root / "tables_generator" / "param" / "TmatParam_Cgg_default").unlink()
        r = run_script(script_root, "--default", "--band", "C")
        assert r.returncode != 0
        assert "step(s) failed" in r.stdout
        assert (script_root / "tables_generator" / "tables" / "default" / "TmatCoefInt_ICE3_1M_Crr").exists()

    @pytest.mark.parametrize("args", [["--hydro", "rr", "--band", "C"], ["--default", "--band", "C", "--arf", "CNST"]],
                             ids=["--hydro", "--arf"])
    def test_removed_edit_mode_options_are_rejected(self, script_root, args):
        """The unfinished edit mode (--hydro, --arf, ...) was removed: its options are refused with a message."""
        r = run_script(script_root, *args)
        assert r.returncode != 0 and "Unknown option" in r.stdout


@needs_fortran
class TestGenerateTablesEndToEnd:
    """6. generate_tables.sh with the real Tmat and TmatInt

    The script is run on a copy of the repository where the default C-band parameter files of three types
    are replaced by small versions, so the whole chain (parameter file, Tmat, TmatInt) runs in seconds.
    """

    TYPES = {"rr": dict(), "gg": dict(type="gg", DIEL="RY19dry", ARfunc="RYdg", ARcnst=0.8),
             "wg": dict(type="wg", DIEL="LBwetgr", ARfunc="RYwg", ARcnst=0.8, canting="analyt", sigbeta=40.0,
                        Fwmax=1.0, Fwstep=0.5)}

    @pytest.fixture
    def small_repo(self, tmp_path, binaries):
        root = tmp_path / "repo"
        gen = root / "tables_generator"
        (gen / "src").mkdir(parents=True)
        (gen / "tables").mkdir()
        (gen / "param").mkdir()
        for f in PARAM.glob("*_constants.txt"):
            shutil.copy(f, gen / "param")
        for h, over in self.TYPES.items():
            p = default_param(band="C", Tcmin=-4.0, Tcmax=4.0, Tcstep=2.0, ELEVmax=90.0, ELEVstep=90.0,
                              expDmin=-1.0, expDmax=0.9, expDstep=0.1, **{**dict(type=h), **over})
            (gen / "param" / f"TmatParam_C{h}_default").write_text(format_param(p))
        shutil.copy(REPO / "generate_tables.sh", root)
        shutil.copy(binaries[0], gen / "src" / "Tmat")
        shutil.copy(binaries[1], gen / "src" / "TmatInt")
        return root

    def test_parallel_run_gives_identical_tables(self, small_repo, report):
        """The tables produced with -j 3 (Tmat split over 3 processes) are byte-identical to those of a
        sequential run (-j 1), for rain, graupel and wet graupel (2 elevations, 5 temperatures)."""
        tables = small_repo / "tables_generator" / "tables" / "default"
        results = {}
        for jobs in ("1", "3"):
            shutil.rmtree(tables, ignore_errors=True)
            r = run_script(small_repo, "--default", "--band", "C", "-j", jobs)
            results[jobs] = {f.name: f.read_bytes() for f in tables.glob("TmatCoef*")}
            assert r.returncode == 1  # the other 5 types have no parameter file in this small copy
            assert "5 step(s) failed" in r.stdout
        report.table("Files compared", ["file", "lines", "identical"],
                     [[n, results["1"][n].count(b"\n"), results["1"][n] == results["3"].get(n)]
                      for n in sorted(results["1"])])
        assert len(results["1"]) == 6  # TmatCoefDiff + TmatCoefInt for 3 types
        assert results["1"] == results["3"]
