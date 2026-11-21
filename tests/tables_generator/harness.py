"""Helpers to build and run the Fortran table generator in an isolated sandbox.

The two Fortran programs use hard-coded relative paths:
  * Tmat    (run from <root>/src) reads  ../param/tmp_config and ../param/ICE3_constants.txt
            and writes ../tables/<type>/TmatCoefDiff_<band><type>
  * TmatInt (argv[1] = <root>/src) reads <root>/param/<MICRO>_constants.txt and
            <root>/tables/<repo>/TmatCoefDiff_<band><type>, and writes
            <root>/tables/<repo>/TmatCoefInt_<MICRO>_<MOMENT>_<band><type>
so every test builds a throw-away <root> tree that mimics tables_generator/.

Only the standard library is used, so the tests run without the operadar
Python environment.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
GEN = REPO / "tables_generator"
SRC = GEN / "src"
PARAM = GEN / "param"

# Order and fixed-column layout expected by the READ statements of Tmatrix.f
# (key, Fortran edit descriptor of the value). The prefix "key:" is skipped with nX.
PARAM_LAYOUT = [
    ("type", "A2"), ("band", "A2"), ("canting", "A6"), ("sigbeta", "F4.1"),
    ("DIEL", "A7"), ("ARfunc", "A4"), ("ARcnst", "F4.2"), ("DSTYfunc", "A4"),
    ("Frim", "F5.1"), ("LAM", "F5.1"), ("ELEVmin", "F4.1"), ("ELEVmax", "F4.1"),
    ("ELEVstep", "F4.1"), ("Tcmin", "F5.1"), ("Tcmax", "F5.1"), ("Tcstep", "F5.1"),
    ("expDmin", "F5.4"), ("expDmax", "F5.4"), ("expDstep", "F5.2"),
    ("Fwmin", "F4.2"), ("Fwstep", "F4.2"), ("Fwmax", "F4.2"),
]
MICRO_KEYS = ["ICE3", "LIMA", "ICJW", "LIMC"]

# Wavelengths (mm) used in the *_default parameter files
LAM_MM = {"L": 190.4, "S": 106.2, "C": 53.2, "X": 31.9, "Ku": 20.0, "K": 15.0, "Ka": 8.4, "W": 3.15}


def lapack_dir() -> Path | None:
    """Directory holding liblapack.a and librefblas.a (env LAPACK_DIR or src/libFortran/lapack)."""
    candidates = [os.environ.get("LAPACK_DIR"), str(SRC / "libFortran" / "lapack")]
    for c in candidates:
        if c and (Path(c) / "liblapack.a").exists() and (Path(c) / "librefblas.a").exists():
            return Path(c)
    return None


def build(outdir: Path) -> tuple[Path, Path]:
    """Compile Tmat and TmatInt into outdir (needs gfortran + LAPACK)."""
    outdir.mkdir(parents=True, exist_ok=True)
    for f in ("Tmatrix.f", "ampld.par.f", "TmatrixIntegration.f90"):
        shutil.copy(SRC / f, outdir / f)
    lap = lapack_dir()
    subprocess.run(["gfortran", "-O2", "-w", "Tmatrix.f", f"-L{lap}", "-llapack", "-lrefblas", "-o", "Tmat"],
                   cwd=outdir, check=True, capture_output=True)
    subprocess.run(["gfortran", "-O2", "-w", "TmatrixIntegration.f90", "-o", "TmatInt"],
                   cwd=outdir, check=True, capture_output=True)
    return outdir / "Tmat", outdir / "TmatInt"


def format_param(p: dict) -> str:
    """Write a TmatParam file with the same text layout as the files in param/."""
    def fmt(key, val):
        if isinstance(val, str):
            return val
        if key in ("ARcnst", "Fwmin", "Fwstep", "Fwmax"):
            return f"{val:.2f}"
        if key in ("expDmin", "expDmax", "expDstep") or key == "LAM":
            return f"{val:g}" if "." in f"{val:g}" else f"{val:.1f}"
        return f"{val:.1f}"
    lines = [f"{k}:{fmt(k, p[k])}" for k, _ in PARAM_LAYOUT]
    lines += [f"{k}:{p.get(k, '1M')}" for k in MICRO_KEYS]
    return "\n".join(lines) + "\n"


def default_param(**over) -> dict:
    p = dict(type="rr", band="S", canting="nocant", sigbeta=0.0, DIEL="Liebe91", ARfunc="BR02", ARcnst=1.0,
             DSTYfunc="RHOX", Frim=1.0, LAM=LAM_MM["S"], ELEVmin=0.0, ELEVmax=0.0, ELEVstep=1.0,
             Tcmin=10.0, Tcmax=10.0, Tcstep=1.0, expDmin=-1.0, expDmax=0.9, expDstep=0.1,
             Fwmin=0.0, Fwstep=0.1, Fwmax=0.0)
    if "band" in over and "LAM" not in over:
        over["LAM"] = LAM_MM[over["band"]]
    p.update(over)
    return p


class Sandbox:
    """A throw-away copy of the tables_generator tree, with compiled binaries."""

    def __init__(self, root: Path, tmat: Path, tmatint: Path):
        self.root = root
        (root / "src").mkdir(parents=True)
        (root / "tables").mkdir()
        shutil.copytree(PARAM, root / "param", ignore=shutil.ignore_patterns("TmatParam_*", "tmp_config"))
        shutil.copy(tmat, root / "src" / "Tmat")
        shutil.copy(tmatint, root / "src" / "TmatInt")

    def run_tmat(self, param: dict) -> Path:
        typeh, band = param["type"], param["band"]
        (self.root / "param" / "tmp_config").write_text(format_param(param))
        (self.root / "tables" / typeh).mkdir(exist_ok=True)
        r = subprocess.run(["./Tmat"], cwd=self.root / "src", capture_output=True, text=True, timeout=int(os.environ.get("TMAT_TIMEOUT", "600")))
        out = self.root / "tables" / typeh / f"TmatCoefDiff_{band}{typeh}"
        if r.returncode != 0 or not out.exists():
            raise RuntimeError(f"Tmat failed ({r.returncode}):\n{r.stderr[-2000:]}")
        return out

    def run_tmatint(self, coefdiff: Path, repo: str, typeh: str, band: str, micro="ICE3", moment="1M",
                    extra_coefdiff: dict[str, Path] | None = None) -> Path:
        """Integrate a TmatCoefDiff table over the PSD (rain table is needed for melting species)."""
        d = self.root / "tables" / repo
        d.mkdir(exist_ok=True)
        shutil.copy(coefdiff, d / coefdiff.name)
        for name, path in (extra_coefdiff or {}).items():
            shutil.copy(path, d / name)
        r = subprocess.run([str(self.root / "src" / "TmatInt"), str(self.root / "src"), repo, typeh, band, micro,
                            moment], capture_output=True, text=True, timeout=1800)
        out = d / f"TmatCoefInt_{micro}_{moment}_{band}{typeh}"
        if r.returncode != 0 or not out.exists():
            raise RuntimeError(f"TmatInt failed ({r.returncode}):\n{r.stderr[-2000:]}")
        return out


def read_coefdiff(path: Path) -> list[dict]:
    """Rows of a TmatCoefDiff file (3 header lines, whitespace separated)."""
    lines = path.read_text().splitlines()
    names = lines[2].split()
    return [dict(zip(names, map(float, l.split()))) for l in lines[3:] if l.strip()]


def read_coefdiff_header(path: Path) -> dict:
    lines = path.read_text().splitlines()
    return dict(zip(lines[0].split(), lines[1].split()))


def read_coefint(path: Path) -> tuple[dict, list[dict]]:
    """(header parameters, rows) of a TmatCoefInt file (';' separated)."""
    lines = path.read_text().splitlines()
    header = dict(zip(lines[0].split(";"), [v.strip() for v in lines[1].split(";")]))
    names = lines[2].split(";")
    rows = [dict(zip(names, map(float, l.split(";")))) for l in lines[3:] if l.strip()]
    return header, rows
