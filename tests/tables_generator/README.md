# Tests of the lookup-table generator

These tests check `tables_generator/` (the Fortran programs `Tmat` and `TmatInt`), its
parameter files and `generate_tables.sh`. They do **not** touch your `tables/` folder:
everything runs on very small tables in a temporary directory, in about 10 seconds.

## Run

From the repository root, in the operadar environment (needs `pytest`; `matplotlib` for the figures):

```bash
LAPACK_DIR=/path/to/lapack pytest tests/tables_generator
```

`LAPACK_DIR` is the folder containing `liblapack.a` and `librefblas.a` (the same as in
`tables_generator/src/makefile_*`). If it is not given, `tables_generator/src/libFortran/lapack`
is used. Without gfortran or LAPACK the Fortran tests are skipped and only the
parameter-file and script tests run.

## Read the report

At the end of the run pytest prints

```
------------------------- table-generator test report --------------------------
open /.../tests/tables_generator/report/test_report.html
```

Open that file in a browser. It is self-contained, so you can send it to someone.
It has one section per part of the generator:

1. **Parameter files**: format and consistency of `param/TmatParam_*` and `*_constants.txt`.
2. **Tmat**: scattering by one particle, compared with known limits (Rayleigh, spheres,
   vertical incidence, no canting, D^6 law, density below solid ice...).
3. **TmatInt**: integration over the size distribution, compared with exact formulas
   (Z and rain rate for Marshall-Palmer rain), trends with content and concentration.
4. **PSD subroutine**: size distribution compared with Marshall-Palmer.
5. **generate_tables.sh**: options and file handling, using fake Tmat/TmatInt.

Each test is a card that gives what is checked, the result and, for the physics tests,
a figure of what was compared. The shaded areas show the range where the numbers are checked.

| Result | Meaning |
|---|---|
| PASSED | works as expected |
| FAILED / ERROR | something is broken: read the error in the card |
| KNOWN BUG | a documented bug that is not fixed yet (the reason is shown) |
| BUG FIXED? | a known-bug test now passes: remove its `xfail` marker in the test file |
| SKIPPED | not run (missing gfortran or LAPACK) |

Options:
* `--tg-report=other/path.html` writes the report elsewhere; `--tg-report=` disables it.
* `-v` lists every test in the terminal as well.

## Add a test

Tests are in `test_tables_generator.py`, grouped in classes (one class = one report section).
The first line of the docstring is shown as "What is checked". To add a figure, request the
`report` fixture:

```python
def test_something(self, sandbox, report):
    """One sentence saying what is checked."""
    rows = tmat_rows(sandbox, band="C", expDmin=0.0, expDmax=0.8)
    with report.figure("caption") as ax:
        ax.plot([r["D"] for r in rows], [r["zdr"] for r in rows])
    assert ...
```

Make the figure before the `assert`, so it also appears when the test fails.
`harness.py` holds the helpers that build the programs, write a parameter file
(`default_param(...)`), run `Tmat`/`TmatInt` and read their output files.
