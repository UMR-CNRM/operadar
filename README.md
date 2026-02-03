# OPERADAR — Observation oPErator for polarimetric RADAR

OPERADAR computes synthetic dual-polarization radar variables from model fields on a 3D grid.
It is developed at CNRM (National Centre for Meteorological Research, France).

Supported inputs
- AROME `.fa` files (via epygram)
- MesoNH `.nc` files (netCDF, MesoNH format)

Computed outputs (netCDF)
- Polarimetric radar variables: Zh, Zdr, Kdp, rhohv
- Specific attenuations Ah, Av
- Meteorological fields (temperature, contents, concentrations)
- Grid coordinates and metadata

Important links
- Wiki: Installation tutorial and lookup table generation (see repository Wiki)
- Lookup table generator: tables_generator executable included in the project (see wiki page "Generation of the lookup tables")

Quick install summary
- Full conda-based installation (recommended): see `README_INSTALLATION.md` for full instructions.
- Pip / virtualenv installation is possible but requires installing system libs and ecCodes beforehand.

Quick start (after installation)
1. Generate lookup tables (if you don't already have them). See the Wiki: "Generation of the lookup tables".
2. Prepare a configuration file — copy `configFiles/template.py` and adapt parameters.
3. Run the forward operator on one file:
```bash
# Example: run operadar on a single filename and configuration
python -m operadar <input_filename_without_path> <config_name> [--append] [--verbose]

# Or using helper script:
./exec_operadar.sh -f <input_filename> -c <config_name> [--verbose] [--append]
```

Installation (short)
- Recommended: use the provided `environment.yml` (conda-forge) that installs binary dependencies (ecCodes, proj, geos, etc.) and then pip-installs Python packages (including `epygram`).
- Alternate: use `requirements.txt` with a Python virtualenv, but ensure ecCodes and system geospatial libs are installed first.

Files added by the env/installation PR
- `requirements.txt` — pip dependency list (includes `epygram`)
- `environment.yml` — conda environment (Python 3.11) with ecCodes and system libs
- `Dockerfile` — example container using mambaforge
- `.github/workflows/ci.yml` — CI smoke test that builds env and runs a basic import
- `README_INSTALLATION.md` — detailed install instructions (conda + pip)
- `TROUBLESHOOTING.md` — FAQ/troubleshooting (epygram / ecCodes / cartopy)

Testing
- After installing environment, run:
```bash
python -c "import operadar; import epygram; print('operadar OK, epygram', getattr(epygram,'__version__','unknown'))"
```

Contributing
- Fork, create a topic branch, make changes, push and open a PR.
- If you add dependencies or change install instructions, please update `requirements.txt` / `environment.yml` and CI.

License
- This software is governed by the CeCILL-C license (see LICENSE.txt).
