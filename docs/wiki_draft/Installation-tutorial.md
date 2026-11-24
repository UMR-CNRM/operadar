# OPERADAR — Installation Guide

OPERADAR requires **Python 3.12 or 3.13** and must be installed in an isolated virtual environment.
This guide provides **copy-paste safe**, step-by-step instructions.

> Python 3.14 is not supported yet: the library that reads the AROME `.fa` files (`falfilfa4py`, used by
> `epygram`) has no package for it, so `pip install` refuses Python 3.14 with
> *"requires a different Python"*. Python 3.10 and 3.11 are not supported either.

## Prerequisites

- **Python 3.12 or 3.13** (verify with `python3 --version`). If your system Python is another version, see
  [Getting Python 3.12](#getting-python-312) below.
- **venv** (included with Python)
- **Git**
- No compiler is needed for OPERADAR itself: all the dependencies (including `epygram` and the FA library)
  are installed as ready-made packages by `pip`. A Fortran compiler is only needed to generate the lookup
  tables (see [Generation of the lookup tables](https://github.com/UMR-CNRM/operadar/wiki/Generation-of-the-lookup-tables)).
- For HPC at Météo France (belenos), if `git clone` fails with a certificate error, disable the certificate
  check for GitHub only (and not for every server), in the `.gitconfig` file:
```bash
[http "https://github.com/"]
        sslVerify = false
```

---
0. **For belenos only**
   ```bash
   module load python/3.12.12
   ```

1. **Create a directory for virtual environments** (safe to re-run):
   ```bash
   mkdir -p $HOME/virtual_envs
   cd $HOME/virtual_envs
   ```

2. **Create the virtual environment** (use `python3.12` or `python3.13` if `python3` is another version):
   ```bash
   python3 -m venv operadar_env
   ```

3. **Activate the environment**:
   ```bash
   source operadar_env/bin/activate
   ```

4. **Verify the activation** (it should point to `.../operadar_env/bin/python`, version 3.12 or 3.13):
   ```bash
   which python
   python --version
   ```

5. **Upgrade pip**:
   ```bash
   pip install -U pip setuptools wheel
   ```

6. **Clone the repository**:
   ```bash
   cd $HOME
   git clone https://github.com/UMR-CNRM/operadar.git
   cd operadar
   ```

7. **Install OPERADAR in editable mode** (about 2 minutes):
   ```bash
   pip install --prefer-binary -e .
   ```
   The editable mode (`-e`) is **required**: operadar reads the configuration files in `configFiles/` of the
   cloned repository, and changes to the source code are immediately active.
   To also run the tests, use `pip install --prefer-binary -e ".[dev]"` (adds `pytest`).

8. **Verify the installation**:
   ```bash
   python -c "import operadar, epygram; print('OPERADAR OK, epygram', epygram.__version__)"
   ```
   Optionally, run the tests from the root of the repository (after `pip install -e ".[dev]"`):
   ```bash
   pytest
   ```

9. **Deactivate when finished** (optional):
    ```bash
    deactivate
    ```

## Quick Verification

To verify your installation in a new terminal later:

1. **Activate your environment**:
   ```bash
   source $HOME/virtual_envs/operadar_env/bin/activate
   ```

2. **Test the import**:
   ```bash
   python -c "import operadar; print(operadar.__file__)"
   ```

---

## <a name="getting-python-312"></a>Getting Python 3.12

If `python3 --version` gives another version (e.g. 3.14 on recent Linux distributions) and no
`python3.12` / `python3.13` command is available:

- **On an HPC**: load a module, e.g. `module load python/3.12.12` (belenos), and check `module avail python`.
- **With uv** (no administrator rights needed): install [uv](https://docs.astral.sh/uv/), then replace steps 2
  to 5 by
  ```bash
  uv venv --python 3.12 operadar_env     # downloads Python 3.12 if needed
  source operadar_env/bin/activate
  uv pip install pip                     # so that the "pip" commands of the next steps work
  ```
- **With conda**: `conda create -n operadar_env python=3.12` then `conda activate operadar_env`, and continue
  at step 5. Do not mix conda and venv: install everything with `pip` inside this conda environment.

## Notes

- **Editable install** (`-e .`): installs OPERADAR in development mode, so changes to the source code take
  effect immediately without reinstalling. It is also needed to find the `configFiles/` directory.
- **Run operadar from the root of the repository**: the paths of the configuration files
  (`./modelFiles/`, `./operadarFiles/`, `./tables_generator/tables/`) are relative to it.
- **Lookup tables**: they are not included in the repository (about 17 GB for the David2026PhD
  configuration). Generate them, or ask for them (see the [Home](https://github.com/UMR-CNRM/operadar/wiki) page).
- **Dependency conflicts**: install OPERADAR in its own environment. If you add other packages to it, `pip check`
  lists the conflicts, if any.


## Running OPERADAR

See the [Execution Guide](https://github.com/UMR-CNRM/operadar/wiki/Execution-of-operadar) for detailed usage instructions.
