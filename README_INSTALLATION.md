# OPERADAR — Quick installation (operadar_env)

Minimal, copy/paste instructions to create an environment named `operadar_env` and install OPERADAR.
Choose one of the two options below.

Option A — venv + pip (no conda)
```bash
# 1) Create and activate a virtual environment named operadar_env
python3.11 -m venv operadar_env
source operadar_env/bin/activate   # on Windows PowerShell: .\operadar_env\Scripts\Activate.ps1

# 2) Upgrade packaging tools and install Python dependencies from requirements.txt
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt

# 3) Install the package in editable mode (development)
pip install -e .

# 4) Quick smoke test
python -c "import operadar; print('operadar import OK')"

# When done, deactivate
deactivate
```

Option B — conda / mamba (recommended when available)
```bash
# Option 1: create a minimal conda env and use pip inside it
conda create -n operadar_env python=3.11 -y
conda activate operadar_env

# Option 2: if you have an environment.yml, create the env from it (may install additional system-level packages)
# mamba env create -f environment.yml -n operadar_env
# or
# conda env create -f environment.yml -n operadar_env
# conda activate operadar_env

# Install Python dependencies and the package (use pip inside the conda env)
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install -e .

# Quick smoke test
python -c "import operadar; print('operadar import OK')"

# Deactivate when done
conda deactivate
```

Notes
- These instructions assume Python 3.11 is available on your system.
- Environment name used in examples: `operadar_env`.
- If a package install fails with `pip`, try the conda option and install heavy dependencies (e.g., netCDF / geospatial packages) via conda first.
- After installation, run OPERADAR with a configuration from `configFiles/template.py`:
  ```bash
  python -m operadar <input_filename_without_path> <config_name> [--append] [--verbose]
  ```
- For more detailed troubleshooting, see `TROUBLESHOOTING.md`.