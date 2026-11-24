# Tmatrix lookup tables
To run, operadar needs pre-computed tmatrix lookup tables.
You can either create these tables yourself following the [lookup tables wiki](https://github.com/UMR-CNRM/operadar/wiki/Generation-of-the-lookup-tables) or request them to clotilde.augros@meteo.fr. The location of the tmatrix tables (directory) has to be written in the configuration file (see below).  

# Configuration file
The user must set up, beforehand, a configuration file based on the template provided `./configFiles/template.py`. Please, do not modify the template directly. Instead, make a copy of the file in `./configFiles/` and name it differently. The configuration file options are explained [here](https://github.com/UMR-CNRM/operadar/wiki/Simulation-configuration).

Ready-to-use examples are provided in `./configFiles/`:

| File | Model | Radar |
|---|---|---|
| `test_AROME_Cband.py` | AROME `.fa`, subdomain over the south-east of France | C band, horizontal beam at each grid point (`cnst_angle = 0`) |
| `test_AROME_Wband_vertical.py` | AROME `.fa`, same subdomain | W band, radar at the ground looking up (`cnst_angle = 90`), with attenuation (`Zh_att`) |
| `test_MesoNH_Cband.py` | MesoNH `.nc`, subdomain around a cyclone near La Réunion | C band, horizontal beam at each grid point |

The paths in the configuration files (`./modelFiles/...`, `./operadarFiles/...`, `./tables_generator/tables/...`) are relative to the directory where operadar is launched: **run operadar from the root of the repository** (`cd $HOME/operadar`).

# Activate your environment
See the [Installation tutorial](https://github.com/UMR-CNRM/operadar/wiki/Installation-tutorial) section. If not already activated, run :
```
$ source $HOME/virtual_envs/operadar_env/bin/activate
```
If you forgot the name of your environment, you can list the available environments with :
```
find $HOME -name "*activate" -type f
```

# Quick execution
For example, from the root of the repository:
```
$ cd $HOME/operadar
$ python -m operadar "historic.arome.franmg-01km30+0018:00.fa" test_AROME_Cband
```

OPERADAR comes with a command line interface. This method is suitable for one file at a time and runs in a terminal. To show the help :
```
>>> $ python -m operadar --help

usage: python3 -m operadar [-h] [--append] [--verbose] filename config

positional arguments:
  filename    Only the name of the input file, not the entire path. Use the
              config file to provide the path to access the file.
  config      Name of the configuration file. Before running the code, you
              need to create a configuration file in ./configFiles/ based on
              the template provided.

options:
  -h, --help  show this help message and exit
  --append    Append the computed dual-polarimetric fields into the input file
              (default value: False)
  --verbose   If argument provided, show more details if activated (default
              value: False)
```

A bash script can also be used for the execution:
```
>>> $ ./exec_operadar.sh 
----------------------------------------------------------------
 This radar forward operator is developped at the CNRM, France. 
----------------------------------------------------------------

Usage: ./exec_operadar.sh -f FILENAME -c CONFIG [--verbose]

  -f FILENAME : Only the filename. Please use the config file to provide the path to access the file.
  -c CONFIG   : Before running the code, you need to create a configuration file in ./configFiles/ based on the template provided.
  --verbose   : Optional, show more details if activated (default value: False)
  --append    : Optional, append the computed dual-polarimetric fields into the input file (default value: False)
                /!\ If this option is activated, no netcdf output file will be saved. Only working for AROME .fa files
```


# Call `operadar()`
`operadar()` can be called inside any python program (run from the root of the repository). This is best adapted
to working with multiple files, or with varying parameters across files. A complete example is
`operadar_multi.py`.

```python
import datetime as dt
import pandas as pd
from operadar.forward_operator import operadar
from operadar.load_config_file import load_configuration_file

config = load_configuration_file('test_AROME_Cband')     # a file of ./configFiles/, with or without .py
read_tables, tables = True, {}
time, end = pd.Timestamp('2025-08-31 18:00'), pd.Timestamp('2025-08-31 19:00')
while time <= end:
    fname = f"historic.arome.franmg-01km30+00{time.strftime('%H:%M')}.fa"
    # the tables are read at the first call only, then passed from one call to the next
    read_tables, tables = operadar(filename=fname,
                                   configuration=config,
                                   read_tables=read_tables,
                                   tables_content=tables,
                                   )
    time += dt.timedelta(minutes=5)
```

**Every option of the configuration file can be overridden** by the argument of the same name when calling
`operadar()` (e.g. `radar_band='X'`, `out_dir_path='./operadarFiles/X/'`); the arguments left to `None` keep the
value of the configuration file.

# `operadar()` arguments

| Argument | Status | Description |
| -------- | ------ | ----------- |
| `filename` | <span style='color:red'>mandatory</span> | Only the name of the file (the folder is `in_dir_path`). |
| `configuration` | <span style='color:red'>mandatory</span> | Configuration loaded with `load_configuration_file('name')`. |
| `read_tables` | <span style='color:green'>optional</span> | Read the lookup tables (`True` at the first call). `operadar()` returns `(read_tables, tables_content)` to pass to the next call. Defaults to `True`. |
| `tables_content` | <span style='color:green'>optional</span> | Content of the lookup tables returned by the previous call (see [To go further](#2gofurther)). |
| `in_dir_path` | <span style='color:green'>optional</span> | Folder of the input files (config: `input_filePath`). |
| `out_dir_path` | <span style='color:green'>optional</span> | Folder of the output files (config: `output_filePath`). |
| `tables_path` | <span style='color:green'>optional</span> | Folder of the lookup tables (config: `path_tables`). |
| `modelname` | <span style='color:green'>optional</span> | `'Arome'` or `'MesoNH'` (config: `model`). |
| `real_case` | <span style='color:green'>optional</span> | `True` or `False` (MesoNH idealized cases). |
| `microphysics_scheme` | <span style='color:green'>optional</span> | `'ICE3'`, `'ICJW'`, `'LIMA'` or `'LIMC'` + an optional name extension (e.g. `'LIMA_noHail'`); only the four first characters select the tables. |
| `hydrometeorMoments` | <span style='color:green'>optional</span> | Dictionary `{'hydrometeor': number of moments}` (config: `hydrometeors_moments`). |
| `cloud_water_over` | <span style='color:green'>optional</span> | `'land'` or `'sea'`: cloud water tables `cl` or `cs`. |
| `subDomain` | <span style='color:green'>optional</span> | `[lon_min,lon_max,lat_min,lat_max]` (real case), `[i_min,i_max,j_min,j_max]` (idealized case) or `None` (all the points). |
| `mixed_phase_parametrization` | <span style='color:green'>optional</span> | `'T_pos'`, `'Fw_pos'` or `'Fw_posg'` (see [Simulation configuration](https://github.com/UMR-CNRM/operadar/wiki/Simulation-configuration)). |
| `save_netcdf_single_hydrometeor` | <span style='color:green'>optional</span> | Also save one file per hydrometeor type. |
| `dpol2add` | <span style='color:green'>optional</span> | List of variables among `['Zh','Zdr','Kdp','Rhohv','Ah','Av']`. |
| `scattering_method` | <span style='color:green'>optional</span> | `'Tmatrix'` (only method currently used). |
| `compute_attenuation` | <span style='color:green'>optional</span> | Compute `Zh_att`, attenuated reflectivity of a vertically pointing radar (needs `'Ah'` in `dpol2add`). |
| `radar_altitude` | <span style='color:green'>optional</span> | Altitude of the vertically pointing radar (m): 0 at the ground, or the aircraft altitude. |
| `radar_band` | <span style='color:green'>optional</span> | `'L'`, `'S'`, `'C'`, `'X'`, `'Ku'`, `'K'`, `'Ka'` or `'W'`. |
| `radarloc` | <span style='color:green'>optional</span> | Radar location for the beam geometry (MesoNH only): `'center'`, `[lat_radar,lon_radar]` or `None`. |
| `distmax_rad` | <span style='color:green'>optional</span> | Maximum range (m) when `radarloc` is set. |
| `cnst_angle` | <span style='color:green'>optional</span> | Elevation (°) when `radarloc=None`: 0 for a horizontal beam, 90 for a vertically pointing radar at each grid point. |
| `get_more_details` | <span style='color:green'>optional</span> | Print more details. Defaults to `False`. |
| `append_in_file` | <span style='color:green'>optional</span> | Write the fields into the input AROME `.fa` file instead of a netCDF file. Defaults to `False`. |

The sensitivity options (`sensitivity_threshold`, `Zmin_1km`) are only set in the configuration file.

# <a name="2gofurther"></a>To go further
- `radar_band` and `read_tables`: the lookup tables are read for one radar band. As long as the band stays the
  same, pass `read_tables` and `tables_content` returned by the previous call. To change the band, set
  `read_tables=True` again:
   ```python
   import itertools
   fname_list = ['historic.arome.franmg-01km30+0006:00.fa', 'historic.arome.franmg-01km30+0006:05.fa']
   for band in ['C', 'S', 'X']:
       read_tables, tables = True, {}
       for fname in fname_list:
           read_tables, tables = operadar(filename=fname, configuration=config, radar_band=band,
                                          read_tables=read_tables, tables_content=tables,
                                          out_dir_path=f'./operadarFiles/{band}band/')
   ```
- `in_dir_path`: files from different folders, same band:
   ```python
   read_tables, tables = True, {}
   for fpath, fname in itertools.product(['xp_arome/GOV2/', 'xp_arome/GOV5/'], fname_list):
       read_tables, tables = operadar(filename=fname, configuration=config, in_dir_path=fpath,
                                      read_tables=read_tables, tables_content=tables)
   ```
