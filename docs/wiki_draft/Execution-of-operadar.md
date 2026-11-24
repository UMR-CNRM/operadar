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
>>> $ python3 -m operadar --help

usage: python3 -m operadar [-h] [--append] [--verbose] filename config

OPERADAR

positional arguments:
  filename    Only the name of the input file, not the entire path. Please use the config file to provide the path to access the file.
  config      Name of the configuration file. Before running the code, you need to create a configuration file in ./configFiles/ based on the template provided.

options:
  -h, --help  show this help message and exit
  --append    Append the computed dual-polarimetric fields into the input file (default value: False).
  --verbose   If argument provided, show more details if activated (default value: False).
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
`operadar()` can be called inside any python program. This method is best adapted to users who want to work with multiple files or with varying parameters/configuration across file(s). Here is a way of doing so.

1) In a novel python script, import `operadar` function :
   ```python
   from operadar.forward_operator import operadar
   ```

2) You will need to load the configuration file before calling operadar so the module can access or update the settings. 
   ```python
   from operadar.load_config_file import load_configuration_file
   ```

3) Let's pretend you are working with two different study cases, that requires to read two different config files, and two different microphysics. Then : 
   ```python
   cases = [('conv','2022-06-20 03:00','2022-06-20 11:00'),
            ('strat','2024-11-19 06:00','2024-11-19 22:00'),
           ]
   microphysics_schemes = ['ICE3','LIMA']
   
   for case_type,begin,end in cases :
   for micro in microphysics_schemes :
       # Read tables and loop over the datetime (again) when changing microphysics and case stuy     
       read_tables = True 
       dict_tables = {}
       ech = pd.to_datetime(begin, format="%Y-%m-%d %H:%M")
       end = pd.to_datetime(end, format="%Y-%m-%d %H:%M")
       # Reload the corresponding configuration in the python environment
       config = load_configuration_file(f'conf_{micro}_{case_type}.py')
       # Execute operadar for each arome file
       while ech <= end :
           time = ech.strftime('%H:%M')
           print('\n','------------------------------------------',case_type,micro,time,'------------------------------------------')
           fname = f'historic.arome.franmg-01km30+00{time}.fa'
           read_tables, dict_tables = operadar(filename=fname,
                                               configuration=config,
                                               read_tables=read_tables,
                                               tables_content = dict_tables,
                                               get_more_details=False,
                                              )
           ech += dt.timedelta(minutes=5)
      ```
    Please note that **you can override all arguments** of the configuration file if they are provided as input arguments of `operadar`.
# `operadar()` arguments
This software has been designed so the user can loop on multiple files and eventually with varying configurations over files that share common config parameters. Thus, all the configuration parameters can be overwritten when calling `operadar()`. When arguments are optional, the default value is always red in the last copied configuration file :

| Parameter      | Status | Description |
| -----------    | ----------- | ----------- |
| `filename`     | <span style='color:red'>mandatory</span>  | Only the name of the file. |
| `modelname`    |<span style='color:green'>optional</span>| Can be either `'Arome'`, `'MesoNH'` or `'WRF'`. |
| `read_tables` |<span style='color:green'>optional</span>| Option to save computing time within a loop. Needs to always be `True` for the first iteration. Defaults to `True`. |
| `in_dir_path`  |<span style='color:green'>optional</span>| Path where the input files are stored. |
| `out_dir_path` |<span style='color:green'>optional</span>| Path to store the output files. |
| `tables_path` |<span style='color:green'>optional</span>| Path where the Tmatrix lookup tables are stored. |
| `microphysics_scheme` |<span style='color:green'>optional</span>| Can be `ICE3`, `ICE4` or `LIMA` + a name extension (e.g. `LIMA_noHail` or `ICE3_CIBU_moins`), which is optional. Please note that only the three first characters are used to select the right scheme in the lookup tables. Then, the microphysics and correct computations are handles with the `hydrometeorMoments` argument. |
| `hydrometeorMoments` |<span style='color:green'>optional</span>| Dictionary of form `{'hydrometeor_key' : number of moment}` corresponding to the microphysics scheme. |
| `radar_band`   |<span style='color:green'>optional</span>| Available bands : `'C'`, `'X'`, `'S'`, `'W'`, `'L'`, `'Ka'`, or `'K'` |
| `radarloc`     |<span style='color:green'>optional</span>| Location of the radar to emulate radar geometry. Either `'center'` (i.e. center of the grid) or a `[lat_radar,lon_radar]` coordinate. |
| `distmax_rad`  |<span style='color:green'>optional</span>| Maximum radar radius to compute pseudo-observations. |
|`tables_content`|<span style='color:green'>optional</span>| Dictionary containing the lookup tables parameters. Used to pass the dictionary throughout loop iterations. Please, also read the [To go further](#2gofurther) section. |
| `mixed_phase_parametrization` |<span style='color:green'>optional</span>| Can be either :<ul><li>`'T_pos'` : the species content is transferred to the melting species only at positive temperatures.</li><li>`'Fw_pos'` : the rain and graupel content are emptied and transferred into the wet graupel content within the melting layer.</li></ul><ul><li>`'Fw_posg'` : only the graupel content is emptied and transferred to the wet graupel content within the melting layer.</li></ul>  |
|`subDomain` |<span style='color:green'>optional</span>|This argument can be used to reduce the size of the output file, so you can chunk the output into smaller domain and thus lower the size of the outputs. Or, you want to work on different areas over the same file. Should be defined like `[lon_min,lon_max,lat_min,lat_max]` or set to `None` to work on all grid points.|
|`get_more_details`|<span style='color:green'>optional</span>| Boolean to print more details and computation steps. |



# <a name="2gofurther"></a>To go further
- `radar_band` and `read_tables` : the Tmatrix lookup tables are red once, for the radar band set in the config file, at the beginning of the code. To save computation time, it is not necessary to read the tables as long as the radar band remains the same, and thus, `read_tables` is automatically set to `False` after the first iteration that produces an output file. If you want to change the radar band over the iterations, you must also set `read_tables=True`. Based on the tutorial :
   ```python
      import itertools

      fname_list = ['historic.arome.franmg-01km30+0006:00.fa',
                    'historic.arome.franmg-01km30+0006:05.fa',
                   ]
      radar_band = ['C','S','X']
      read_tables = [True]

      # Creating the unique combinations
      combinations = list(itertools.product(fname_list,radar_band,read_tables))
   ```
   will produce
   ```
   >>> $ python3 operadar_multi.py
   historic.arome.franmg-01km30+0006:00.fa C True
   historic.arome.franmg-01km30+0006:00.fa S True
   historic.arome.franmg-01km30+0006:00.fa X True
   historic.arome.franmg-01km30+0006:05.fa C True
   historic.arome.franmg-01km30+0006:05.fa S True
   historic.arome.franmg-01km30+0006:05.fa X True
   ```
   and can then be used like 
   ```python
      for fname,band,read_tmat in combinations :
         read_tmat,dict_Tmatrix = operadar(filename=fname,
                                           read_tables=read_tmat,
                                           radar_band=band,
                                           get_more_details=True,
                                           )

   ```
- `out_dir_path` and `in_dir_path` : similarly, one may want to store the output files in different folders depending on the radar band...
   ```python
      for band in radar_band :
         read_tables=True
         dict_Tmatrix={}
         outPath = f'/home/my_path/output/{band}band_simus/'
         for fname in fname_list :
            read_tables,dict_Tmatrix = operadar(filename=fname,
                                                 read_tables=read_tmat,
                                                 radar_band=band,
                                                 Tmatrix_params=dict_Tmatrix,
                                                 out_dir_path=outPath,
                                                )

   ```
   ... or provide files from different folders :
   ```python
      fpath_list = ['xp_arome/GOV2/',
                    'xp_arome/GOV5/',
                   ]
      fname_list = ['historic.arome.franmg-01km30+0006:00.fa',
                    'historic.arome.franmg-01km30+0006:05.fa',
                   ]
      combinations = list(itertools.product(fpath_list,fname_list))
      band='X'
      read_tables=True
      dict_Tmatrix={}
      for fpath,fname in combinations :
         read_tables,dict_Tmatrix = operadar(filename=fname,
                                             read_tables=read_tmat,
                                             radar_band=band,
                                             Tmatrix_params=dict_Tmatrix,
                                             in_dir_path=fpath,
                                             )

   ```

