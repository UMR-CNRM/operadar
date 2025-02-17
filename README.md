# Radar forward operator (operadar)
This radar forward operator is developped by the researchers of the National Centre for Meteorological Research (CNRM).

The software is designed to compute synthetic dual-polarization variables (Z<sub>H</sub>, Z<sub>DR</sub>, K<sub>DP</sub>, Rho<sub>HV</sub>) in a 3D grid.
It is adapted to read and manipulate AROME and MesoNH model files and currently requires a repository of T-matrix lookup tables.

### INPUT
* Model file (`.fa` for AROME files or `.nc` for MesoNH)
* Configuration file
* T-matrix lookup tables (created beforehand with the T-matrix generation code). For colleagues at the CNRM, the tables can be found on belenos at `/home/augros/TmatCoefInt_SCXW/`

### OUTPUT : netcdf file with
* X/Y grid horizontal coordinates (1D)
* model pressure levels vertical coordinates (1D)
* latitude/longitude fields (2D)
* dual-pol variables : $Z_{H}$ , $Z_{DR}$ , $K_{DP}$ , $\rho_{HV}$ (3D)
* temperature, altitude, contents and concentrations fields (3D)
* dataset attributes such as the horizontal resolution, the options used to compute the mixed-phase or simulate the radar data.

# Installation
## On Belenos: 
create a .gitconfig file (on your home directory) with:

```
[http]
       sslVerify = false
```

## As a package (not sure yet of the procedure)
In a virtual environment (recommended) or in your base environment
```
pip install git+https://github.com/UMR-CNRM/operadar.git
```


# How to run ( WIKI )

## Configuration file
The user must set up, beforehand, a configuration file based on the template provided `./configFile/conf_template.py`. Please, do not modify the template directly. Instead, make a copy of the file and name it differently. 

## Quick execution in a terminal (suitable for 1 file at a time)
The script `exec_operadar.sh` wraps the execution of the Python code. To show the help :
```
>>> $ ./exec_operadar.sh -h

----------------------------------------------------------------
 This radar forward operator is developped at the CNRM, France. 
----------------------------------------------------------------

Usage: ./exec_operadar.sh -f FILENAME -c CONFIG [--verbose]

  -f FILENAME : Only the filename. Please use the config file to provide the path to access the file.
  -c CONFIG   : Before running the code, you need to create a configuration file in ./configFiles/ based on the template provided.
  --verbose   : Optional, show more details if activated (default value: False)

```
## Inside another Python program (ideal for multiple file)
The user can work with multiple files, if desired. The paths, output directories and other parameters must be specified in the configuration file. Here is a way of doing so.
### Simple tutorial
1) In another python code (let say `operadar_multi.py`), call `operdar()` function :
   ```python
   from operadar.forward_operator import operadar
   ```

2) The configuration file need to be copied beforehand under the name and location `operadar/operadar_conf.py`, so all files of the module can access the settings.
   ```python
   import os
   configFileName = 'my_conf_file' # based on the template provided
   os.system(f'cp ./configFiles/{configFileName}.py operadar/operadar_conf.py')  
   ```

3) The user must then create a list of file. It can be done by hand or within a loop if the number of file is huge.
   ```python
   fname_list = ['historic.arome.franmg-01km30+0006:00.fa',
                 'historic.arome.franmg-01km30+0006:05.fa',
                 'historic.arome.franmg-01km30+0006:10.fa',
                 'historic.arome.franmg-01km30+0006:15.fa',
                 'historic.arome.franmg-01km30+0006:20.fa',
                 'historic.arome.franmg-01km30+0006:25.fa',
                ]
   ```

4) You are ready to call `operadar()` ! Don't forget to initiate the arguments passed to the function :
   ```python
      # Parameters initialization
      read_tmatrix = True
      dict_Tmatrix = {}

      for fname in fname_list :
         read_tmatrix,dict_Tmatrix = operadar(filename=fname,
                                              read_tmatrix=read_tmatrix,
                                              Tmatrix_params=dict_Tmatrix,
                                              get_more_details=True
                                              )
      ```

## To go further
This code has been designed so the user can loop on multiple files and eventually with varying configurations over files that share common config parameters. Thus, some of the parameters can be overwritten when calling `operadar()` :
- `radar_band` and `read_tmatrix` : the Tmatrix lookup tables are red once, for the radar band set in the config file, at the beginning of the code. To save computation time, it is not necessary to read the tables as long as the radar band remains the same, and thus, `read_tmatrix` is automatically set to `False` after the first iteration that produces an output file. If you want to change the radar band over the iterations, you must also set `read_tmatrix=True`. Based on the tutorial :
   ```python
      import itertools

      fname_list = ['historic.arome.franmg-01km30+0006:00.fa',
                    'historic.arome.franmg-01km30+0006:05.fa',
                   ]
      radar_band = ['C','S','X']
      read_tmatrix = [True]

      # Creating the unique combinations
      combinations = list(itertools.product(fname_list,radar_band,read_tmatrix))
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
                                           read_tmatrix=read_tmat,
                                           radar_band=band,
                                           get_more_details=True,
                                           )

   ```
- `out_dir_path` : similarly, one may want to store the output files in different folders depending on the radar band for exemple.
   ```python
      for band in radar_band :
         read_tmatrix=True
         outPath = f'/home/my_path/output/{band}band_simus/'
         for fname in fname_list :
            read_tmatrix,dict_Tmatrix = operadar(filename=fname,
                                                 read_tmatrix=read_tmat,
                                                 radar_band=band,
                                                 out_dir_path=outPath,
                                                )

   ```
- `subDomain` : This argument can be used to reduce the size of the output file, so you can chunk the output into smaller domain and thus lower the size of the outputs. Or, you want to work on different areas over the same file. 



# To Do
- [ ] find a way to read the ice concentration in AROME file ?
- [ ] add the radar geometry option (with elevations and beam filtering with gaussian)
- [ ] link this project with tmatrix DPOLSIMUL git (which produces the required tables !)
- [ ] put Tmatrix tables online so that every one can get them
- [ ] put one MesoNh and one AROME example test files


## TODO CLOE
- [x] mettre au propre README
- [ ] maj toutes les docstrings
- [ ] configurer fichier gitignore : operad_conf, tmp obj, etc
- [ ] read_mesonh : changer les sorties (enlever CC et CCI pour Nc)
- [ ] calcul explicite de la concentration pour les espèces 1 moment


# Contributing
If you wish to contribute to the project, first, fork the code to create your own copy of the project (see https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project#creating-your-own-copy-of-a-project).
<br>NOTE : If you want to be part of the main developpers (ask clotilde.augros@meteo.fr)

Then, clone the repository in a dedicated folder.
```
mkdir my_folder
cd my_folder
git clone https://github.com/UMR-CNRM/operadar.git
```

Before making changes to the project, you should create a new branch and check it out. Try to be self explanatory for the name of the branch.
```
git branch my_branch
git checkout my_branch
```

In the installation folder, install operadar as a package with
```
pip install -e .
```

To ask for integration of your modifications, please open a pull request following https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project#making-a-pull-request.

# License
??