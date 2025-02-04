# Radar forward operator (operad)
This radar forward operator is developped by the researchers of the GMME/PRECIP team at the National Centre for Meteorological Research (CNRM).

This software is designed to compute synthetic dual-polarization variables (Z<sub>H</sub>, Z<sub>DR</sub>, K<sub>DP</sub>, Rho<sub>HV</sub>) in a 3D grid.
It is adapted to read and manipulate AROME and MesoNH model files and requires a repository of T-matrix lookup tables.

* INPUT
  * T-matrix lookup tables (created beforehand with the T-matrix generation code). For colleagues at the CNRM, the tables can be found on belenos at `/home/augros/TmatCoefInt_SCXW/`
  * Model file (`.fa` for AROME files or `.nc` for MesoNH)
* OUTPUT : netcdf file with
  * latitude/longitude fields (2D), and X/Y 1D coordinates
  * model level coordinates (1D)
  * dual-pol variables ($Z_{H}$ , $Z_{DR}$ , $K_{DP}$ , $\rho_{HV}$)
  * temperature, altitude, contents and concentrations fields

# Installation
## On Belenos: 
create a .gitconfig file (on your home directory) with:

```
[http]
       sslVerify = false
```

## As a package (à MAJ)
```
git clone https://github.com/UMR-CNRM/operadar.git
```

* If you want to suggest slight modifications
  1. fork the code (select: Fork, create a new fork)
  2. get the name of your new fork repository by clicking on the Code button (copy the https url) 
     e. g. :
     ```
     git clone https://github.com/augros/operadar.git 
     ```
* If you want to be part of the main developpers (ask clotilde.augros@meteo.fr)


# How to run
`operad` is intended to work on one model file and takes only one argument. The user must set up, beforehand, a configuration file based on the examples furnished in the `./src/operad/configFile/` repository.

## Quick execution in a terminal (suitable for 1 file)
Please, run `exec_operad.sh` in a terminal with the following arguments and in this order :
1) the file path (extension included)
2) the datetime in YYYYmmddHHMM format
3) the configuration file
```
>>> $ ./exec_operad file_path config_file_name
```
## With `main.py` program
The user can work with multiple files, if desired. The paths, output directories and other parameters must be specified in the configuration file.
1) **Working with a unique file** <br> 
   In `main.py` just specify the file name and the configuration file arguments.
   ```python
   operad(filename='my_file_name.extension', config='my_config_file.py')
   ```
2) **Working with multiple files at a time** <br> 
   `filename` argument also accepts `list` type. The user must create a list of files before calling `operad()`. If the list of file is huge, parallelization can be activated with `parallel=True` (default to `False`).
<!---
PENSER À AJOUTER DANS OPERAD :
+ ARGUMENT SWITCH POUR NE PAS REPETER TMATRICE -> The switch argument will be equal to False after one iteration, so the T-matrix tables are read once. 
+ VERIF DU TYPE DE FILENAME : SI = LISTE ALORS BOUCLE, SINON LANCEMENT OPERAD DIRECT
+ ARGUMENT PARALLEL POUR PARALLELISER LE TRAITEMENT DES FICHIERS
-->
   ```python
   my_list_of_file_names = ['toto1.extension','toto2.extension','toto3.extension']
   operad(filename=my_list_of_file_names, config = 'my_config_file.py')
   ```
   

## Configuration file

1) Create a study case file in ./study_cases/  (eg: CORSE_Arome.csv or CORSE_MesoNH.csv)
   ==> you copy one of the examples available that best fits to your simulation

 start_time;end_time;radar_id_list;radar_band;run_model;latmin;latmax;lonmin;lonmax;

* start_time;end_time => period of your study (each time corresponds to 1 Arome or MesoNH file, the time step is defined in the configuration file below)
* radar_id_list => not used in this version 
* radar_band => to specify which tables to read (S, C or X band)
* run_model => to specify at what time the simulation starts (for MesoNH) or the Arome run (00, 03 ...)
* latmin;latmax;lonmin;lonmax => to restrict a specific region where to compute the radar variables (only for Arome)
        
2) Create a configuration file in ./configFiles/ (eg: conf_AROME_ICE3_CORSEbe.py or conf_MesoNH_ICE3_CORSEbe.py)

You need to change:
<!-- A MODIFIER :
Il n'y a plus htypes marqué en dur dans le fichier de config, c'est intégré au programme
-->
* htypes_model and list_types_tot

   - ICE3 or LIMA (without hail) :
     - `htypes_model   = ['vv','cc','rr','ii','ss','gg']` # model's hydrometeors related variables
     - `list_types_tot = ['rr','ii','ss','gg','wg']`    # idem **+ wet hydrometeors (computed)**
   - ICE4 or LIMA **with** hail :
     - `htypes_model   = ['vv','cc','rr','ii','ss','gg','hh']`
     - `list_types_tot = ['rr','ii','ss','gg','wg','hh','wh']`
* directories and file name options (T-matrix table and model file directories)
* forward operator options:
  save_netcdf = True (keep netcdf format)
  step = dt.timedelta(hours=1) or dt.timedelta(minutes=15) or dt.timedelta(minutes=5) ... depending of the time step between 2 consecutive model files in your simulation
* radar options:
   * distmax_rad (m) => maximum range (m) until which variables are computed (with the default option corresponding to a radar located in the center of the domain)
   * radarloc="center" (center of the domain by default, no other option yet)

3) Run with exec_operad.sh
   
takes 4 arguments in this order :
- 1 Arome or MesoNH
- 2 date into the yyyymmdd format or "all"
- 3 microphysics scheme name in capital letter (ICE3, LIMA, ICE4, LIMAAG)
- 4 Config file specifying directories and forward operator options

  Before running, select the option with or without nohup in exec_operad.sh (directly with python3 -i if you need to debub, or with nohup if the script is running well)

Examples :

   `>>> ./exec_operad.sh MesoNH 20220818 ICE3 conf_MesoNH_ICE3_CORSEbe.py`
   
   `>>> ./exec_operad.sh MesoNH 20220818 LIMA conf_MesoNH_LIMA_CORSEbe.py`
   
   `>>> ./exec_operad.sh Arome 20220818 ICE3 conf_Arome_ICE3_CORSEbe.py`

   `>>> ./exec_operad.sh MesoNH 20220818 LIMAAG conf_MesoNH_LIMH_CORSEbe.py`

## To Do
- [ ] find a way to read the ice concentration in AROME file ?
- [ ] add the radar geometry option (with elevations and beam filtering with gaussian)
- [ ] link this project with tmatrix DPOLSIMUL git (which produces the required tables !)
- [ ] put Tmatrix tables online so that every one can get them
- [ ] put one MesoNh and one AROME example test files


# TODO CLOE
- [] mettre au propre README
- [] configurer fichier gitignore : operad_conf, tmp obj, etc
- [] read_mesonh : changer les sorties (enlever CC et CCI pour Nc)
- [] calcul explicite de la concentration pour les espèces 1 moment