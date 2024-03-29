# operadar
Computes dual-pol variables (Z<sub>H</sub>, Z<sub>DR</sub>, K<sub>DP</sub>, Rho<sub>HV</sub>) in the 3D model grid for Arome or MesoNH model using existing T-matrix tables
* INPUT  : T-matrix tables and model file (Arome .fa or MesoNH netcdf)
  
Tmatrix tables directory at CNRM on belenos:  /home/augros/TmatCoefInt_SCXW/
* OUTPUT : netcdf file with $lat/lon$ (or $X/Y$) + $Z_{H}$ , $Z_{DR}$ , $K_{DP}$ , $\rho_{HV}$ , $T$, and altitude for each model level

## How to get the code
For belenos: 
create a .gitconfig file (on your home directory) with:

```
[http]
       sslVerify = false
```

* If you only need to run the code:
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


## How to run operadar
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
