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
## On a server 
NOTE : It is easier to work with virtual environments. Here we present the installation with anaconda, but using any module to create a virtual environment will do the trick. The installation can also be realised with `pip`.
### Additional steps (only for HPCs at Météo-France)
Create a .gitconfig file (on your home directory) with :

```
[http]
       sslVerify = false
```
And deactivate the pre-installed conda environment :
```
$ condam
```
### Procedure
1) Create a virtual environment. Let's named it `myenv`. `operadar` requires at least Python 3.11 :
   ```
   $ conda create -n myenv python=3.11
   ```
2) Activate the environment with
   ```
   $ conda activate myenv
   ```
3) Install the following dependencies...
   ```
   $ conda install pandas xarray netCDF4 numpy h5py
   ```
   ...as well as the [EPyGrAM package](https://github.com/UMR-CNRM/EPyGrAM) (only available with `pip`).
   ```
   $ pip install epygram
   ```
   CAREFUL : You may need to install `matplotlib` and `cartopy`, but it seems to generate incompatibilities with `numpy`. To solve the issue, the currently proposed solution is either to remove the plot directory of `operadar`, or update `numpy` to its last version. This is a temporary solution and the issue will be fixed later.
3) Now, clone the `operadar` repository
   ```
   $ git clone https://github.com/UMR-CNRM/operadar.git
   $ cd operadar
   >>> $ tree -L 1
   .
   ├── configFiles
   ├── exec_operadar.sh
   ├── operadar
   ├── pyscattering
   ├── README.md
   ├── requirements.txt
   └── setup.py
   ```
4) Installing the package with `pip install -e .` should create a `operadar.egg-info` directory.
   ```
   >>> $ tree -L 1
   .
   ├── configFiles
   ├── exec_operadar.sh
   ├── operadar
   ├── operadar.egg-info
   ├── pyscattering
   ├── README.md
   ├── requirements.txt
   └── setup.py
   ```

## As a package (not sure yet of the procedure)
In a virtual environment (recommended) or in your base environment
```
pip install git+https://github.com/UMR-CNRM/operadar.git
```


# How to run
The lookup tables should be generated BEFORE running `operadar` (unless you have a direct access to the tables).

## Generation of the Tmatrix lookup tables
To work properly, the forward operator needs to have access to the Tmatrix lookup tables. The codes that generate the tables are in Fortran, and requires the installation of a Fortran compiler.
### 1. Fortran compiler installation
For this installation, we will use the `lapack` library.
1) First, clone the last version of `lapack` ([check here](https://github.com/Reference-LAPACK/lapack)).
   ```
   cd tmatrix_generator/src
   mkdir my_lapack
   cd my_lapack
   git clone https://github.com/Reference-LAPACK/lapack.git
   ```
   **NOTE** : You can install `lapack` wherever you like but need to check that the path specified in the `makefile` corresponds to your installation.
2) To compile the library, you need a `make.inc` file. Multiple examples of such a configuration file are available in the `INSTALL` folder, or you can simply use the `make.inc.example`. Still in the `lapack/` directory do :
   ```
   cp make.inc.example make.inc
   make
   ```
   When the installation is completed, you should see on your terminal :
   ```
               -->   LAPACK TESTING SUMMARY  <--
      Processing LAPACK Testing output found in the TESTING directory
   SUMMARY             	nb test run numerical error   other error  
   ================   	===========	=================	================  
   REAL             	   1569648		0	(0.000%)	      0	(0.000%)	
   DOUBLE PRECISION	   1570470		0	(0.000%)	      0	(0.000%)	
   COMPLEX          	   980455		0	(0.000%)	      0	(0.000%)	
   COMPLEX16         	1030797		0	(0.000%)	      0	(0.000%)	

   --> ALL PRECISIONS	5151370		0	(0.000%)	      0	(0.000%)
   ``` 
3) Go back to `src/` and open both `makefile` and `makeTmatInt`. Make sure the `LIBS` variable points to the previously installed version of `lapack` :
   ```shell
   LIBS = -L ./my_lapack -llapack -lrefblas
   ```
### 2. Creation of the lookup tables for discrete hydrometeor diameters
This is a mandatory step. These tables will be used later for the integration of the diameters over the particle size distribution (PSD), in the second Fortran code.
1) Make sure you are in the `tmatrix_generator/` directory.
2) The parametrization for each hydrometeor type (density, axis ratio, etc) is defined in a txt file that you may need to modify accordingly to your needs : `tmatrix_generator/param/TmatParam_(radarBand)(hydrometeorType)`
3) To specify which radar band and hydrometeor types the lookup tables should be created, open `tmatrix_generator/src/Tmatrix.f`, search for `DO idtype` and `DO idband`. You can loop over all the band and hydrometeor types (`DO idxxx=minNumber,maxNumber`), or just a few combination of them, or only a single band and/or hydrometeor type (`DO idxxx=sameNumber,sameNumber`). In the following example, tables will be generated for C band only and for rain, snow, pristine ice, graupel, and wet graupel hydrometeor types.
   ```fortran
   !=======  Loop over radar frequency bands
      DO idbande=2,2
          IF (idbande .EQ. 1) bande='S'
          IF (idbande .EQ. 2) bande='C'
          IF (idbande .EQ. 3) bande='X' ! LAM=31.9 mm
          IF (idbande .EQ. 4) bande='W' ! LAM=3.15 mm (Rasta) 
	      IF (idbande .EQ. 5) bande='K' ! LAM=8.40 mm Ka (C3IEL)
      
   !======= Loop over hydrometeor types
      DO idtype=1,5!3,3 !2,2 !1,1 !2,2
         IF (idtype .EQ. 1) typeh='rr'
         IF (idtype .EQ. 2) typeh='ss'
         IF (idtype .EQ. 3) typeh='ii'
         IF (idtype .EQ. 4) typeh='gg'
         IF (idtype .EQ. 5) typeh='wg' !wet graupel
         IF (idtype .EQ. 6) typeh='tt' ! MARY: cloud water over land
         IF (idtype .EQ. 7) typeh='mm' ! MARY: cloud water over sea
         IF (idtype .EQ. 8) typeh='hh' ! hail
         IF (idtype .EQ. 9) typeh='wh' ! wet hail
         IF (idtype .EQ. 10) typeh='ws' !wet snow
   ```
4) To compile the `Tmatrix.f` code and create the `Tmat` executable, just do :
   ```
   make
   ```
   It should take a few moment to compile the code. Once it's done, check that `Tmat` is executable, and if not do :
   ```
   chmod u+x Tmat
   ```

5) You can now generate the lookup tables (for discrete diameters) for different hydrometeor types and different radar band ! 
   ```
   Tmat
   ```
   Each time you modify `Tmatrix.f`, you will need to update the `Tmat` executable with :
   ```
   make clean
   make
   ```
### 3. Creation of the lookup tables for discrete hydrometeor contents (integration over the PSD)
The previous step is mandatory, otherwise tables won't be created.
1) Open `tmatrix_generator/src/TmatInt.f90` and search for `DO idtype` and `DO idband`.
2) As in step 2.3, modify the loops to match those of `Tmatrix.f`. Given the previous example, we will have to loop over C band only and to loop over rain, snow, pristine ice, graupel, and wet graupel hydrometeor types.
3) Search for the `CCLOUD` variable and update the wanted microphysics scheme for the integration over the PSD by changing the name.
   ```fortran
   CCLOUD="ICE3" ! Possible microphysics : "ICE3", "LIMA" or "LIMT"
   ```
4) No more modifications are required, you can now compile the code and create another executable (don't forget to check if it is executable) :
   ```
   make -f makeTmatInt
   ```
5) You can now generate the lookup tables for discrete hydrometeor contents ! 
   ```
   TmatInt
   ```
   Each time you modify `TmatInt.f90`, you will need to update the `TmatInt` executable with :
   ```
   make -f makeTmatInt clean
   make 
   ```

### Executable (useful for sensitivity tests)
An executable is provided to facilitate table generation for a single hydrometeor type but different parameters combinations. You still need to complete step 1.5 (select one or multiple radar band but **only one** hydrometeor type), and step 7.
```
>>> $ ./execTmat.sh --help
----------------------------------------
 Executable to create the lookup tables 
----------------------------------------

Usage: ./execTmat.sh -hydro HYDRO -af ARfunc -av ARvalue -c CANTING -dsty DSTYfunc -riming RIMING -diel DIELfunc

  -hydro HYDRO   : rr, ii, gg, ss, tt, wg, hh, wh 
  -af ARfunc     : AUds, CNST, BR02, RYdg, RYwg
  -av ARvalue    : any value.
  -c CANTING     : any value.
  -dsty DSTYfunc : BR07, RHOX
  -riming RIMING : any value starting from 1 (1=unrimed)
  -diel DIELfunc : Liebe91, RY19dry, LBwetgr, MGwMA08
```

## Forward operator configuration file
The user must set up, beforehand, a configuration file based on the template provided `./configFile/template.py`. Please, do not modify the template directly. Instead, make a copy of the file and name it differently. 

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
The user can work with multiple files, if desired. Here is a way of doing so.
### Simple tutorial
1) In another python code (let say `operadar_multi.py`), call `operadar()` function :
   ```python
   from operadar.forward_operator import operadar
   ```

2) The configuration file need to be copied beforehand under the name and location `operadar/operadar_conf.py`, so the module can access the settings.
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

### `operadar()` parameters
This code has been designed so the user can loop on multiple files and eventually with varying configurations over files that share common config parameters. Thus, all the configuration parameters can be overwritten when calling `operadar()`. When arguments are optional, the default value is always read in the last copied configuration file :

| Parameter      | Status | Description |
| -----------    | ----------- | ----------- |
| `filename`     | <span style='color:red'>mandatory</span>  | Only the name of the file. |
| `modelname`    |<span style='color:green'>optional</span>| Can be either `'Arome'` or `'MesoNH'`. |
| `read_tmatrix` |<span style='color:green'>optional</span>| Option to save computing time within a loop. Needs to always be `True` for the first iteration. Defaults to `True`. |
| `in_dir_path`  |<span style='color:green'>optional</span>| Path where the input files are stored. |
| `out_dir_path` |<span style='color:green'>optional</span>| Path to store the output files. |
| `tmatrix_path` |<span style='color:green'>optional</span>| Path where the Tmatrix lookup tables are stored. |
| `microphysics_scheme` |<span style='color:green'>optional</span>| Can be `ICE3`, `ICE4` or `LIMA` + a name extension (e.g. `LIMA_noHail` or `ICE3_CIBU_moins`), which is optional. Please note that only the three first characters are used to select the right scheme in the lookup tables. Then, the microphysics and correct computations are handles with the `hydrometeorMoments` argument. |
| `hydrometeorMoments` |<span style='color:green'>optional</span>| Dictionary of form `{'hydrometeor_key' : number of moment}` corresponding to the microphysics scheme. |
| `radar_band`   |<span style='color:green'>optional</span>| Available bands : `'C'`, `'X'`, `'S'`, `'W'` or `'K'` |
| `radarloc`     |<span style='color:green'>optional</span>| Location of the radar to emulate radar geometry. Either `'center'` (i.e. center of the grid) or a `[lat_radar,lon_radar]` coordinate. |
| `distmax_rad`  |<span style='color:green'>optional</span>| Maximum radar radius to compute pseudo-observations. |
|`Tmatrix_params`|<span style='color:green'>optional</span>| Dictionary containing the Tmatrix tables parameters. Used to pass the dictionary throughout loop iterations. Please, also read the [To go further](#2gofurther) section. |
| `mixed_phase_parametrization` |<span style='color:green'>optional</span>| Can be either :<ul><li>`'T_pos'` : the species content is transferred to the melting species only at positive temperatures.</li><li>`'Fw_pos'` : the rain and graupel content are emptied and transferred into the wet graupel content within the melting layer.</li></ul><ul><li>`'Fw_posg'` : only the graupel content is emptied and transferred to the wet graupel content within the melting layer.</li></ul>  |
|`subDomain` |<span style='color:green'>optional</span>|This argument can be used to reduce the size of the output file, so you can chunk the output into smaller domain and thus lower the size of the outputs. Or, you want to work on different areas over the same file. Should be defined like `[lon_min,lon_max,lat_min,lat_max]` or set to `None` to work on all grid points.|
|`get_more_details`|<span style='color:green'>optional</span>| Boolean to print more details and computation steps. |



## <a name="2gofurther"></a>To go further
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
- `out_dir_path` and `in_dir_path` : similarly, one may want to store the output files in different folders depending on the radar band...
   ```python
      for band in radar_band :
         read_tmatrix=True
         dict_Tmatrix={}
         outPath = f'/home/my_path/output/{band}band_simus/'
         for fname in fname_list :
            read_tmatrix,dict_Tmatrix = operadar(filename=fname,
                                                 read_tmatrix=read_tmat,
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
      read_tmatrix=True
      dict_Tmatrix={}
      for fpath,fname in combinations :
         read_tmatrix,dict_Tmatrix = operadar(filename=fname,
                                             read_tmatrix=read_tmat,
                                             radar_band=band,
                                             Tmatrix_params=dict_Tmatrix,
                                             in_dir_path=fpath,
                                             )

   ```


# To Do
- [ ] find a way to read the ice concentration in AROME file ?
- [ ] add the radar geometry option (with elevations and beam filtering with gaussian)
- [ ] link this project with tmatrix DPOLSIMUL git (which produces the required tables !)
- [ ] put Tmatrix tables online so that every one can get them
- [ ] put one MesoNh and one AROME example test files


## To do Cloe
- [x] mettre au propre README
- [x] maj toutes les docstrings
- [x] configurer fichier gitignore : operad_conf, tmp obj, etc
- [ ] read_mesonh : changer les sorties (enlever CC et CCI pour Nc)
- [ ] calcul explicite de la concentration pour les espèces 1 moment
- [ ] check radar band consistency if operdar used within a loop
- [ ] add a way to reinject the compute dpol fields into the input file


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
