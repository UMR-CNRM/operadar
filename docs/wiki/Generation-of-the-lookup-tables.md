To work properly, the forward operator needs to have access to lookup tables. The codes that generate the tables are in Fortran, and require the installation of a Fortran compiler. There is one table per hydrometeor type and per frequency band (combinations among those available, see [here](https://github.com/UMR-CNRM/operadar/wiki#available-radar-bands)).

# 1. Fortran compiler installation
For this installation, we will use the `lapack` library.
1) First, clone the last version of `lapack` ([check here](https://github.com/Reference-LAPACK/lapack)).
   ```
   cd tables_generator/src
   mkdir libFortran
   cd libFortran
   git clone https://github.com/Reference-LAPACK/lapack.git
   ```
2) To compile the library, you need a `make.inc` file. Multiple examples of such a configuration file are available in the `INSTALL` folder, or you can simply use the `make.inc.example`. Go to the `lapack/` directory and do :
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
3) Go back to `src/` and open both `makefile_createTable` and `makefile_integrationOverPsd`. Make sure the `LIBS` variable points to the previously installed version of `lapack` :
   ```shell
   LIBS = -L ./libFortran/lapack -llapack -lrefblas
   ```
   :pencil2: NOTE : You can install `lapack` wherever you like, but you need to check that the path specified in both files corresponds to your installation.

# 2. Creation of the lookup tables

Now that a compiler is installed, in `tables_generator/src/`, compile the Fortran codes thanks to the makefiles. To do so, run these two commands :
```
make -f makefile_createTable
make -f makefile_integrationOverPsd
```
Two executable named `Tmat` and `TmatInt` should have appeared. Once it's done, check that they are both executable. If not, use `chmod u+x Tmat` and `chmod u+x TmatInt`.  
:warning: The executables are not stored in git: compile them once after cloning, and again after any change of the Fortran sources (e.g. after a `git pull`). If they are missing, `generate_tables.sh` stops at once and prints these two commands.  
These executable are launched by the shell script `generate_tables.sh` located at the root of the project. Thus, tables can be generated for different frequency bands, either with the default configuration or your personal configuration. They are stored under the `tables_generator/tables/` directory.
<!--- In `generate_tables.sh`, the creation of the tables is parallelized, and each output of the loop over the hydrometeor types is redirected to a `{band}_{hydrometeor}.log` file in the `./log/` directory. ---> 

Tables are created in a two-step process.
1) For a range of temperature, elevation, liquid water fraction (if applicable) and diameters, scattering coefficients of each hydrometeor class are computed and arranged in a first table by `Tmat`.
2) This first table is red by `TmatInt` and an integration over a given particle size distribution (PSD) is performed, so new tables are created as a function of the hydrometeor content. Here, the PSD is inherited from one of the two French microphysics schemes : ICE3 (1-moment ; [Pinty and Jabouille, 1998](http://195.83.22.22/mesonh/dir_publication/pinty_jabouille_ams_ccp1998.pdf)) or LIMA (2-moments ; [Vié et al. (2016)](https://doi.org/10.5194/gmd-9-567-2016) and [Taufour et al. (2024)](https://doi.org/10.1002/qj.3283)). Currently, tables are generated twice (once with ICE3 PSD and once with LIMA PSD).

## a. Predefined configurations
These modes create the tables for a given frequency band and for all hydrometeor types (some might be unavailable depending on the chosen band), from predefined parameter files that **should not** be modified. Two configurations are available (see [Lookup table configurations](https://github.com/UMR-CNRM/operadar/wiki/Lookup-table-configurations) for their options and their differences):
* `--David2026PhD` (all bands: L, S, C, X, Ku, K, Ka, W), the **default** configuration: `--default` is an alias for it;
* `--David2025AMT` (C band): the options of [David et al. (2025)](https://doi.org/10.5194/amt-18-3715-2025), no canting and spherical pristine ice.

Both cover elevations from 0° to 90° every 10° (ground, vertically pointing and spaceborne radars).

To generate the default tables in C band, at the root of the project, run :
```
./generate_tables.sh --default --band C -j 8
```
or
```
nohup ./generate_tables.sh --default --band C -j 8 > resuCband 2>&1 &
```
and for another band or configuration, e.g. `./generate_tables.sh --default --band W -j 8` or `./generate_tables.sh --David2025AMT --band C -j 8`.
The option `-j N` runs `Tmat` and `TmatInt` as N parallel processes (blocks of temperatures; the tables are identical to a run with `-j 1`). Use the number of cores of your machine, e.g. `-j $(nproc)`. With `-j 8`, the whole C band set (ICE3 and LIMA, 10 elevations) takes about 25 minutes for David2025AMT and 30 minutes for David2026PhD (the averaging over the orientations of `canting=orient` costs a little more); in NetCDF, the tables of one band take about 1 GB.

Tables will be stored under `tables_generator/tables/<configuration>/` (e.g. `tables/David2026PhD/`) and named after the hydrometeor type, the frequency band, and the microphysics scheme and number of moments used for the integration over the diameters. For example, `TmatCoefInt_ICE3_1M_Crr.nc` will be the table of scattering coefficients in rain, simulated by a C band radar and with the ICE3 PSD (`TmatCoefInt_LIMA_2M_Crr.nc` with the 2-moment LIMA PSD). The single-particle tables are `TmatCoefDiff_Crr.nc`, and the `Tmat` logs `Tmat_Crr.log`.

### Table format (NetCDF)
`Tmat` and `TmatInt` write text tables, which `generate_tables.sh` converts to NetCDF (`<table>.nc`, about 10 times smaller) before deleting the text files. The conversion needs a Python with `xarray` and `netCDF4` (the operadar environment): activate it before running the script, or give it with `PYTHON=/path/to/python ./generate_tables.sh ...`. The option `--text` keeps the text tables (no conversion).

In a NetCDF table, each column of the former text table is a variable on the grid of the table:
* `TmatCoefInt_*`: dimensions `(Tc, ELEV, Fw, M)` for 1-moment species, `(Tc, ELEV, N, M)` for 2-moment species and pristine ice (M in kg/m³, N in m⁻³);
* `TmatCoefDiff_*`: dimensions `(Tc, ELEV, Fw, D)` (D in mm).

Variables: `sighh`, `sigvv`, `REdeltaco`, `IMdeltaco`, `zhh`, `zdr`, `kdp`, `rhohv`, `Ah`, `Av` (T-matrix), the same with the suffix `R` (Rayleigh), `RR` (rain rate) and `Mfrac` for the integrated tables, and `Dm`, `Deq`, `Deqm`, `Deqr`, `Deqrm` for the single-particle tables. The header parameters (`Tcmin`, `Tcstep`, `ELEVmax`, `expMmin`, `LAM`, `canting`, `ARfunc`...) are global attributes. Values are stored in single precision (as computed by the Fortran programs) and compressed.

The tables can be read and plotted directly with `xarray`, e.g.
```python
import xarray as xr
ds = xr.open_dataset("tables_generator/tables/David2026PhD/TmatCoefInt_ICE3_1M_Crr.nc")
ds.zdr.sel(Tc=10, ELEV=0, Fw=0).plot(x="M", xscale="log")
ds.kdp.sel(Tc=10, Fw=0).sel(M=1e-3, method="nearest").plot()   # Kdp as a function of the elevation
```
or looked at with `ncdump -h` or `ncview`. operadar reads only the variables it needs. `operadar/read/table_io.py` converts tables between the two formats (`python operadar/read/table_io.py to-netcdf TABLE` converts a text table, e.g. from an older version; `to-text TABLE.nc` writes a single-particle table back to text, which `generate_tables.sh` does automatically when `TmatInt` needs it).

:warning: Existing tables are kept: `generate_tables.sh` only computes the missing ones. To regenerate a table (e.g. after an update of the code or of the parameter files), delete the old `TmatCoefDiff_*` and `TmatCoefInt_*` files (`.nc` or text) first.

:pencil2: NOTE : Hail and wet hail tables are not computed at Ku, K, Ka and W bands for now (`SKIP_TYPES` in `generate_tables.sh`). Dry graupel and hail tables stop at 0 °C (solid ice cannot be warmer); operadar uses the 0 °C values above. Wet graupel and wet hail tables go up to 30 °C.


## b. Customized configuration
To increase the program's flexibility, you are allowed to modify, for each hydrometeor class and frequency band combination, the parameters used to compute the tables. These files are stored in the `tables_generator/param/` folder, and named as follows : `TmatParam_{band}{hydrometeor}`. **Please do not modify** the `TmatParam_*_David2025AMT` and `TmatParam_*_David2026PhD` files.  
If you take a closer look at the `TmatParam_Crr` file, you will see that all options are filled in. The base files `TmatParam_{band}{hydrometeor}` are the default configuration (David2026PhD) at all bands, so a new configuration only needs the options you want to change. You can refer to the next section to learn more about the available parameterizations.

To use these customized files (e.g. for band C), run
```
./generate_tables.sh --newConf nameOfMyNewConf --band C
```
The base files are copied to `TmatParam_{band}{hydrometeor}_nameOfMyNewConf` (kept for traceability) and the tables are stored under `tables_generator/tables/nameOfMyNewConf/`, still named `TmatCoefInt_{microphysics}_{moments}_{band}{hydrometeor}.nc`.

# Parameterization options and references
`TmatParam_{band}{hydrometeor}` are used to launch the creation of the tables for individual hydrometeor type and frequency band. All files contain the following parameterization options that you may want to adjust for a customized configuration of the radar forward operator. Please note that some options are only available for specific types of hydrometeor.

| Key | Description | Used for | Options |
| :---: | --- | :---: | --- |
|`canting`|Turn on/off the oscillation of the particle|all hydrometeors|<ul><li>`nocant`= no canting</li><li>`orient`= canting by averaging the scattering over the orientations of the particles (Gaussian distribution of the axis around the vertical, of standard deviation `sigbeta`, uniform azimuth), valid at any elevation (recommended)</li><li>`analyt`= canting based on the analytical solutions of Ryzhkov et al. 2011, eq. 20: only valid at low elevation (a warning is printed above 10°)</li></ul>|
|`sigbeta`|Standard deviation of the canting angle.<br>If canting is turned off, this parameter will have no effect|all hydrometeors|Any float value (in degree)|
|`DIEL`|Complex dielectric function|all hydrometeors|<ul><li>Only applicable for `rr` :</li><ul><li>`Liebe91` [Liebe et al. (1991)](https://doi.org/10.1007/BF01008897)</li></ul><li>Applicable for `ss`, `gg`, `ii`, `hh` :</li><ul><li>`RY19dry` [Ryzhkov et al. (2019)](https://doi.org/10.1007/978-3-030-05093-1_4)</li></ul><li>Applicable for wet hydrometeors (`wg`, `wh`) :</li><ul><li>`LBwetgr` [Le Bastard (2019)](https://theses.hal.science/tel-04169882v1)</li></ul><li>Applicable for all hydrometeors :</li><ul><li>`MGwMA08` [Matrosov (2008)](https://doi.org/10.1109/TGRS.2008.915757) using the Maxwell Garnett weighted rule</li></ul></ul>|
|`ARfunc`| Axis ratio function |all hydrometeors|<ul><li>Only applicable for `rr`:</li> <ul><li>`BR02` [Brandes et al. (2002)](https://doi.org/10.1175/1520-0450(2002)041%3C0674:EIREWA%3E2.0.CO;2)</li></ul><li>Only applicable for `ss`:</li><ul><li>`AUds` [Augros et al. (2016)](https://doi.org/10.1002/qj.2572)</li></ul><li>Applicable for `gg` and `hh` :</li> <ul><li>`RYdg` [Ryzhkov et al. (2011)](https://doi.org/10.1175/2010JAMC2363.1)</li></ul><li>Applicable for wet hydrometeors (`wg`, `wh`) :</li><ul><li>`RYwg` [Ryzhkov et al. (2011)](https://doi.org/10.1175/2010JAMC2363.1)</li></ul><li>Applicable for all hydrometeors :</li> <ul><li>`CNST` constant axis ratio</li></ul></ul>|
|`ARcnst`| Axis ratio value |all hydrometeors|Any float value (no unit)|
|`DSTYfunc`|Density function|all hydrometeors|<ul><li>Only applicable for `ss` :<ul><li>`ZA05`= [Zawadzki et al. (2005)](https://doi.org/10.1175/JAS3563.1) --> the density formula is multiplied by `Frim` </li> <li>`BR07`= [Brandes et al. (2007)](https://doi.org/10.1175/JAM2489.1) climatological relationship</li> <li>`LS15`= [Leinonen and Szyrmer (2015)](https://doi.org/10.1002/2015EA000102) --> the selection of their density model (A or B) is based on `Frim` parameter. </li></ul></li><li>Applicable for all hydrometeors :<ul><li>`RHOX`= classic density calculation (mass/volume)</li></ul></li></ul>|
|`Frim`|Degree of riming|`ss`|This parameter is used in conjunction with the `DSTYfunc` (when applicable):<ul><li>Unrimed snow :</li><ul><li>`Frim=1` + `DSTYfunc=ZA05` (green plain curve in <a href="#rho-D_fig">Figure 1</a>) or `DSTYfunc=LS15` (red plain curve in <a href="#rho-D_fig">Figure 1</a>)</li></ul><li>Rimed snow :</li><ul><li>`Frim > 1` + `DSTYfunc=ZA05` (e.g. green dashed curve for `Frim=20` in <a href="#rho-D_fig">Figure 1</a>)</li><li>`1<Frim<3` + `DSTYfunc=LS15` : use the formula of model A (red dashed curve in <a href="#rho-D_fig">Figure 1</a>)</li><li>`Frim>3` + `DSTYfunc=LS15` : use the formula of model B (red pointed curve, not visible in <a href="#rho-D_fig">Figure 1</a>)</li></ul></ul>|
|`LAM`|Wavelength|all hydrometeors|Any float value (in mm) that matches with the `band` argument (second line in the configuration file).|

:hammer: **N.B.** The parameters `ELEV`, `Tc`, `expD`, `Fw` + `min`, `max` or `step` (e.g. `ELEVmin` `ELEVmax` `ELEVstep`) are used to compute the scattering coefficients over a range of elevations, temperatures, diameters, and liquid water fraction. For instance, `ELEVmin=0` `ELEVmax=90` `ELEVstep=90` will produce tables with only `0` or `90°` elevation angles (simulate horizontal or vertical pointing radar). 

The steps of the predefined configurations were chosen with a grid accuracy study, so that the linear interpolation of operadar stays within 0.05 dB for Zh, 0.02 dB for Zdr and 1 % for Kdp and Ah (elevation step 10°; diameter step `expDstep` 0.02; temperature step 5 °C for ice, snow, rain and cloud (2 °C for Ka-band rain), 10 °C for graupel and hail, 2 °C for wet graupel/hail (1 °C at K, Ka, W); Fw step 0.01). See `tables_generator/tools/grid_study/SUMMARY.md` before changing them. The integrated tables use a content step of 0.02 in log10(M) and, for 2-moment species, a concentration step of 0.1 in log10(N).

***

<a name="rho-D_fig"> <img src="https://github.com/user-attachments/assets/938d2079-e345-4049-842b-9b7fffd81e94"> </a>
<p align=center>
<strong>Figure 1 : Different density-diameter relationships implemented in OPERADAR. Diameters (x-axis) are in mm.</strong>
<br>This figure has been realized with the online plot tool Desmos.</p>

***

# 3. Checking the tables with the plot tools
The scripts of the `plot_tools/` folder plot the tables (run them from this folder; figures are written to `plot_tools/IMG/`). The full guide, with all the options, is [plot_tools/README.md](https://github.com/UMR-CNRM/operadar/blob/master/plot_tools/README.md).

| Script | What it shows | Configuration |
|---|---|---|
| `plot_tables_2D.py` | 2D maps of the full tables: content × temperature, content × Fw, content × N, diameter × temperature, diameter × Fw | command line |
| `plot_D_M_varpol_alltypes.py` | curves of all types on one figure (one band), as a function of the content M or of the diameter D | variables at the top of the file |
| `plot_D_M_varpol_Rayleigh.py` | one panel per type, one curve per band, T-matrix and Rayleigh | variables at the top of the file |
| `plot_theoretical_varpol.py` | sensitivity to the table options (axis ratio, density, dielectric function, Fw, N...) | command line |
| `plot_2Dmaps.py` | maps of the simulated radar variables (operadar outputs) | variables at the top of the file |

**Start with the 2D maps** to check a new set of tables:
```
cd plot_tools
python plot_tables_2D.py --band C                                   # ICE3, all types, David2026PhD tables
python plot_tables_2D.py --band C --micro LIMA --types rr cl cs ii  # 2-moment LIMA species
python plot_tables_2D.py --band W --config David2025AMT --types ss --vars Zh Zdr --Tc -20
python plot_tables_2D.py --band C --types rr wg --elev 60             # at 60° elevation
```
For each type they show Zh, Zdr, Kdp, rhohv, Ah and the difference with the Rayleigh approximation (Mie effects). The maps of the integrated tables also show `Mfrac`, a variable of the `TmatCoefInt_*` tables: the fraction of the content M carried by particles inside the diameter range of the table, computed from the size distribution. **Hatched areas (Mfrac < 0.99) are where the size distribution is cut by the table, so the values there are too low.** In realistic conditions this only happens for ICE3 snow above about 6 g/m³, for very heavy ice crystals and for drizzle-size LIMA cloud drops.

Then use the curves for a closer look, e.g. in `plot_D_M_varpol_alltypes.py` set
```python
micro="ICE3"
band,color_list,typeh_list ='C',['b','k','c','g','gray','gray','r'],['rr','cl','ss','ii','gg','wg','hh']
listplot=['M','D']
listvar=['Zh','Zdr','Kdp']
```
and run `python plot_D_M_varpol_alltypes.py`. When a requested temperature, Fw or concentration is not on the table grid, the scripts use the closest value and print a warning.
