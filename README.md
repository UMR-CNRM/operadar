# OPERADAR "Observation oPErator for polarimetric RADAR"
OPERADAR is developed at CNRM (National Centre for Meteorological Research) 
to compute synthetic dual-polarization radar variables from model fields in a 3D grid.

The computed variables include reflectivity $Z_{H}$ , differential reflectivity $Z_{DR}$ , specific differential phase shift $K_{DP}$ , cross-correlation coefficient $\rho_{HV}$ and specific attenuation at horizontal and vertical polarizations $A_{H}$ and $A_{v}$ (3D).

It is adapted to read AROME and MesoNH atmospheric model files, and the available options (mass-density relations, particle size distributions) are consistent with the microphysics options currently available for these models.

OPERADAR is used for research applications, mainly for the evaluation of AROME and MesoNH microphysics schemes in the polarimetric radar space.
It currently requires a repository of lookup tables. These tables contain the scattering coefficients to compute the polarimetric variables, either in the Rayleigh or Mie scattering regime (with the T-matrix method). 

<span style='color:DodgerBlue;font-weight:bold'>------- NEW SINCE JUNE 2025 ! -------</span>
<br>OPERADAR software now includes an executable to generate the lookup tables. The scattering coefficients stored in the tables are computed both in the Rayleigh approximation and using the T-matrix method (see [Mishchenko and Travis 1994](https://www.sciencedirect.com/science/article/pii/0030401894907315?via%3Dihub)).

### INPUT
* Model file (`.fa` for AROME files or `.nc` for MesoNH)
* Configuration file
* Lookup tables (created beforehand with the table generator executable).
* For colleagues at CNRM, the default configuration of the tables can be found on HPC at `/home/augros/TmatCoefInt_SCXKaW/`

### OUTPUT : netcdf file with
* X/Y grid horizontal coordinates (1D)
* model pressure levels vertical coordinates (1D)
* latitude/longitude fields (2D)
* dual-pol variables : $Z_{H}$ , $Z_{DR}$ , $K_{DP}$ , $\rho_{HV}$ (3D)
* specific attenuation : $A_{H}$ , $A_{v}$ (3D)
* temperature, altitude, contents and concentrations fields (3D)
* dataset attributes such as the horizontal resolution, the options used to compute the mixed-phase or simulate the radar data.

# How to install
The installation of `operadar` should be preferred in a virtual environment. The procedure is detailled [in the Wiki](https://github.com/UMR-CNRM/operadar/wiki/Installation-tutorial).


# How to run
1) To work properly, lookup tables should be generated **before** running `operadar` (unless you have a direct access to the tables).
Steps to create the tables are available in [Generation of the lookup tables](https://github.com/UMR-CNRM/operadar/wiki/Generation-of-the-lookup-tables).<br>*--> If you already have access to the tables, you may skip this step.*
2) A configuration file is red by `operadar` to compute the polarimetric fields, according to your chosen parameterization. See [Simulation configuration](https://github.com/UMR-CNRM/operadar/wiki/Simulation-configuration) for extended explanations.
3) Finally, a shell executable is provided to launch the forward operator for one file at a time. However, `operadar` can be used within other algorithm, to loop on multiple files and eventually with varying configurations over files that share certain common parameters. Explanations are in [Execution of the forward operator](https://github.com/UMR-CNRM/operadar/wiki/Execution-of-the-forward-operator)

# Contributing
If you wish to contribute to the project, first, fork the code to create your own copy of the project (see https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project#creating-your-own-copy-of-a-project).

Before making changes to the project, you should create a new branch and check it out. Try to be self explanatory for the name of the branch. Then, install operadar as a package with the aforementioned procedure.
To ask for integration of your modifications, please open a pull request following https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project#making-a-pull-request.

# License
This software is governed by the open-source CeCILL-C license under French law, cf. LICENSE.txt. Downloading and using this code means that you have had knowledge of the CeCILL-C license and that you accept its terms.

