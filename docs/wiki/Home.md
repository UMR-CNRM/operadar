# Description of the software
OPERADAR "Observation oPErator for polarimetric RADAR" is developed at CNRM (National Centre for Meteorological Research) 
to compute synthetic dual-polarization radar variables from model fields in a 3D grid.

The computed variables include reflectivity $Z_{H}$ , differential reflectivity $Z_{DR}$ , specific differential phase shift $K_{DP}$ , cross-correlation coefficient $\rho_{HV}$ and specific attenuation at horizontal and vertical polarizations $A_{H}$ and $A_{v}$ (3D).

It is adapted to read AROME and MesoNH atmospheric model files, and the available options (mass-density relations, particle size distributions) are consistent with the microphysics options currently available for these models.

OPERADAR is used for research applications, mainly for the evaluation of AROME and MesoNH microphysics schemes in the polarimetric radar space.

# OPERADAR in a glimpse
###  INPUTs
To work properly `operadar` requires the following elements or files : 
* Model file (`.fa` for AROME files or `.nc` for MesoNH)
* Configuration file
* Lookup tables (created beforehand with the table generator executable). 

### OUTPUTs
The output fields are stored in a `.nc` file. An `--append` option is available to write the dual-pol fields into the input model file (see [Execution of the forward operator](https://github.com/UMR-CNRM/operadar/wiki/Execution-of-operadar) for more).
* X/Y grid horizontal coordinates (1D)
* model pressure levels vertical coordinates (1D)
* latitude/longitude fields (2D)
* dual-pol variables : $Z_{H}$, $Z_{DR}$, $K_{DP}$, $\rho_{HV}$ (3D)
* specific attenuation : $A_{H}$, $A_{V}$ (3D)
* temperature, altitude, contents and concentrations fields (3D)
* dataset attributes such as the horizontal resolution, the options used to compute the mixed-phase or simulate the radar data.

To check and analyse the outputs (CFADs as in David et al. 2025, quick look and sanity checks, cross-sections,
comparison of two runs, storm swaths), see [Plot tools](https://github.com/UMR-CNRM/operadar/wiki/Plot-tools).

### How to get the tables
* you can generate the tables with `generate_tables.sh` (see [Generation of the lookup tables](https://github.com/UMR-CNRM/operadar/wiki/Generation-of-the-lookup-tables) and the predefined [Lookup table configurations](https://github.com/UMR-CNRM/operadar/wiki/Lookup-table-configurations))
* at Météo-France you can access the default tables on HPC or Lustre (contact clotilde.augros@meteo.fr)

### Available radar bands :satellite: 
| Name | Wavelength | Frequency |
|:----:|:----------:|:---------:|
|`W`| 3.16 mm | 95 GHz | 
|`Ka`| 8.40 mm | 35,75 GHz |
|`K`| 15.0 mm | 20 GHz |
|`Ku`| 20.0 mm | 15 GHz |
|`X`| 31.9 mm | 9,4 GHz | 
|`C`| 53.2 mm | 5,638 GHz |
|`S`| 106.2 mm | 2,824 GHz |
|`L`| 190.4 mm | 1.57542 GHz |

### Available hydrometeor types :snowflake:
| Shortcut | Corresponding type |
|:--------:|:------------------:|
|`cl`| cloud water over land |
|`cs`| cloud water over sea |
|`rr`| rain |
|`gg`| graupel |
|`hh`| hail |
|`ss`| snow |
|`ii`| pristine ice |
|`wh`| wet hail |
|`wg`| wet graupel |

# How to cite

If you use `operadar` in a published work, please cite the following papers:

* Augros, C., Caumont, O., Ducrocq, V., Gaussiat, N., & Tabary, P. (2016). Comparisons between S, C, and X band polarimetric radar observations and convective-scale simulations of HyMeX first special observing period. Quarterly Journal of the Royal Meteorological Society, 142, Issue S1, 347–362. [https://doi.org/10.1002/qj.2572](https://doi.org/10.1002/qj.2572)

* David, Cloé, et al. “Improved Simulation of Thunderstorm Characteristics and Polarimetric Signatures with LIMA Two-Moment Microphysics in AROME.” Atmospheric Measurement Techniques, vol. 18, no. 15, Aug. 2025, pp. 3715–45. Copernicus Online Journals, https://doi.org/10.5194/amt-18-3715-2025.


You should also cite the Fortran 77 T-matrix code that is employed to generate the scattering coefficients with the T-matrix method:

* Mishchenko, M. I., & Travis, L. D. (1994). T-matrix computations of light scattering by large spheroidal particles. Optics Communications, 109(1), 16–21. [http://www.sciencedirect.com/science/article/pii/0030401894907315](http://www.sciencedirect.com/science/article/pii/0030401894907315)



# Major updates

## April 2026
first release v1.0.0
include attenuation (gaz + hydrometeors) integration along the beam (for vertical pointing radar)
improve code robustness (addition of test functions, improvement of the hypercube interpolation function for LUT) 

## June 2025
OPERADAR software now includes an executable to generate the lookup tables. The scattering coefficients stored in the tables are computed both in the Rayleigh approximation and using the T-matrix method (see [Mishchenko and Travis 1994](https://www.sciencedirect.com/science/article/pii/0030401894907315?via%3Dihub)).

## April 2025
OPERADAR software now handles Meso-NH model files
