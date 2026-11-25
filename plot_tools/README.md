# Plot tools

Scripts to look at the lookup tables and at the operadar outputs. Run them from this folder (`cd plot_tools`);
figures are written to `plot_tools/IMG/` (not tracked by git). They need `numpy`, `pandas`, `matplotlib`, `xarray`
and `scipy` (all installed with operadar); only `plot_2Dmaps.py` needs `cartopy` and `epygram` (`cartopy` is
optional for `plot_max_over_time.py`). Every command-line script has a help: `python <script> -h`.

**Lookup tables**

| Script | What it shows | How it is configured |
|---|---|---|
| [`plot_tables_2D.py`](#plot_tables_2dpy-2d-maps-of-the-full-tables) | 2D maps of the full tables (all temperatures, contents, Fw, N) | command line |
| [`plot_D_M_varpol_alltypes.py`](#plot_d_m_varpol_alltypespy-all-types-on-one-figure) | curves of all types on one figure, one band | variables at the top of the file |
| [`plot_D_M_varpol_Rayleigh.py`](#plot_d_m_varpol_rayleighpy-several-bands-t-matrix-and-rayleigh) | one panel per type, several bands, T-matrix vs Rayleigh | variables at the top of the file |
| [`plot_theoretical_varpol.py`](#plot_theoretical_varpolpy-sensitivity-to-the-table-options) | sensitivity to the table options (axis ratio, density, Fw...) | command line |

**Operadar outputs** (see [Analysis of the operadar outputs](#analysis-of-the-operadar-outputs))

| Script | What it shows | How it is configured |
|---|---|---|
| [`compute_cfad.py`](#cfads-as-in-david-et-al-2025-amt) + [`plot_cfad.py`](#cfads-as-in-david-et-al-2025-amt) | CFADs of Zh, Zdr, Kdp, contents and Fw inside storm cells, simulations vs observations (David et al. 2025, AMT, Fig. 6 and 8) | command line |
| [`plot_output_check.py`](#plot_output_checkpy-quick-look-and-sanity-checks) | quick look and sanity checks of one output: maps, joint distributions Zh-Zdr/Kdp, variables vs temperature, warnings | command line |
| [`plot_cross_section.py`](#plot_cross_sectionpy-vertical-cross-sections) | vertical cross-sections of the radar variables and contents, with isotherms | command line |
| [`plot_compare_outputs.py`](#plot_compare_outputspy-compare-two-outputs) | differences between two outputs of the same model file (tables, configurations, bands) | command line |
| [`plot_max_over_time.py`](#plot_max_over_timepy-maximum-over-several-times) | maximum of a variable over several times (swaths) with the cell contours | command line |
| [`plot_2Dmaps.py`](#plot_2dmapspy-maps-of-operadar-outputs) | maps of the simulated variables at model levels (needs cartopy and epygram) | variables at the top of the file |
| `draft/` | former scripts, kept for reference (see `draft/README.md`) | |

Tables are read from `../tables_generator/tables/<config>/`, in NetCDF (`<table>.nc`; text tables of older
versions are still read):
- `TmatCoefDiff_<band><type>`: single-particle values as a function of the diameter D (axis "D" in the scripts);
- `TmatCoefInt_<micro>_<n>M_<band><type>`: values integrated over the size distribution, as a function of the
  content M (axis "M"). The third parameter `P3` is Fw for wet species and the concentration N for 2-moment
  species and pristine ice; `Mfrac` is the fraction of M inside the diameter range of the table.

The scripts read the tables with `utils/tables.py` (which uses `operadar/read/table_io.py`) as a table with the
same columns as the former text files. The NetCDF tables can also be plotted directly with xarray, e.g.
`xr.open_dataset("TmatCoefInt_ICE3_1M_Crr.nc").zdr.sel(Tc=10, ELEV=0, Fw=0).plot(x="M", xscale="log")`.

Types: `rr` rain, `cl`/`cs` cloud water (land/sea), `ii` pristine ice, `ss` snow, `gg` graupel, `hh` hail,
`wg`/`wh` wet graupel/hail. Bands: L, S, C, X, Ku, K, Ka, W (hail and wet hail are not computed at Ku, K, Ka, W).

When a requested temperature, Fw or concentration is not on the table grid, the scripts use the closest grid
value and print a warning such as `/!\ P3=800 is not in the table, using the closest value P3=794.0`.

---

## `plot_tables_2D.py`: 2D maps of the full tables

The best way to check a whole table at once. For each type it draws Zh, Zdr, Kdp, rhohv, Ah and
Zh - Zh Rayleigh (Mie effects) as colour maps:

| Figure | x | y | fixed |
|---|---|---|---|
| `Int_..._M-Tc.png` | content M | temperature | Fw (wet species) or N (2 moments, ice) |
| `Int_..._M-Fw.png` (wet species) | M | Fw | temperature |
| `Int_..._M-N.png` (2 moments, ice) | M | N | temperature |
| `Diff_..._D-Tc.png` | diameter D | temperature | Fw |
| `Diff_..._D-Fw.png` (wet species) | D | Fw | temperature |

The integrated maps also show `Mfrac`. **Hatched areas (Mfrac < 0.99) are where the size distribution is cut
by the table**: values there are too low (or meaningless) and are left out of the colour scales. Fields that do
not vary (e.g. rhohv = 1 and Zdr = 0 for spheres) are written as a constant instead of a colour map.

```bash
python plot_tables_2D.py --band C                                   # ICE3, all types, David2026PhD tables
python plot_tables_2D.py --band C --micro LIMA --types rr cl cs ii  # 2-moment species of LIMA
python plot_tables_2D.py --band W --config David2025AMT --types ss --vars Zh Zdr --Tc -20
python plot_tables_2D.py --band C --types wg --Fw 0.5 --elev 30     # wet graupel at 30 deg elevation
```

Options: `--band`, `--config` (table folder: `David2026PhD` by default, L, S, C; `David2025AMT` for all bands), `--micro` (ICE3, LIMA, ICJW, LIMC), `--types`,
`--vars` (Zh Zdr Kdp Rhohv Ah Zh-ZhR), `--elev` (elevation, 0 to 90 every 10 deg; default 0), `--Tc` (temperature of the
M-Fw, M-N and D-Fw maps), `--Fw` (Fw of the M-Tc maps of wet species, default 0.2), `--N` (N of the M-Tc maps of
2-moment species and ice), `--out`. Figures go to `IMG/tables2D/<config>/<band>/`. Missing tables are skipped.
The C band set takes about a minute.

---

## `plot_D_M_varpol_alltypes.py`: all types on one figure

Curves of Zh (and/or Zdr, Kdp, Rhohv, Ah, Av) as a function of the content M or of the diameter D, for all
types of one band, at one temperature per type. Edit the variables at the top of the file, then run
`python plot_D_M_varpol_alltypes.py`.

| Variable | Meaning |
|---|---|
| `micro` | microphysics scheme of the M tables (ICE3, LIMA...) |
| `band, color_list, typeh_list` | band, one colour per type, types to plot (e.g. add `'hh'`) |
| `listplot` | `['M']`, `['D']` or `['M','D']` (one figure each) |
| `listvar` | variables, one panel each (e.g. `['Zh','Zdr','Kdp']`) |
| `TmatOption` | table folder per band (`David2026PhD` for L, S, C; `David2025AMT` for all bands) |
| `T_dict` | temperature of each type (hail: 0 °C at most) |
| `Fw_list` | Fw of the wet graupel/hail curves |
| `Nii`, `N2M` | concentration of pristine ice and of the other 2-moment species (LIMA rain, cloud) |
| `ELEVsel` | elevation in degrees (0 to 90, every 10): 0 for low ground scans, 90 for vertically pointing radars |

Output: `IMG/DistTmat<band><micro><config>T...varpol_alltypes<config>.png`.

## `plot_D_M_varpol_Rayleigh.py`: several bands, T-matrix and Rayleigh

One panel per type, one colour per band, with the Rayleigh approximation as dashed lines (`plotR = True`). Useful
to see where Mie effects start for each band. Edit the top of the file (`micro`, `TmatOption`, `band_list`,
`typeh_list`, `plot_list`, `var_list`, `ELEVsel`, `T_dict`, `Fw_list`, `expN_list` for 2-moment species) and run
`python plot_D_M_varpol_Rayleigh.py`. The tables of every band of `band_list` must exist in the folder given by
`TmatOption` (e.g. `'W':'David2025AMT'` needs `tables/David2025AMT/`).

## `plot_theoretical_varpol.py`: sensitivity to the table options

Compares tables computed with different options (axis ratio, canting, density function, riming, dielectric
function, Fw, concentration). It reads one folder per set of options:
```
../tables_generator/tables/<type>/<ARfunc>_AR<ARvalue>_CANT<sigbeta>_<DSTYfunc>_Frim<Frim>_<DIEL>/
    e.g. tables/rr/BR02_AR1.00_CANT0.0_RHOX_Frim1.0_Liebe91/TmatCoefInt_ICE3_1M_Crr
```
(create each set with `generate_tables.sh --newConf <name>` and move or link the tables there).
The reference options are `base_cfg` in `utils/utils_plot.py`. Run from this folder:
```bash
PYTHONPATH=.. python plot_theoretical_varpol.py rr C --axe M
PYTHONPATH=.. python plot_theoretical_varpol.py wg C --Fw 0.1 0.6 0.9
PYTHONPATH=.. python plot_theoretical_varpol.py rr C --axe M --micro LIMA --Nc 100 1000 10000
PYTHONPATH=.. python plot_theoretical_varpol.py ss C --DSTYfunc RHOX BR07 --ref
```
See `python plot_theoretical_varpol.py -h` for all options (`--ARvalue`, `--ARfunc`, `--CANTING`, `--DSTYfunc`,
`--FRIM`, `--DIEL`, `--combine`, `--invert`, `--method Tmatrix|Rayleigh|Both`). Figures go to
`IMG/<type>_<axis>xaxis_<band>band/`.

## `plot_2Dmaps.py`: maps of operadar outputs

Maps of the simulated variables (Zh, Zdr, Kdp, Ah, Rhohv, Zh_att) from operadar NetCDF files of real cases
(the configuration must have `real_case=True`), at chosen model levels and as the column maximum. Edit
`dataDir`, `outputDir`, `variables`, `levels` and the domain at the top of the file, then run
`python plot_2Dmaps.py`. Needs `xarray`, `cartopy` and `epygram`.

---

# Analysis of the operadar outputs

These scripts read the operadar NetCDF outputs directly (`dpolvar_<model>_<micro>_<band>band_<date>.nc`, with
`Zh`, `Zdr`, `Kdp`..., `Contents`, `T` and `Alt` on the model levels), for AROME and Meso-NH, real or idealised
cases. Model levels are interpolated on regular heights where needed. They also read **gridded observations**:
any NetCDF file with a vertical dimension `z` (or `height`, `alt`) in metres, possibly a `time` dimension, and
the variables `zh`, `zdr`, `kdp`, `rhohv` (lower or upper case), e.g. radar data interpolated on a Cartesian
grid with Py-ART. The common code is in `utils/operadar_output.py` and `utils/cfad.py`.

## CFADs as in David et al. (2025, AMT)

Two steps: `compute_cfad.py` reads all the files of one dataset (simulation or observations) and saves the
distribution of each variable at each height in a small statistics file; `plot_cfad.py` plots several datasets
together. Computing is the long step (about 1 s per AROME file of 130 x 120 points, 25 s and 3.5 GB of memory
for a Meso-NH file of 1200 x 1200 points); the plots are then instantaneous and can be changed freely.

```bash
# 1) statistics inside the storm cores, for each dataset
python compute_cfad.py "/path/AROME_ICE3/dpolvar_*.nc" --label "AROME ICE3" --out CFAD/ICE3_core.nc
python compute_cfad.py "/path/AROME_LIMA/dpolvar_*.nc" --label "AROME LIMA" --out CFAD/LIMA_core.nc
python compute_cfad.py "/path/OBS/*.nc" --label obs --core-dbz 48 --out CFAD/obs_core.nc
# 2) figures (IMG/cfad/)
python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc                    # Fig. 6 of the paper
python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc --figs freq --diff # frequency CFADs
```

**Where the values are taken** (`compute_cfad.py --mask`), from the column maximum of Zh, as in the paper
(Sect. 3.3): `core` (default) = contiguous areas with Zmax >= 40 dBZ (`--core-dbz`) and at least 5 grid points
(`--min-core`), inside an envelope; `envelope` = contiguous areas with Zmax >= 25 dBZ (`--envelope-dbz`) and at
least 20 points (`--min-envelope`); `zh` = every point with Zh >= `--zh-local` (3D criterion); `all`; `var` =
where a variable of the files is > 0 (`--mask-var`, e.g. a cell or ZDR column mask from a tracking algorithm).
With `--cells-from FILE...` (one per input file, same grid) the cores/envelopes are those of other files, e.g. the
C band output for a W band output of the same run, whose Zh rarely reaches 40 dBZ.
The paper uses **48 dBZ for the observations and 40 dBZ for the simulations**. The radar variables are only
used where Zh >= `--zh-min` (0 dBZ); the contents everywhere in the mask (zeros included, as in the paper); the
liquid water fraction Fw of the wet graupel/hail where there is some (it is rebuilt from the saved contents,
which is only possible with `mixed_phase_parametrization = 'Fw_posg'`). Heights go from 0 to 15 km every
500 m (`--hmin`, `--hmax`, `--dz`), above the lowest model level (`--height agl`, close to the height above
ground, as in the paper) or above sea level (`--height asl`). `--box LONMIN LONMAX LATMIN LATMAX` keeps a
sub-domain.

**Figures** (`plot_cfad.py --figs`): `summary` (default, layout of Fig. 6: radar variables of all datasets,
contents of each simulation, Fw), `radar`, `contents`, `fw`, `freq` (2D histograms, one figure per variable,
`--diff` adds the difference with the first dataset), `all`. Options: `--vars Zh Zdr Kdp Rhohv`,
`--extremes` (Q5 and Q95), `--xlim VAR MIN MAX`, `--hmax 12`, `--min-count 20` (heights with fewer values are
hidden), `--contents-log`, `--freq-norm level|total|count`, `--bin VAR WIDTH`, `--log`, `--labels`, `--colors`
(observations are grey by default, then orange, green...). **Statistics files with the same label are added**,
so cases can be computed separately (or in parallel) and plotted together.

How the CFADs are computed: the statistics files hold, for each variable and each height, the histogram of the
values with fine bins (0.1 dBZ, 0.002 dB, 0.002 deg/km, 3 % for the contents; exact zeros counted apart). The
quantiles are taken from these histograms: they agree with the exact percentiles within one fine bin (checked
on an AROME case: < 0.05 dBZ, < 0.01 dB or deg/km, < 1 % for the contents).

## `plot_output_check.py`: quick look and sanity checks

To run after each new simulation. Prints a summary of each variable (valid points, NaN fraction, percentiles,
contents) with warnings (Zh < -30 dBZ: radar sensitivity not applied; points with precipitation but no Zh:
values outside the tables?; Zdr < -1 or > 8 dB, Kdp < -1 deg/km, Rhohv > 1 or < 0.7), and draws the column
maximum of each radar variable and content, the joint distributions Zh-Zdr, Zh-Kdp (Zh-Rhohv) in rain and in ice,
and the median and Q5-Q95 of each variable as a function of temperature.
```bash
python plot_output_check.py ../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
python plot_output_check.py "../operadarFiles/*.nc" --figs maps joint --t-rain 4 --t-ice -10
```
Figures and the text summary go to `IMG/check/`. For large files, percentiles and joint distributions use a
random subsample of 5 million points (`--max-points`).

## `plot_cross_section.py`: vertical cross-sections

Radar variables (and with `--contents` the hydrometeor contents) along a line, with isotherms, and a map of the
column maximum of Zh showing the line. Default: west-east line through the maximum of Zh.
```bash
python plot_cross_section.py FILE                                     # through the strongest cell
python plot_cross_section.py FILE --points 2.3 45.1 2.8 45.9 --contents --hmax 13
python plot_cross_section.py FILE --y 60 --vars Zh Zdr Kdp Rhohv --isotherms 0 -10 -20 -40
```
Other positions: `--x KM`, `--ij I1 J1 I2 J2` (grid indices). Figures go to `IMG/sections/`.

## `plot_compare_outputs.py`: compare two outputs

For two outputs computed from the same model file (old and new tables, two configurations, two bands,
T-matrix and Rayleigh): statistics of B - A by temperature range (rain, melting layer, ice), scatter plots,
maps of the column maximum of A, B and B - A, and B - A as a function of temperature.
```bash
python plot_compare_outputs.py old/dpolvar_X.nc new/dpolvar_X.nc --labels "old tables" "new tables"
```
Figures and statistics go to `IMG/compare/`.
If one output was computed on a subdomain (`subDomain`) and the other on the whole domain, the larger one
is cut to the subdomain automatically; `--box LONMIN LONMAX LATMIN LATMAX` compares only a lon/lat box.

## `plot_max_over_time.py`: maximum over several times

Maximum of a variable over all the files (one per time) or all the times of a gridded file, at a height
(`--height 2`, km) or over the column (default), with the contour Zh = 40 dBZ of each time (`--contour`): storm
swaths and cell tracks, e.g. to compare the paths of simulated and observed cells.
```bash
python plot_max_over_time.py "../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20220818_0*.nc" --var Zh
python plot_max_over_time.py OBS/20220818_obs.nc --var Zdr --height 2 --box 4 11 41 45
```
Figures go to `IMG/maxtime/`.
