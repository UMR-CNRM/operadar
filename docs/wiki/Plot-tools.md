The folder `plot_tools/` contains scripts to check and analyse the **operadar outputs** (CFADs as in David et
al. 2025, quick look, cross-sections, comparison of two runs, storm swaths) and to look at the **lookup
tables**. This page explains how to use them; all options are listed in
[plot_tools/README.md](https://github.com/UMR-CNRM/operadar/blob/master/plot_tools/README.md) and by
`python <script> -h`.

# 1. Before starting

* The scripts only need the Python packages installed with operadar (`numpy`, `pandas`, `matplotlib`, `xarray`,
  `scipy`). `cartopy` is optional (coastlines in `plot_max_over_time.py`); only the older `plot_2Dmaps.py` also
  needs `epygram`.
* Run them from the folder `plot_tools/` (`cd plot_tools`). Figures go to `plot_tools/IMG/<kind>/` (not tracked
  by git), or to the folder given with `--out`.
* Input files: the NetCDF files written by operadar (`dpolvar_<model>_<micro>_<band>band_<date>.nc`), AROME or
  Meso-NH, real or idealised cases. File names accept shell patterns between quotes, e.g.
  `"../operadarFiles/AROME/dpolvar_*.nc"`.
* To compare with **observations**, the radar data must be on a 3D Cartesian grid (e.g. interpolated with
  Py-ART): one NetCDF file with a vertical dimension `z` (or `height`, `alt`) in metres, optionally a `time`
  dimension, and the variables `zh`, `zdr`, `kdp`, `rhohv` (lower or upper case). lat/lon (2D) are needed only for
  maps in lon/lat and `--box`.

| Script | Use |
|---|---|
| `compute_cfad.py` + `plot_cfad.py` | CFADs of Zh, Zdr, Kdp, contents and Fw inside storm cells, simulations vs observations |
| `plot_output_check.py` | quick look and sanity checks of a new output |
| `plot_cross_section.py` | vertical cross-sections with isotherms |
| `plot_compare_outputs.py` | differences between two outputs of the same model file |
| `plot_max_over_time.py` | maximum over several times (swaths) with the cell contours |
| `plot_tables_2D.py`, `plot_D_M_varpol_*.py`, `plot_theoretical_varpol.py` | lookup tables (see section 4) |

# 2. CFADs as in David et al. (2025, AMT)

CFADs (Contoured Frequency by Altitude Diagrams, Yuter and Houze 1995) show the distribution of a variable at
each height. Fig. 6 and 8 of [David et al. (2025)](https://doi.org/10.5194/amt-18-3715-2025) show, inside the
storm cores of all cases and time steps, the median (line) and the Q25-Q75 interval (shaded) of Zh, Zdr and Kdp
(observations and simulations), of the hydrometeor contents and of the liquid water fraction of the wet
graupel. The same figures are made in two steps.

### Step 1: statistics of each dataset (`compute_cfad.py`)

One command per dataset (one simulation, or the observations), over all its files:
```bash
python compute_cfad.py "/path/AROME_ICE3/dpolvar_*.nc" --label "AROME ICE3" --out CFAD/ICE3_core.nc
python compute_cfad.py "/path/AROME_LIMA/dpolvar_*.nc" --label "AROME LIMA" --out CFAD/LIMA_core.nc
python compute_cfad.py "/path/OBS/*.nc" --label obs --core-dbz 48 --out CFAD/obs_core.nc
```
The output is a small NetCDF file with, for each variable and each height, the histogram of the values (fine
bins). Time: about 1 s per AROME file of 130 x 120 points; 25 s and 3.5 GB of memory for a Meso-NH file of
1200 x 1200 x 82 points.

**Storm cells** are found as in the paper (Sect. 3.3) from the column maximum of Zh, at each time step:

| `--mask` | Values taken in | Options (default) |
|---|---|---|
| `core` (default) | storm cores: contiguous areas with Zmax >= 40 dBZ, inside an envelope | `--core-dbz 40`, `--min-core 5` grid points |
| `envelope` | storm cells: contiguous areas with Zmax >= 25 dBZ | `--envelope-dbz 25`, `--min-envelope 20` |
| `zh` | every point with Zh >= threshold (3D criterion) | `--zh-local 40` |
| `all` | every point | |
| `var` | where a variable of the files is > 0 (2D or 3D), e.g. a cell or ZDR column mask from a tracking | `--mask-var NAME` |

:warning: The paper uses **48 dBZ for the observations** and **40 dBZ for the simulations** (simulated
reflectivities are lower): give `--core-dbz 48` for the observations.

Other options:
* `--zh-min 0`: the radar variables are only used where Zh >= 0 dBZ; the contents are used everywhere in the mask
  (zeros included, as in the paper).
* `--hmin 0 --hmax 15 --dz 0.5` (km): heights of the CFADs. Model levels are interpolated linearly at these
  heights, **above the lowest model level** (`--height agl`, close to the height above ground used in the paper)
  or above sea level (`--height asl`). Observations are taken at the nearest height of their grid.
* `--vars Zh Zdr Kdp Rhohv`: radar variables (default: Zh, Zdr, Kdp, Rhohv found in the files);
  `--zh-var Zh_att` to find the cells with the attenuated reflectivity.
* `--box LONMIN LONMAX LATMIN LATMAX`: sub-domain.
* `--cells-from FILE...`: take the cores/envelopes from other files on the same grid (one per input file), e.g.
  the C band output for the W band output of the same run (at W band Zh rarely reaches 40 dBZ).
* The CFAD of **Fw** (liquid water fraction of the wet graupel and wet hail) is rebuilt from the saved contents,
  Fw = rr / (rr + wg + wh). This is only possible with `mixed_phase_parametrization = 'Fw_posg'` (with
  `'Fw_pos'` the rain of the melting layer is moved into the wet graupel); otherwise it is skipped with a message.

### Step 2: figures (`plot_cfad.py`)

```bash
python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc          # layout of Fig. 6
python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc --figs radar --extremes
python plot_cfad.py CFAD/obs_core.nc CFAD/ICE3_core.nc CFAD/LIMA_core.nc --figs freq --diff
```

| `--figs` | Figure |
|---|---|
| `summary` (default) | Fig. 6 layout: (a, b, c) radar variables of all datasets, (d, e...) contents of each simulation, last panel Fw |
| `radar` | radar variables only; `--extremes` adds Q5 (dashed) and Q95 (dotted) |
| `contents` | hydrometeor contents, one panel per simulation (`--contents-log` for a log axis) |
| `fw` | Fw of the wet graupel: median, Q25-Q75 and Q95 (dotted) |
| `freq` | frequency CFADs (2D histograms), one figure per variable; `--diff`: difference with the first dataset (put the observations first); `--log`; `--freq-norm level` (default, % of the values of each height), `total` or `count`; `--bin Zh 1` |
| `all` | all of them |

Useful options: `--vars`, `--xlim Zh 0 65` (repeatable, `M` for the contents), `--hmax 12`, `--labels`,
`--colors` (observations are grey by default, then orange, green, blue...), `--min-count 20` (heights with fewer
values are hidden), `--name`, `--out`.

**Several cases**: run `compute_cfad.py` once per case (possibly in parallel) with the **same `--label`**:
`plot_cfad.py` adds the statistics files that have the same label.

**Accuracy**: the quantiles are computed from histograms with fine bins (0.1 dBZ, 0.002 dB, 0.002 deg/km, 3 % for
the contents, exact zeros counted apart). On an AROME case they agree with the exact percentiles within
0.05 dBZ, 0.01 dB or deg/km and 1 % for the contents.

# 3. Other figures to check and analyse the outputs

### Quick look and sanity checks: `plot_output_check.py`
To run after each new simulation:
```bash
python plot_output_check.py ../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
```
* A text summary (also saved in `IMG/check/<file>_check.txt`): valid points, NaN fraction and percentiles of each
  radar variable, maximum and frequency of each content, Fw, and **warnings**: Zh below the default
  sensitivity of the band (-30 dBZ L to X, -40 Ku to Ka, -50 W: radar sensitivity not applied, see
  `sensitivity_threshold` in [Simulation configuration](https://github.com/UMR-CNRM/operadar/wiki/Simulation-configuration)),
  points at exactly Zh = 0 dBZ (false echoes of operadar versions before Nov. 2026), points with precipitation
  (contents >= 0.1 g/m3) but no Zh (outside the tables? those with Ah were removed by the sensitivity after
  attenuation and are only reported), Zdr < -1 or > 8 dB, Kdp < -1 deg/km (expected at Ka and W band: only
  reported), Rhohv > 1 or < 0.7.
* Maps of the column maximum of the radar variables and of the contents.
* Joint distributions Zh-Zdr and Zh-Kdp (Zh-Rhohv) in rain (T >= 4 °C) and in ice (T <= -10 °C), with the median:
  in rain, Zdr and Kdp must increase with Zh.
* Median and Q5-Q95 of each variable as a function of temperature: melting layer signature, temperature
  dependence of the tables.

### Vertical cross-sections: `plot_cross_section.py`
```bash
python plot_cross_section.py FILE                                  # west-east line through the strongest cell
python plot_cross_section.py FILE --points 2.3 45.1 2.8 45.9 --contents
python plot_cross_section.py FILE --y 60 --vars Zh Zdr Kdp Rhohv --isotherms 0 -10 -20 -40
```
Radar variables (and contents with `--contents`) along the line, isotherms 0, -20 and -40 °C, and a small map
showing the line. Other positions: `--x KM`, `--ij I1 J1 I2 J2`.

### Compare two outputs: `plot_compare_outputs.py`
For two outputs of the **same model file** (old and new tables, two configurations, two bands, T-matrix and
Rayleigh):
```bash
python plot_compare_outputs.py old/dpolvar_X.nc new/dpolvar_X.nc --labels "old tables" "new tables"
```
Statistics of B - A (mean, RMSE, correlation) for all points and by temperature range (rain, melting layer,
ice), scatter plots, maps of the column maximum of A, B and B - A, B - A as a function of temperature.
If one output was computed on a subdomain (`subDomain`) and the other on the whole domain, the larger one
is cut to the subdomain automatically; `--box LONMIN LONMAX LATMIN LATMAX` compares only a lon/lat box.

### Storm swaths: `plot_max_over_time.py`
```bash
python plot_max_over_time.py "../operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20220818_0*.nc" --var Zh
python plot_max_over_time.py OBS/20220818_obs.nc --var Zdr --height 2 --box 4 11 41 45
```
Maximum of a variable over all the times, at a height or over the column, with the Zh = 40 dBZ contour of each
time (`--contour`): to compare the paths and intensity of simulated and observed cells.

# 4. Lookup tables

| Script | Use |
|---|---|
| `plot_tables_2D.py` | 2D maps of a whole table (content x temperature, Fw or N), with the areas where the size distribution is cut by the table (best overall check of a table) |
| `plot_D_M_varpol_alltypes.py` | curves of all types for one band, as a function of the content M or diameter D |
| `plot_D_M_varpol_Rayleigh.py` | one panel per type, several bands, T-matrix vs Rayleigh (Mie effects) |
| `plot_theoretical_varpol.py` | sensitivity to the table options (axis ratio, density, canting...) |

```bash
python plot_tables_2D.py --band C                                  # ICE3, all types, David2026PhD tables
python plot_tables_2D.py --band W --config David2025AMT --types ss --vars Zh Zdr --Tc -20
```
See [plot_tools/README.md](https://github.com/UMR-CNRM/operadar/blob/master/plot_tools/README.md) for the
options of each script and [Lookup table configurations](https://github.com/UMR-CNRM/operadar/wiki/Lookup-table-configurations).

# Reference

David, C., et al. (2025). Improved simulation of thunderstorm characteristics and polarimetric signatures with
LIMA two-moment microphysics in AROME. Atmos. Meas. Tech., 18, 3715–3745,
[https://doi.org/10.5194/amt-18-3715-2025](https://doi.org/10.5194/amt-18-3715-2025).
