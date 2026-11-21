# Plot tools

Scripts to look at the lookup tables and at the operadar outputs. Run them from this folder (`cd plot_tools`);
figures are written to `plot_tools/IMG/` (not tracked by git). They need `numpy`, `pandas` and `matplotlib`
(and `xarray`, `cartopy`, `epygram` for `plot_2Dmaps.py`).

| Script | What it shows | How it is configured |
|---|---|---|
| [`plot_tables_2D.py`](#plot_tables_2dpy-2d-maps-of-the-full-tables) | 2D maps of the full tables (all temperatures, contents, Fw, N) | command line |
| [`plot_D_M_varpol_alltypes.py`](#plot_d_m_varpol_alltypespy-all-types-on-one-figure) | curves of all types on one figure, one band | variables at the top of the file |
| [`plot_D_M_varpol_Rayleigh.py`](#plot_d_m_varpol_rayleighpy-several-bands-t-matrix-and-rayleigh) | one panel per type, several bands, T-matrix vs Rayleigh | variables at the top of the file |
| [`plot_theoretical_varpol.py`](#plot_theoretical_varpolpy-sensitivity-to-the-table-options) | sensitivity to the table options (axis ratio, density, Fw...) | command line |
| [`plot_2Dmaps.py`](#plot_2dmapspy-maps-of-operadar-outputs) | maps of the simulated radar variables (operadar output) | variables at the top of the file |
| `draft/` | CFADs and other work in progress | |

Tables are read from `../tables_generator/tables/<config>/`:
- `TmatCoefDiff_<band><type>`: single-particle values as a function of the diameter D (axis "D" in the scripts);
- `TmatCoefInt_<micro>_<n>M_<band><type>`: values integrated over the size distribution, as a function of the
  content M (axis "M"). The column `P3` is Fw for wet species and the concentration N for 2-moment species
  and pristine ice; `Mfrac` is the fraction of M inside the diameter range of the table.

Types: `rr` rain, `cl`/`cs` cloud water (land/sea), `ii` pristine ice, `ss` snow, `gg` graupel, `hh` hail,
`wg`/`wh` wet graupel/hail. Bands: L, S, C, X, Ku, K, Ka, W (hail is not computed at Ku, K, Ka, W).

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
python plot_tables_2D.py --band C                                   # ICE3, all types, default tables
python plot_tables_2D.py --band C --micro LIMA --types rr cl cs ii  # 2-moment species of LIMA
python plot_tables_2D.py --band W --types ss --vars Zh Zdr --Tc -20 # W-band snow, two variables
python plot_tables_2D.py --band C --config David2026PhD --types wg --Fw 0.5
```

Options: `--band`, `--config` (table folder, default `default`), `--micro` (ICE3, LIMA, ICJW, LIMC), `--types`,
`--vars` (Zh Zdr Kdp Rhohv Ah Zh-ZhR), `--elev` (closest table elevation, default 0), `--Tc` (temperature of the
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
| `TmatOption` | table folder per band (`default`, `David2026PhD`...) |
| `T_dict` | temperature of each type (hail: 0 °C at most) |
| `Fw_list` | Fw of the wet graupel/hail curves |
| `Nii`, `N2M` | concentration of pristine ice and of the other 2-moment species (LIMA rain, cloud) |
| `ELEVsel` | 0 for ground radars, 90 for vertically pointing radars |

Output: `IMG/DistTmat<band><micro><config>T...varpol_alltypes<config>.png`.

## `plot_D_M_varpol_Rayleigh.py`: several bands, T-matrix and Rayleigh

One panel per type, one colour per band, with the Rayleigh approximation as dashed lines (`plotR = True`). Useful
to see where Mie effects start for each band. Edit the top of the file (`micro`, `TmatOption`, `band_list`,
`typeh_list`, `plot_list`, `var_list`, `ELEVsel`, `T_dict`, `Fw_list`, `expN_list` for 2-moment species) and run
`python plot_D_M_varpol_Rayleigh.py`. The tables of every band of `band_list` must exist in the folder given by
`TmatOption` (e.g. `'C':'vertical'` needs `tables/vertical/`).

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
