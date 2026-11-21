This page summarises the predefined configurations of the lookup tables (`David2026PhD`, `David2025AMT`), what each one is for, how they differ, and the points still to be checked. All values come from the parameter files `tables_generator/param/TmatParam_<band><type>_<config>` (state of November 2026).

# 1. Overview

| Configuration | Bands | Purpose | Command |
|---|---|---|---|
| `David2026PhD` (**default**) | L, S, C | Options of C. David's PhD thesis (2026, configuration NEW of Table 4.3): canting of all ice and rain particles, oblate pristine ice (AR 0.2) and snow (AR 0.6) | `./generate_tables.sh --David2026PhD --band C -j 8` or `./generate_tables.sh --default --band C -j 8` |
| `David2025AMT` | L, S, C, X, Ku, K, Ka, W | Options of [David et al. (2025), AMT](https://doi.org/10.5194/amt-18-3715-2025) (configuration REF of the thesis): no canting, spherical pristine ice | `./generate_tables.sh --David2025AMT --band X -j 8` |
| *base files* `TmatParam_<band><type>` | all | Templates copied by `--newConf <name>` for a new configuration (see section 5) | `./generate_tables.sh --newConf myConf --band C` |

`--default` is an alias for `--David2026PhD`: it uses the `TmatParam_*_David2026PhD` files and writes the tables to `tables_generator/tables/David2026PhD/`. For the bands without David2026PhD files (X, Ku, K, Ka, W) the script stops and points to `--David2025AMT`.

All configurations cover **elevations from 0° to 90° every 10°**, so the same tables serve ground radars at any elevation, vertically pointing radars and spaceborne radars (there is no separate `vertical` configuration any more). Hail and wet hail are not computed at Ku, K, Ka and W bands, and wet hail is not computed at all for now (`HYDRO_LIST` in `generate_tables.sh`).

Wavelengths (`LAM`, mm): L 190.4, S 106.2, C 53.2, X 31.9, Ku 20.0, K 15.0, Ka 8.4, W 3.15.

# 2. Options common to both configurations

| Type | Dielectric function | Axis ratio function | Density | Diameter range (mm) |
|---|---|---|---|---|
| `rr` rain | `Liebe91` | `BR02` (Brandes et al. 2002) | water | 0.01-10 |
| `cl`, `cs` cloud | `Liebe91` | spheres | water | 0.01-2 |
| `ii` pristine ice | `RY19dry` | `CNST` (section 3) | `RHOX` (m(D) law, capped at 916 kg/m³) | 0.01-10 (0.01-9.55 for David2025AMT at Ku band: the T-matrix of a 10 mm ice sphere does not converge at Ku band) |
| `ss` snow | `RY19dry` | `AUds` or `CNST` (section 3) | `RHOX` | 0.01-50 |
| `gg` graupel | `RY19dry` | `RYdg` (Ryzhkov et al. 2011, eq. 13): 1 - 0.02 D for D ≤ 10 mm, 0.8 above (note: David et al. 2025 write "set to 0.8", but this law was used) | `RHOX` | 0.1-50 (0.05-50 at Ku-W) |
| `hh` hail | `RY19dry` | `RYdg`, 0.8 for D > 10 mm | `RHOX` | 0.1-100 |
| `wg` wet graupel | `LBwetgr` (Le Bastard 2019) | `RYwg` (Ryzhkov et al. 2011) | `RHOX` + melt water | 0.1-50 (0.05-31.6 at Ku-W) |

Temperature ranges (°C):

| Type | David2025AMT | David2026PhD |
|---|---|---|
| rain, cloud | -20 to 40 | -60 to 40 |
| pristine ice, snow | -75 to 10 | -75 to 10 |
| graupel | -40 to 10 | -60 to 10 (L, C), -40 to 10 (S) |
| hail | -40 to 0 | -40 to 30 (section 5) |
| wet graupel | -40 to 10 | -40 to 10 (C, S), -60 to 30 (L) |

operadar uses the values at the closest table bound outside these ranges.

Grid steps: elevation 10°; diameter 0.02 in log10(D); temperature 5 °C (ice, snow, rain, cloud; 2 °C for Ka-band rain), 10 °C (graupel, hail), 2 °C for wet graupel (1 °C at K, Ka, W); Fw 0.01. They were chosen with a grid accuracy study (interpolation errors ≤ 0.05 dB on Zh, ≤ 0.02 dB on Zdr (0.06 dB along elevation), ≤ 1 % on Kdp), see `tables_generator/tools/grid_study/SUMMARY.md`.

Particle size distributions: ICE3 (1 moment) for all types; LIMA with 2 moments for rain, cloud and pristine ice.

# 3. Options that differ: canting and particle shape

`nocant` = no canting (all particles aligned, horizontal major axis); a number = standard deviation (°) of the angle between the symmetry axis and the vertical, with `canting=orient`: the scattering is averaged over the orientations (Gaussian distribution around the vertical, uniform azimuth), which is valid at any elevation. The former option `canting=analyt` (analytical angular moments of Ryzhkov et al. 2011, eq. 20) is only valid at low elevation; it gives the same results as `orient` at 0° but overestimates rhohv (and the Zh of wet graupel) at higher elevations. AR = axis ratio (minor/major); for `AUds` and `RYdg` the value is the axis ratio of large particles (D > 8 mm for snow, D > 10 mm for graupel and hail), small particles being closer to spheres.

Both configurations use the same options at all their bands (checked by the tests).

| Type | `David2025AMT` (REF) | `David2026PhD` (NEW) |
|---|---|---|
| rain | nocant | 10° |
| pristine ice | nocant, **spheres** (AR 1.0) | 10°, plates AR 0.2 |
| snow | nocant, `AUds` 0.75 | 40°, `CNST` 0.6 |
| graupel | nocant | 40° |
| hail (not in the publications: as graupel) | nocant | 40° |
| wet graupel | nocant | 10° to 40° depending on Fw: σ = σg + Fw (σr − σg) (thesis eq. 4.1) |
| cloud | spheres | spheres |

# 4. Impact of these choices

Values at 0° elevation (ICE3 PSD), David2025AMT → David2026PhD (same parameters measured with the November 2026 tables; C-band snow and hail not compared).

| Type | Conditions | C band | L band |
|---|---|---|---|
| rain | 10 °C, 1 g/m³ | Zdr 1.36 → 1.24 dB, Kdp 0.62 → 0.56 °/km | Zdr 1.29 → 1.18 dB |
| rain | 10 °C, 5 g/m³ | Zdr 3.02 → 2.75 dB, Kdp 5.9 → 5.4 °/km, rhohv 0.967 → 0.973 | Zdr 2.24 → 2.03 dB |
| pristine ice | -30 °C, 0.1 g/m³, N = 1000 m⁻³ | Zdr 0.00 → 1.18 dB, Kdp 0 → 0.064 °/km | Zdr 0.00 → 1.18 dB |
| graupel | 0 °C, 5 g/m³ | Zdr 0.07 → 0.02 dB, Kdp 0.11 → 0.03 °/km | same |
| hail | -10 °C, 5 g/m³ | | Zdr 0.88 → 0.22 dB, Kdp 0.33 → 0.09 °/km |
| wet graupel | 0 °C, Fw 0.5, 5 g/m³ | Zdr 2.30 → 1.28 dB, Kdp 26 → 15 °/km, rhohv 1.000 → 0.996 | Zdr 2.21 → 1.24 dB |

Zh changes by less than 0.4 dB in all cases. The choices that matter most are the shape of pristine ice (Zdr 0 dB for spheres, about 1.2 dB for plates of AR 0.2) and the canting of wet graupel and hail (Zdr and Kdp divided by about 2 to 4). Without canting, rhohv stays at 1 for ice species.

Elevation: Zdr and Kdp decrease with elevation (rain at C band, 5 g/m³: Zdr 3.02 dB at 0°, 2.22 dB at 30°, 0.70 dB at 60°, 0 at 90°).

# 5. Points still to check

1. **Hail up to 30 °C and wet hail up to 40 °C in David2026PhD**: solid ice is computed at 0 °C above 0 °C, so these temperatures only add identical rows; `Tcmax 0` as in David2025AMT gives the same results.
2. **The base files used by `--newConf`** are neither David2025AMT nor David2026PhD (e.g. `TmatParam_Crr` has canting 10°, `TmatParam_Xii` has `Frim 0.0`), and they still have the old elevations. A user who creates a new configuration should start from a copy of one of the predefined configurations. *Proposal: make the base files a copy of David2026PhD (L, S, C) and David2025AMT (other bands).*
3. **ICE3 snow at high content**: the snow size distribution goes beyond the table maximum (50 mm) above 2-3 g/m³ (fraction of the mass inside the table 0.91 at 3 g/m³, 0.62 at 5 g/m³ at C band). If simulations have more snow than that, increase `expDmax` for snow.
4. **Size of the tables**: with 10 elevations, the tables are 5 times larger than with 2 (94 GB for David2025AMT with its 8 bands, 39 GB for David2026PhD with L, S, C; 5 to 15 GB per band). operadar reads each table entirely (all 26 columns); reading only the columns it needs would divide the memory by about 3.

# 6. Note on the rhohv column of the tables

In cells with an extremely small signal (Zh below about -150 dBZ), the `rhohv` and `rhohvR` columns of the `TmatCoefInt_*` files contain `NaN` or `Infinity` (single-precision underflow of σhh·σvv). operadar is not affected: it reads σhh, σvv and δco and computes rhohv itself. Only plots of the `rhohv` column of the tables see these values (`plot_tables_2D.py` ignores them).
