# Verification of operadar before the push of `fix-table-generator-bugs` (November 2026)

The whole package (new code and new lookup tables) is tested on the AROME file of the repository before the branch is pushed.
The new outputs are compared with the outputs of the former code (`master`) and the former tables.

| | |
|---|---|
| Model file | `modelFiles/AROME/historic.arome.franmg-01km30+0018:00.fa` (ICE3, 31 Aug. 2025 18 UTC, lon 2-4° E, lat 45-46.5° N) |
| Radar | C band, T-matrix, mixed phase `Fw_posg`, vertical geometry, attenuation on |
| Outputs | `operadarFiles/verif_2026-11-23/<run>/dpolvar_Arome_ICE3_Cband_20250831_1800.nc` |
| Figures | `docs/reports/figures_2026-11_operadar/*_final.png`, statistics `stats_final.json` |

**Reproduce** (from the repository root, environment with epygram):

```bash
python -m operadar historic.arome.franmg-01km30+0018:00.fa _verif_new_AMT     # same for the other configs
V=operadarFiles/verif_2026-11-23; f=dpolvar_Arome_ICE3_Cband_20250831_1800.nc
FIGTAG=_final python docs/reports/figures_2026-11_operadar/make_figures.py former=$V/old_code_old_tables/$f \
   new_oldtables=$V/new_code_old_tables/$f new_AMT=$V/new_code_AMT/$f new_PhD=$V/new_code_PhD/$f
```

## 1. Runs

| Label | Code | Tables | Config file |
|---|---|---|---|
| `former` | `master` | former tables (`tables/default_old`) | (master); bit-identical to the run of 20 Nov. |
| `new_oldtables` | branch (c8c545e) | former tables | `_verif_new_oldtables.py` |
| `new_AMT` | branch | `David2025AMT` (regenerated, C band) | `_verif_new_AMT.py` |
| `new_PhD` | branch | `David2026PhD` (C band) | `_verif_new_PhD.py` |
| `new_PhD_sens` | branch | `David2026PhD` + sensitivity threshold | `_verif_new_PhD_sens.py` |

All runs end without error (2-3 min each).

## 2. Summary

| Topic | Result |
|---|---|
| Refactoring (NetCDF reader, float64, reading only the needed columns) | Bit-identical to `master` when the bugs below are left in place |
| Bugs found and fixed (c8c545e) | (1) Interpolation corner index `np.floor` → `np.rint`: Zh too low by up to ~0.5 dB on ~30 % of the points, also on `master`. (2) Zh = 0 dBZ written where linear Z = 0 (6497 false echoes in the former output) → NaN. (3) rhohv of mixtures summed per species → covariances summed first, capped at 1. (4) Optional config keys: `test_AROME`/`test_MesoNH` did not load |
| New code, former tables | Zh +0.06 dB on average (p99 0.28 dB), from the interpolation fix; Zdr, Kdp, Ah unchanged (< 0.01) |
| **David2025AMT tables** | **Reproduce the former outputs**: Zh +0.08 dB on average (rain +0.06, ice +0.10, p99 0.55 dB), Zdr and Kdp differences < 0.01 dB / 0.01 °/km on 99 % of the points. The extra ice difference comes from the ice permittivity fix and the wider snow diameter range |
| **David2026PhD tables** | Physical differences, as expected from the thesis: ice Zdr +0.6 dB on average (median up to 2.5 dB at 12 km in the cores: oblate plates/aggregates, thesis p. 102 and 107), ice Zh +0.37 dB, ice Kdp +0.04 °/km; melting layer Zdr -0.16 dB and Kdp -0.12 °/km (canting of wet graupel); rain Zdr -0.04 dB |
| False 0 dBZ echoes | 6497 in `former`, **0** in all new runs |
| Sensitivity threshold | `new_PhD_sens` masks exactly the points with Zh_att < -30 dBZ (at 1 km, range-dependent); no other change |
| Convective cores (column max Zh >= 40 dBZ) | 2885 (former), 2902, 2905, 2840 columns: the core selection is stable |
| Not tested | Wet hail (no hail in this AROME file), MesoNH, real radar geometry |

**Conclusion:** the new code with the David2025AMT tables reproduces the former operadar to within 0.1 dB in Zh,
and the differences are explained by the bug fixes. The David2026PhD configuration gives the expected
differences in ice and in the melting layer. Nothing blocks the push.

## 3. Differences per region (run − former, all points with a valid Zh)

| Run | Variable | Rain (T>0) | Melting layer | Ice (T<0) |
|---|---|---|---|---|
| new_oldtables | Zh (dB) | +0.05 | +0.05 | +0.06 |
| new_AMT | Zh (dB) | +0.06 | +0.05 | +0.10 |
| new_PhD | Zh (dB) | +0.05 | -0.00 | +0.37 |
| new_oldtables | Zdr (dB) | +0.00 | +0.00 | +0.00 |
| new_AMT | Zdr (dB) | +0.00 | +0.00 | +0.00 |
| new_PhD | Zdr (dB) | -0.04 | -0.16 | +0.61 |
| new_oldtables | Kdp (°/km) | +0.00 | +0.01 | +0.00 |
| new_AMT | Kdp (°/km) | +0.00 | +0.00 | +0.00 |
| new_PhD | Kdp (°/km) | -0.01 | -0.12 | +0.04 |

Ah differences are below 0.001 dB/km everywhere. Full statistics (mean absolute difference, 1st/99th
percentiles, maximum) in `stats_final.json`.

## 4. Figures

| Figure | Content |
|---|---|
| `profiles_cores_final.png` | Median and Q25-Q75 profiles of Zh, Zdr, Kdp in the convective cores (David et al. 2025, AMT, Fig. 6). `former`, `new_oldtables` and `new_AMT` overlap; `new_PhD` differs above the melting layer (Zdr) |
| `cfad_{Zh,Zdr,Kdp}_final.png` | Frequency CFADs in the cores and differences with `former` (a few % around the melting layer for `new_AMT`, up to 25 % in ice for `new_PhD`) |
| `scatter_*_final.png` | Point-by-point comparison with `former` |
| `maps_Zmax_final.png`, `section_final.png`, `contents_cores_final.png` | Maps of column max Zh, cross-section, hydrometeor contents in the cores |

## 5. Lookup tables, state on 24 Nov.

All tables of the queue of 23 Nov. are regenerated (David2025AMT C; David2026PhD L, S, C, X, Ku, K, Ka, W)
except **W band snow**: T-matrix does not converge for snowflakes larger than ~91 mm at 3.15 mm wavelength
(`CONVERGENCE IS NOT OBTAINED FOR NPN1=160`). Tests at -75 and +10 °C: `expDmax` 1.96 converges, 1.98 fails at
+10 °C, 2.0 fails everywhere. `TmatParam_Wss` and `TmatParam_Wss_David2026PhD` now use `expDmax:1.96`
(91 mm instead of 100 mm; the former W tables used 1.7 = 50 mm). The impact should be very small,
since snowflakes larger than 91 mm are very rare in the ICE3/LIMA size distributions (not checked
quantitatively).
The W snow tables were regenerated with `expDmax:1.96` on 24 Nov. (exit 0), and the test
`test_snow_tables_reach_100_mm` now requires 1.96 at W band.

## 6. Update of 24 Nov.: W band, MesoNH

Details and figures in `2026-11_operadar_output_analysis.md`, sections 6-8.

* **AROME, C and W band from the .fa file** (current code, David2026PhD, attenuation): C band identical to
  `new_PhD` with the sensitivity threshold; Zh_att changes exactly like Zh compared with the 20 Nov. run (the
  attenuation code is unchanged). W band: all differences with C band are those expected in the Mie regime
  (Zh saturation, Zdr of rain ~0, negative Kdp in rain, strong attenuation).
* **MesoNH reading was broken, also on master** (invalid annotation in `with_netCDF4.py`, functions called
  without their module in `mesonh.py`): no MesoNH file could be read. Fixed, with three other bugs of the same
  reader (last row/column of the domain and of the subdomain dropped, water vapour `qv` not read). New tests
  `tests/test_read_mesonh.py`.
* **MesoNH 3BL13** (cyclone near La Réunion, subdomain 53-60° E, 24-17° S): checks pass, differences with the
  20 Nov. output only in the melting layer (Zdr -0.08 dB, Kdp -0.07 °/km). The whole 1202 x 1202 x 82 domain
  does not fit in memory (23 GB available).
* The former outputs (AROME of 20 Nov., MesoNH) contain false echoes at exactly 0 dBZ (6497 and 1 365 628
  points); `plot_output_check.py` now detects them.

## 7. Update of 24 Nov. (evening): ICJW and LIMC, MesoNH ice concentration

Two more bugs, also on master:

* **1-moment ice of ICJW and LIMC was interpolated with the wrong variable** (`formats_data.Fw_or_Nc`): only
  `'ICE3'` and a misspelled `'IJW'` used the concentration for 1-moment pristine ice, so ICJW and LIMC ice was
  looked up with the liquid fraction Fw (≈ 0) on the concentration axis of the tables. Ice now always uses the
  concentration, for all schemes. Test in `tests/test_dualpol_variables.py`.
* **MesoNH 1-moment ice concentration multiplied by the air density**: `CIT` is already in m⁻³, but the reader
  multiplied it by ρ (correct only for the 2-moment concentrations `CRAIN`/`CICE`/`CCLOUD`, in kg⁻¹). The
  concentrations are now chosen from the number of moments of each species, not from the scheme name: 1-moment
  ice = `CIT` (m⁻³, or 800 m⁻³ if absent), 2-moment species = C × ρ. LIMC ("LIMA chaud") is 2-moment for
  cloud water and rain only; its ice is 1-moment. Tests in `tests/test_read_mesonh.py`.

**Impact of the CIT fix on MesoNH 3BL13 (ICE3, C band, same subdomain as section 6):** the ice concentration
increases by a factor 3.3 (median, = 1/ρ aloft), the ice crystals get smaller and Zh decreases with altitude:

| Altitude | 0-4 km | 4-6 km | 6-8 km | 8-10 km | 10-12 km | 12-14 km | 14-16 km |
|---|---|---|---|---|---|---|---|
| ΔZh (dB) | 0.00 | -0.18 | -0.57 | -1.16 | -2.55 | -4.65 | -6.13 |

Mean over all points: Zh -1.41 dB, Zdr +0.07 dB, Kdp +0.002 °/km; ice-dominated points: Zh -5.4 dB.
The output of section 6 was computed before this fix (output with the fix:
`operadarFiles/MesoNH/2026-11-24_CITfix/`).

**ICJW and LIMC check (C band).** Integrated tables were computed for ICJW and LIMC from the David2026PhD C-band
scattering tables (test tables in `tables_generator/tables/_test_ICJW_LIMC/`, not delivered). The four schemes
were run on the same AROME fields (ICE3 output of the .fa file; for the 2-moment species the concentrations
are derived from the ICE3 size distributions, ice at 800 m⁻³), one output per species:

* ICJW = ICE3 and LIMC = LIMA bit for bit for cloud, rain, ice, graupel and wet graupel;
* ICJW snow = LIMC snow bit for bit;
* the snow differs, as expected: ICJW/LIMC use the Wurtz (2023) snow (slope λ function of T only, N0 ∝ content),
  ICE3/LIMA the ICE3 snow (λ from the content). At very small contents (~10⁻⁴ g m⁻³, -55 °C) the ICE3 snow
  gives Zh of -100 to -400 dBZ (NaN on 1831 points where Z underflows), the Wurtz snow about -30 dBZ; at high
  contents (1-10 g m⁻³) the Wurtz snow is weaker (median 14 dBZ instead of 32 dBZ).

ICJW on the real .fa file (C band, David2026PhD scattering, sensitivity threshold;
`operadarFiles/AROME/fa_ICJW_2026-11-24/`) vs ICE3 (section 6): contents and concentrations identical, Zh +0.20 dB
on average (0 below 2 km, +0.6 dB at 2-4 km, -1.5 dB at 6-8 km, +1.4 to +2.1 dB above 8 km), Zdr -0.05 dB,
Kdp +0.003 °/km; 22 172 more points above the sensitivity threshold (weak snow). Points with Zh > 40 dBZ unchanged.
