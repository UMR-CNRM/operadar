# Grid accuracy study: results and proposed steps (2026-11-21)

Method: tables built with the real Tmat/TmatInt (ice permittivity fix included), operadar's linear
interpolation emulated on a coarser grid and compared with a finer reference.
Tolerances: Zh 0.05 dB, Zdr 0.02 dB, rhohv 0.002, Kdp and Ah 1 % (floors 0.01 deg/km, 0.001 dB/km),
in cells with true Zh > -40 dBZ whose size distribution is covered by the table (Mfrac >= 0.99).
Raw results: `results/*.jsonl` (dstep_SCX2, dstep_high2, tc2, fw4, m2, m_lima, cov2, cov_lima).

## Proposed steps (hail: L-X only)

| | current | proposed | evidence (worst error at the proposed step) |
|---|---|---|---|
| expDstep rr, cl, cs | 0.02 | 0.02 | unchanged |
| expDstep ii, ss, gg (Ku-W) | 0.005 / 0.001 / 0.002 | 0.02 | Zh <= 0.017 dB (W ss 0.02) |
| expDstep hh (L-X) | 0.01 | 0.02 | Zh 0.004 dB (S-X); 0.05 also passes |
| expDstep wg | 0.005 (L-X), 0.002 (Ku-W) | 0.02 | Zh <= 0.005 dB |
| Tcstep ii, ss | 1 | 5 | Zh <= 0.004 dB |
| Tcstep gg, hh | 1 | 10 | Zh <= 0.006 dB |
| Tcstep rr, cl, cs | 1 | 5 (Ka rr: 2) | Ah <= 0.5 %; Ka rr Kdp 1.1 % at 5 |
| Tcstep wg | 1 | 2 (L-Ku), 1 (K, Ka, W) | at K-W, 2 gives Kdp 1-12 % |
| Fwstep wg | **0.05** | **0.01** | 0.05 gives Zh 0.4-0.7 dB, 0.02 still 0.16 dB (vs 0.01) |
| expMstep (TmatInt) | 0.01 | 0.02 | Zh <= 0.028 dB; Kdp 1.4-2.8 % only for Ka rr/wg |
| expCCstep (TmatInt, 2M) | 0.1 | 0.1 | 0.2 gives Zh 0.12-0.6 dB |

Fw is read as F4.2: 0.01 is the finest possible step.

## Low values (coverage of the size distribution)

The tables only under-estimate (by up to 40 dB) where the mean particle is larger than the table range:
- pristine ice: mean mass > 3e-6 kg (3 mg per crystal), e.g. N = 1000 /m3 with M > 3e-3 kg/m3;
- LIMA cloud: mean mass > 5e-8 kg (mean drop ~0.5 mm, i.e. drizzle, not cloud);
- snow (1M): M > 9.5e-3 kg/m3 (bias up to 6 dB at L band at 3e-2 kg/m3). Note: the ICE3 snow distribution
  is already partly outside the table (D > 50 mm) from 2-3 g/m3 (Mfrac 0.91 at 3 g/m3, 0.62 at 5 g/m3, C band);
  the Zh bias stays small up to ~9.5 g/m3 because these very large, low-density flakes scatter little.
  The 2D maps of plot_tools/plot_tables_2D.py show these areas (hatched where Mfrac < 0.99).
These combinations are not physical for the species, so no unrealistic low values are expected in normal
model output. Rain (ICE3 and LIMA), graupel, hail, wet graupel: no bias anywhere.
Small particles below 10 um (expDmin = -2) are cut, but only in cells below -40 dBZ.
Mfrac slightly above 1 (<= 1.017 in visible cells) where the PSD peak sits at expDmin: integration
error of the first diameter bins, negligible.

## Bugs found by the study
- Ice permittivity (Hufford 1991) used above 0 C: pole at 29.06 C (hail +11 dB at 29 C, no T-matrix
  convergence at Ku). Fixed: ice taken at min(T, 0 C).
- mixed_phase.py: wet hail never computed, hail counted twice when 'wh' enabled; Fw_pos counted rain twice
  with wet hail. Fixed.

## Applied (2026-11-21, approved)
All the proposed steps were applied to every param file (`_default` (now `_David2025AMT`), base, `_David2026PhD`, `_vertical` (removed since);
`.ref` files untouched) and `expMstep=0.02` in TmatrixIntegration.f90. Wet hail (not in the study) follows
wet graupel for Tcstep and Fwstep; its diameter step 0.02 passes at L-X (Zh <= 0.007 dB, Kdp <= 0.4 %, `results/dstep_wh.jsonl`).

## Elevation step (2026-11-21)
Tables now cover 0 to 90 deg (one table for ground, vertically pointing and spaceborne radars; the `vertical`
configuration was removed). Reference 0-90 every 5 deg, C, X, Ka, W bands (`exp_elev.py`, per-range errors with
`analyse_elev.py`, `results/elev_ranges.txt`). Kdp tends to 0 towards 90 deg, so its error is also given relative
to its 0 deg value.

| step | Zh | Zdr | Kdp (relative to 0 deg value) | elevations |
|---|---|---|---|---|
| 10 deg | <= 0.03 dB (wet graupel Ka/W 0.06-0.11) | <= 0.06 dB (rain, wet graupel C/X) | <= 1 % (wet graupel Ka 15 %) | 10 |
| 15 deg | x2 | x2 | x2 | 7 |
| 30 deg | up to 0.6 dB | up to 0.47 dB | up to 6 % (Ka wet graupel 70 %) | 4 |

Chosen: 0-90 deg every 10 deg in all David2025AMT and David2026PhD files. A 5 deg step would bring Zdr within
0.02 dB but doubles the size of the tables again (already x5 with 10 elevations instead of 2).
Before, tables had 0 and 10 deg only (L to X): operadar used the 10 deg values at higher elevations
(e.g. C-band rain at 5 g/m3 and 30 deg: Zdr 2.92 dB instead of 2.22 dB, Kdp +29 %).
