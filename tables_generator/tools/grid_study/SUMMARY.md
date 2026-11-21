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
- snow (1M): M > 9.5e-3 kg/m3 (bias up to 6 dB at L band at 3e-2 kg/m3).
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
All the proposed steps were applied to every param file (`_default`, base, `_David2026PhD`, `_vertical`;
`.ref` files untouched) and `expMstep=0.02` in TmatrixIntegration.f90. Wet hail (not in the study) follows
wet graupel for Tcstep and Fwstep; its diameter step 0.02 passes at L-X (Zh <= 0.007 dB, Kdp <= 0.4 %, `results/dstep_wh.jsonl`).
