"""Elevation-step errors per elevation range, from results/elev_rows.pkl (written by exp_elev.py).
Kdp and Ah errors are given relative to the local value (as for the other axes) and relative to the 0 deg value."""
import pickle, sys
from study import *
rows_all = pickle.load(open("results/elev_rows.pkl", "rb"))
RANGES = [(0, 30), (30, 60), (60, 90)]
for name, rows in rows_all.items():
    if rows is None:
        print(name, "failed"); continue
    ref = as_derived(rows)
    out = []
    for st in (10, 15, 30):
        test = interp_axis(rows, 1, lambda e, st=st: abs(e / st - round(e / st)) < 1e-6)
        cells = []
        for (lo, hi) in RANGES:
            w = dict(Zh=0, Zdr=0, rho=0, Kdp=0, Kdp0=0)
            for k, r in ref.items():
                t = test.get(k)
                if t is None or r is None or r["zref"] < ZMIN or r["mfrac"] < MFRAC_MIN or not (lo < k[1] <= hi):
                    continue
                r0 = ref.get((k[0], 0.0) + k[2:])
                w["Zh"] = max(w["Zh"], abs(t["Zh"] - r["Zh"])); w["Zdr"] = max(w["Zdr"], abs(t["Zdr"] - r["Zdr"]))
                w["rho"] = max(w["rho"], abs(t["rhohv"] - r["rhohv"]))
                w["Kdp"] = max(w["Kdp"], abs(t["Kdp"] - r["Kdp"]) / (abs(r["Kdp"]) + FLOOR["Kdp"]))
                if r0: w["Kdp0"] = max(w["Kdp0"], abs(t["Kdp"] - r["Kdp"]) / (abs(r0["Kdp"]) + FLOOR["Kdp"]))
            cells.append(f"{lo}-{hi}: Zh {w['Zh']:.3f} Zdr {w['Zdr']:.3f} rho {w['rho']:.4f} Kdp {100*w['Kdp']:.0f}%/{100*w['Kdp0']:.1f}%0")
        out.append(f"  step {st}: " + " | ".join(cells))
    print(name); print("\n".join(out))
