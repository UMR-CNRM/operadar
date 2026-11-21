"""Integrated impact of canting=orient (orientation averaging, valid at any elevation) vs canting=analyt
(Ryzhkov et al. 2011 eq. 20, low elevation only) for the David2026PhD C-band configuration."""
import sys, math
from study import *
CASES = [("rr", 10.0, 0.0, [1e-3, 5e-3]), ("ii", -30.0, None, [1e-4]), ("gg", 0.0, 0.0, [1e-3, 5e-3]),
         ("wg", 0.0, 0.5, [1e-3, 5e-3])]
jobs = []
for h, T, fw, _ in CASES:
    p0 = load_default("C", h, "David2026PhD")
    extra = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.5) if h == "wg" else {}
    for cant in ("analyt", "orient"):
        jobs.append(dict(name=f"{h}|{cant}", param=dict(p0, canting=cant, Tcmin=T, Tcmax=T, ELEVmin=0.0, ELEVmax=90.0,
                                                       ELEVstep=30.0, **extra)))
res = run_all(jobs)
for h, T, fw, Ms in CASES:
    a, o = res[f"{h}|analyt"]["rows"], res[f"{h}|orient"]["rows"]
    for M in Ms:
        for e in (0.0, 30.0, 60.0, 90.0):
            def pick(rows):
                ks = [k for k in rows if k[1] == e and (fw is None or abs(k[2] - fw) < 1e-6)]
                if fw is None:  # ice: N closest to 1000 /m3
                    n = min({k[2] for k in ks}, key=lambda x: abs(math.log10(x) - 3)); ks = [k for k in ks if k[2] == n]
                k = min(ks, key=lambda k: abs(math.log10(k[3]) - math.log10(M)))
                return derived(rows[k])
            da, do = pick(a), pick(o)
            print(f"{h} T={T:5.1f} M={M:.0e} E={int(e):2d}: Zh {da['Zh']:6.2f}->{do['Zh']:6.2f}  Zdr {da['Zdr']:6.3f}->{do['Zdr']:6.3f}"
                  f"  Kdp {da['Kdp']:.3g}->{do['Kdp']:.3g}  rhohv {da['rhohv']:.4f}->{do['rhohv']:.4f}", flush=True)
