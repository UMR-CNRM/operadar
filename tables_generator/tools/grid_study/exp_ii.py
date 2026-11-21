import math, sys
from study import *
b = sys.argv[1]; h = sys.argv[2]
p0 = load_default(b, h)
steps = [0.00125, 0.0025, 0.005, 0.01, 0.02]
jobs = [dict(name=str(s), param=dict(p0, expDstep=s, Tcmin=-10.0, Tcmax=-10.0, ELEVmax=p0["ELEVmin"])) for s in steps]
res = run_all(jobs)
ref = as_derived(res[str(steps[0])]["rows"])
for s in steps[1:]:
    t = as_derived(res[str(s)]["rows"])
    # worst Zh error by concentration (P3) where ref Zh > ZMIN, and the M at which it happens
    by = {}
    for k, r in ref.items():
        if r is None or r["zref"] < ZMIN or t.get(k) is None: continue
        e = abs(t[k]["Zh"] - r["Zh"])
        if e > by.get(k[2], (0,))[0]: by[k[2]] = (e, k[3], r["zref"])
    print("step", s, "tmat %.1fs" % res[str(s)]["t_tmat"], " ".join("N=%.0e:%.2fdB@M=%.0e(Z=%.0f)" % (n, *v) for n, v in sorted(by.items())[::10]))
