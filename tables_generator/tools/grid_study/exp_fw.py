"""Experiment 3: liquid-fraction step of wet graupel (operadar interpolates linearly in Fw). Reference step 0.01 (Fwstep is read as F4.2: 0.025 would become 0.03)."""
import json, sys
from study import *
jobs = []
for b in sys.argv[1].split(","):
    p0 = load_default(b, "wg")
    for T in (0.0, 4.0):
        jobs.append(dict(name=f"{b}|{T}", param=dict(p0, Fwstep=0.01, Tcmin=T, Tcmax=T, ELEVmax=p0["ELEVmin"], expDstep=0.02 if b in "LSCX" else 0.01)))
res = run_all(jobs)
for b in sys.argv[1].split(","):
    line = dict(band=b, type="wg", steps={})
    for st in (0.02, 0.05):
        worst = None
        for T in (0.0, 4.0):
            rows = res[f"{b}|{T}"]["rows"]
            e = errors(as_derived(rows), interp_axis(rows, 2, lambda f, st=st: abs(f / st - round(f / st)) < 1e-6))
            worst = e if worst is None else {k: (max(worst[k], e[k]) if k != "ok" else worst[k] and e[k]) for k in e}
        line["steps"][st] = worst
    print(json.dumps(line), flush=True)
