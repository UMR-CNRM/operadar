"""Experiment 2: temperature step, emulating operadar's linear interpolation in Tc from a coarser grid."""
import json, sys
from study import *
bands = sys.argv[1].split(","); types = sys.argv[2].split(",")
jobs, meta = [], []
for b in bands:
    for h in types:
        p0 = load_default(b, h)
        # coarse diameter step: the temperature dependence does not depend on it
        rng = p0["expDmax"] - p0["expDmin"]
        s = max(p0["expDstep"], 0.02 if abs(rng / 0.02 - round(rng / 0.02)) < 1e-6 else p0["expDstep"])
        fw = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.25) if h == "wg" else {}
        jobs.append(dict(name=f"{b}{h}", param=dict(p0, expDstep=s, ELEVmax=p0["ELEVmin"], **fw)))
        meta.append((b, h, p0["Tcmin"], p0["Tcmax"]))
res = run_all(jobs)
for b, h, tmin, tmax in meta:
    rows = res[f"{b}{h}"]["rows"]
    if rows is None:
        print(json.dumps(dict(band=b, type=h, failed=res[f"{b}{h}"]["error"][-300:])), flush=True); continue
    ref = as_derived(rows)
    line = dict(band=b, type=h, steps={})
    for st in (2, 5, 10):
        if abs((tmax - tmin) / st - round((tmax - tmin) / st)) > 1e-6:
            continue
        test = interp_axis(rows, 0, lambda T, st=st, tmin=tmin: abs((T - tmin) / st - round((T - tmin) / st)) < 1e-6)
        line["steps"][st] = errors(ref, test)
    print(json.dumps(line), flush=True)
