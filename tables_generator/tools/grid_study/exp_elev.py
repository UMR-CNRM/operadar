"""Experiment: elevation step. Reference = 0 to 90 deg every 5 deg; candidates 10, 15, 30, 45 deg, emulating
operadar's linear interpolation in elevation (on the linear coefficients, as for the other axes)."""
import json, sys, pickle
from study import *
bands = sys.argv[1].split(","); types = sys.argv[2].split(",")
TC = {"rr": 10.0, "cl": 10.0, "cs": 10.0, "ii": -30.0, "ss": -10.0, "gg": 0.0, "hh": 0.0, "wg": 0.0}
jobs, meta = [], []
for b in bands:
    for h in types:
        p0 = load_default(b, h)
        fw = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.25) if h == "wg" else {}
        p = dict(p0, ELEVmin=0.0, ELEVmax=90.0, ELEVstep=5.0, Tcmin=TC[h], Tcmax=TC[h], **fw)
        jobs.append(dict(name=f"{b}{h}", param=p))
        meta.append((b, h))
res = run_all(jobs)
pickle.dump({k: v['rows'] for k, v in res.items()}, open('results/elev_rows.pkl', 'wb'))
for b, h in meta:
    rows = res[f"{b}{h}"]["rows"]
    if rows is None:
        print(json.dumps(dict(band=b, type=h, failed=True)), flush=True); continue
    ref = as_derived(rows)
    line = dict(band=b, type=h, steps={})
    for st in (10, 15, 30, 45):
        test = interp_axis(rows, 1, lambda e, st=st: abs(e / st - round(e / st)) < 1e-6)
        line["steps"][st] = errors(ref, test)
    print(json.dumps(line), flush=True)
