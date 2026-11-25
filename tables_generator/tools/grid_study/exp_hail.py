"""Experiment 1b: diameter step of hail with a fine reference (0.002), (hail tables are not needed at Ku, K, Ka, W for now)."""
import json, sys
from study import *
bands = [b for b in sys.argv[1].split(",") if b in ("L", "S", "C", "X")]  # hail not needed at Ku, K, Ka, W for now
steps = (0.002, 0.005, 0.01)
jobs = []
for b in bands:
    p0 = load_default(b, "hh")
    for s in steps:
        for T in (-10.0, 0.0):
            jobs.append(dict(name=f"{b}|{s}|{T}", param=dict(p0, expDstep=s, Tcmin=T, Tcmax=T, ELEVmax=p0["ELEVmin"])))
res = run_all(jobs)
for b in bands:
    line = {"band": b, "type": "hh", "current": load_default(b, "hh")["expDstep"], "ref": steps[0], "steps": {}}
    for s in steps[1:]:
        worst = None
        for T in (-10.0, 0.0):
            e = errors(as_derived(res[f"{b}|{steps[0]}|{T}"]["rows"]), as_derived(res[f"{b}|{s}|{T}"]["rows"]))
            worst = e if worst is None else {k: (max(worst[k], e[k]) if k != "ok" else worst[k] and e[k]) for k in e}
        line["steps"][s] = worst
    print(json.dumps(line), flush=True)
