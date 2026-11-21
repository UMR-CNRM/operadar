"""Experiment 1: diameter step. Reference = current step / 2 (or finer); candidates = coarser steps."""
import json, sys, time
from study import *
bands = sys.argv[1].split(",")
types = sys.argv[2].split(",")
jobs, meta = [], []
for b in bands:
    for h in types:
        p0 = load_default(b, h)
        rng = p0["expDmax"] - p0["expDmin"]
        cur = p0["expDstep"]
        steps = sorted({round(s, 4) for s in (cur / 2, cur, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
                        if s >= min(cur / 2, cur) - 1e-9 and s >= 0.001 - 1e-9
                        and abs(rng / s - round(rng / s)) < 1e-6})
        Ts = (-10.0, 0.0) if h not in ("rr", "cl", "cs", "wg", "wh") else (0.0, 20.0) if h not in ("wg", "wh") else (0.0, 4.0)
        fw = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.25) if h in ("wg", "wh") else {}
        for s in steps:
            for T in Ts:
                p = dict(p0, expDstep=s, Tcmin=T, Tcmax=T, ELEVmax=p0["ELEVmin"], **fw)
                jobs.append(dict(name=f"{b}{h}|{s}|{T}", param=p))
        meta.append((b, h, cur, steps, Ts))
t = time.time()
res = run_all(jobs)
out = []
for b, h, cur, steps, Ts in meta:
    ref_s = steps[0]
    line = {"band": b, "type": h, "current": cur, "ref": ref_s, "steps": {}}
    for s in steps[1:]:
        worst = None
        for T in Ts:
            e = errors(as_derived(res[f"{b}{h}|{ref_s}|{T}"]["rows"]), as_derived(res[f"{b}{h}|{s}|{T}"]["rows"]))
            worst = e if worst is None else {k: (max(worst[k], e[k]) if k != "ok" else worst[k] and e[k]) for k in e}
        tm = sum(res[f"{b}{h}|{s}|{T}"]["t_tmat"] for T in Ts)
        line["steps"][s] = dict(worst, t_tmat=tm)
    out.append(line)
    print(json.dumps(line), flush=True)
print("elapsed %.0f s" % (time.time() - t), file=sys.stderr)
