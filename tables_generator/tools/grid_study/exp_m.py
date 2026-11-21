"""Experiment 4: content (expM) and concentration (expCC) steps of TmatInt, emulating operadar's
linear interpolation of the linear coefficients in log10(M) and log10(N)."""
import json, math, sys
from study import *
jobs = []
for b in sys.argv[1].split(","):
    for h in sys.argv[2].split(","):
        p0 = load_default(b, h)
        fw = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.25) if h == "wg" else {}
        T = -10.0 if h in ("ii", "ss", "gg", "hh") else 0.0 if h == "wg" else 10.0
        jobs.append(dict(name=f"{b}{h}", param=dict(p0, Tcmin=T, Tcmax=T, ELEVmax=p0["ELEVmin"], **fw)))
res = run_all(jobs)
for b in sys.argv[1].split(","):
    for h in sys.argv[2].split(","):
        rows = res[f"{b}{h}"]["rows"]
        if rows is None:
            print(json.dumps(dict(band=b, type=h, failed=res[f"{b}{h}"]["error"][-300:])), flush=True); continue
        ref = as_derived(rows)
        zmin = -70.0 if h in ("cl", "cs") else ZMIN  # cloud: judged on attenuation, below -40 dBZ
        line = dict(band=b, type=h, micro=MICRO, steps={})
        for st in (0.02, 0.05, 0.1):
            test = interp_axis(rows, 3, lambda M, st=st: abs((math.log10(M) + 7) / st - round((math.log10(M) + 7) / st)) < 0.1)
            line["steps"]["M%g" % st] = errors(ref, test, zmin)
        if load_default(b, h)[MICRO] == "2M":  # P3 = concentration N
            # interpolation in log10(N): emulate with a log axis on P3
            for st in (0.2, 0.5):
                logged = {(k[0], k[1], math.log10(k[2]), k[3]): v for k, v in rows.items()}
                test = interp_axis(logged, 2, lambda ln, st=st: abs(ln / st - round(ln / st)) < 0.1)
                test = {(k[0], k[1], float(f"{10 ** k[2]:.3g}"), k[3]): v for k, v in test.items()}
                line["steps"]["N%g" % st] = errors(ref, test, zmin)
        print(json.dumps(line), flush=True)
