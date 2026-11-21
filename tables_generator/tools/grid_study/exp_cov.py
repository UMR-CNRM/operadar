"""Experiment 5: coverage of the size distribution by the diameter range of each default table.
Compares the current table with one whose diameter range is extended (expDmin-1, expDmax+ext;
drops of rain/cloud are not extended beyond 10 mm, where the T-matrix does not converge at W band)."""
import json, sys
from study import *
ext_hi = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
jobs, meta = [], []
for b in sys.argv[1].split(","):
    for h in sys.argv[2].split(","):
        p0 = load_default(b, h)
        s = 0.02 if p0["expDstep"] < 0.02 else p0["expDstep"]
        lo, hi = p0["expDmin"], p0["expDmax"]
        if abs((hi - lo) / s - round((hi - lo) / s)) > 1e-6: s = p0["expDstep"]
        fw = dict(Fwmin=0.0, Fwmax=1.0, Fwstep=0.5) if h == "wg" else {}
        T = -10.0 if h in ("ii", "ss", "gg", "hh") else 0.0 if h == "wg" else 10.0
        base = dict(p0, Tcmin=T, Tcmax=T, ELEVmax=p0["ELEVmin"], expDstep=s, **fw)
        jobs.append(dict(name=f"{b}{h}|cur", param=base))
        jobs.append(dict(name=f"{b}{h}|ext", param=dict(base, expDmin=lo - 1.0, expDmax=min(hi + ext_hi, 1.0) if h in ("rr", "cl", "cs") else hi + ext_hi)))
        meta.append((b, h))
res = run_all(jobs)
for b, h in meta:
    if res[f"{b}{h}|cur"]["rows"] is None or res[f"{b}{h}|ext"]["rows"] is None:
        print(json.dumps(dict(band=b, type=h, failed=(res[f"{b}{h}|cur"].get("error") or res[f"{b}{h}|ext"]["error"])[-300:])), flush=True); continue
    cur, ext = as_derived(res[f"{b}{h}|cur"]["rows"]), as_derived(res[f"{b}{h}|ext"]["rows"])
    by_p3 = {}
    for k, c in cur.items():
        e = ext.get(k)
        if c is None or e is None: continue
        by_p3.setdefault(k[2], []).append((k[3], c, e))
    out = dict(band=b, type=h, micro=MICRO, p3={})
    for p3, items in sorted(by_p3.items()):
        items.sort()
        bias = [(M, e["zref"], e["Zh"] - c["Zh"], c["mfrac"], e["mfrac"]) for M, c, e in items]
        visible = [x for x in bias if x[1] > ZMIN]
        worst = max(visible, key=lambda x: x[2], default=None)
        ok = [x for x in bias if x[2] < 0.05]
        out["p3"][p3] = dict(
            worst_visible_bias=worst,   # (M, Zh_true, bias dB, Mfrac cur, Mfrac ext)
            M_ok_min=min((x[0] for x in ok), default=None), M_ok_max=max((x[0] for x in ok), default=None),
            n_visible_biased=sum(1 for x in visible if x[2] > 0.05), n_visible=len(visible))
    print(json.dumps(out), flush=True)
