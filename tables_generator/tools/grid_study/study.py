"""Accuracy-vs-cost study of the lookup-table grids (run with the test venv python).

Every experiment builds small tables with the real Tmat/TmatInt (parallel processes),
emulates what operadar does with them (linear interpolation of the linear coefficients
in Tc, Fw, log10 M) and compares the operadar variables with a finer reference.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import sys
import tempfile
import time
import uuid
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/dir/tests/tables_generator")
from harness import PARAM, PARAM_LAYOUT, Sandbox  # noqa: E402

BIN = Path(os.environ.get("STUDY_BIN", "/dir/tables_generator/src"))
WORK = Path(os.environ.get("STUDY_WORK", "/tmp/grid_study_work"))
TOL = dict(Zh=0.05, Zdr=0.02, rhohv=0.002, Kdp=0.01, Ah=0.01)  # dB, dB, -, relative, relative
FLOOR = dict(Kdp=1.0, Ah=0.1)  # error = |d|/(|ref|+floor) <= 1%: i.e. 1% or 0.01 deg/km (Kdp), 0.001 dB/km (Ah)
ZMIN = -40.0  # errors evaluated where the reference Zh (dBZ) is above this value
COLS = ("Tc", "ELEV", "P3", "M", "sighh", "sigvv", "REdeltaco", "IMdeltaco", "kdp", "Ah", "zhh", "Mfrac")
MICRO = os.environ.get("STUDY_MICRO", "ICE3")  # microphysics scheme of the TmatInt tables
MFRAC_MIN = 0.99  # cells where the PSD is truncated are excluded from step comparisons


def load_default(band, typeh, suffix="default"):
    d = {}
    for l in (PARAM / f"TmatParam_{band}{typeh}_{suffix}").read_text().splitlines():
        k, _, v = l.partition(":")
        d[k] = v
    p = {k: (d[k] if f.startswith("A") else float(d[k])) for k, f in PARAM_LAYOUT}
    for k in ("ICE3", "LIMA", "ICJW", "LIMC"):
        p[k] = d.get(k, "1M")
    return p


def _run(job):
    """Run Tmat + TmatInt for one parameter dict; return {key: values} of the integrated table
    (rows=None and the error message if Tmat or TmatInt fails)."""
    try:
        return _run1(job)
    except Exception as e:
        msg = [l for l in str(e).splitlines() if "CONVERG" in l or "failed" in l or "Timeout" in l]
        print(f"FAILED {job['name']}: {' | '.join(msg[-2:]) or str(e)[:200]}", file=sys.stderr, flush=True)
        return dict(name=job["name"], rows=None, error=str(e)[-500:])


def _run1(job):
    p = job["param"]
    root = WORK / uuid.uuid4().hex
    try:
        sb = Sandbox(root, BIN / "Tmat", BIN / "TmatInt")
        t0 = time.time()
        diff = sb.run_tmat(p)
        t1 = time.time()
        out = sb.run_tmatint(diff, "t", p["type"], p["band"], micro=MICRO, moment=p.get(MICRO, "1M"))
        t2 = time.time()
        lines = out.read_text().splitlines()
        names = lines[2].split(";")
        idx = [names.index(c) for c in COLS]
        rows = {}
        for l in lines[3:]:
            v = l.split(";")
            vals = [float(v[i]) for i in idx]
            rows[(round(vals[0], 2), round(vals[1], 2), float(f"{vals[2]:.3g}"), float(f"{vals[3]:.4g}"))] = vals[4:]
        return dict(name=job["name"], rows=rows, t_tmat=t1 - t0, t_int=t2 - t1)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def run_all(jobs, workers=8):
    WORK.mkdir(parents=True, exist_ok=True)
    import multiprocessing
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("fork")) as ex:
        return {r["name"]: r for r in ex.map(_run, jobs)}


def derived(v):
    sighh, sigvv, re, im, kdp, ah, zhh, mfrac = v
    if sighh <= 0 or sigvv <= 0:
        return None
    return dict(Zh=10 * math.log10(sighh), Zdr=10 * math.log10(sighh / sigvv),
                rhohv=4 * math.pi * math.hypot(re, im) / math.sqrt(sighh * sigvv), Kdp=kdp, Ah=ah, zref=zhh,
                mfrac=mfrac)


def errors(ref, test, zmin=ZMIN):
    """Max error per variable between two dicts of derived values (same keys), where ref Zh > zmin."""
    worst = {k: 0.0 for k in TOL}
    n = 0
    for key, r in ref.items():
        t = test.get(key)
        if t is None or r is None or r["zref"] < zmin or r["mfrac"] < MFRAC_MIN:
            continue
        n += 1
        worst["Zh"] = max(worst["Zh"], abs(t["Zh"] - r["Zh"]))
        worst["Zdr"] = max(worst["Zdr"], abs(t["Zdr"] - r["Zdr"]))
        worst["rhohv"] = max(worst["rhohv"], abs(t["rhohv"] - r["rhohv"]))
        for k in ("Kdp", "Ah"):
            worst[k] = max(worst[k], abs(t[k] - r[k]) / (abs(r[k]) + FLOOR[k]))
    worst["n"] = n
    worst["ok"] = all(worst[k] <= TOL[k] for k in TOL)
    return worst


def as_derived(rows):
    return {k: derived(v) for k, v in rows.items()}


def interp_axis(rows, axis, keep):
    """Emulate a coarser grid along one axis (0 Tc, 2 P3=Fw, 3 M as log10): keep only the values for which
    keep(value) is True (value as in the table, not log10) and rebuild the others by linear interpolation of the linear coefficients
    between the neighbouring kept values (what operadar does). Returns derived values at all points."""
    groups = {}
    for k, v in rows.items():
        other = tuple(x for i, x in enumerate(k) if i != axis)
        groups.setdefault(other, []).append((k[axis], k, v))
    out = {}
    for other, items in groups.items():
        items.sort()
        xs = [(math.log10(x) if axis == 3 else x) for x, _, _ in items]
        kept = [i for i, (x, _, _) in enumerate(items) if keep(x)]
        for i, (x, k, v) in enumerate(items):
            j = max([kk for kk in kept if kk <= i], default=None)
            l = min([kk for kk in kept if kk >= i], default=None)
            if j is None or l is None:
                continue
            if j == l:
                out[k] = derived(v)
                continue
            w = (xs[i] - xs[j]) / (xs[l] - xs[j])
            vj, vl = items[j][2], items[l][2]
            out[k] = derived([a + w * (b - a) for a, b in zip(vj, vl)][:6] + [v[6], v[7]])
    return out
