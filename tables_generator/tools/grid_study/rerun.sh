#!/bin/bash
# Reruns the experiments that failed in run_all.sh (Ku hail at 29 C before NPN1=200, W rain >10 mm, hail timeout).
cd "$(dirname "$0")"
PY=${PY:-python}; R=results; T="rr,ii,ss,gg,hh,wg,cl,cs"; ALL="L,S,C,X,Ku,K,Ka,W"
$PY exp_tc.py   $ALL $T                > $R/tc.jsonl    2> $R/tc.err
$PY exp_cov.py  $ALL $T                > $R/cov.jsonl   2> $R/cov.err
TMAT_TIMEOUT=7200 $PY exp_hail.py X,Ku > $R/hail.jsonl 2> $R/hail.err
echo done > $R/RERUN_DONE
