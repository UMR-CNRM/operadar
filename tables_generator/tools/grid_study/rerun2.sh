#!/bin/bash
# Fixed content/concentration (m) and liquid-fraction (fw) experiments.
cd "$(dirname "$0")"
PY=${PY:-python}; R=results; T="rr,ii,ss,gg,hh,wg,cl,cs"; ALL="L,S,C,X,Ku,K,Ka,W"
$PY exp_m.py  $ALL $T > $R/m2.jsonl  2> $R/m2.err
$PY exp_fw.py $ALL    > $R/fw2.jsonl 2> $R/fw2.err
echo done > $R/RERUN2_DONE
