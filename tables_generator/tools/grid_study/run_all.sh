#!/bin/bash
# Runs the remaining grid-study experiments one after the other (each uses 8 processes).
cd "$(dirname "$0")"
PY=${PY:-python}; R=results; T="rr,ii,ss,gg,hh,wg,cl,cs"; ALL="L,S,C,X,Ku,K,Ka,W"
$PY exp_dstep.py S,C,X rr,ii,ss,gg,hh,wg > $R/dstep_SCX2.jsonl 2> $R/dstep_SCX2.err
$PY exp_hail.py X            > $R/hail.jsonl  2> $R/hail.err
$PY exp_tc.py   $ALL $T                > $R/tc.jsonl    2> $R/tc.err
$PY exp_fw.py   $ALL                   > $R/fw.jsonl    2> $R/fw.err
$PY exp_m.py    $ALL $T                > $R/m.jsonl     2> $R/m.err
$PY exp_cov.py  $ALL $T                > $R/cov.jsonl   2> $R/cov.err
echo done > $R/ALL_DONE
