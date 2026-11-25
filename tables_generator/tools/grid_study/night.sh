#!/bin/bash
# Overnight run: Fw (0.02/0.05 vs 0.01), Tc steps, coverage, then LIMA content/concentration steps and coverage.
# Each step is skipped if its result file is already complete (marker .done), so the script can be restarted.
cd "$(dirname "$0")"
PY=${PY:-python}; R=results; T="rr,ii,ss,gg,hh,wg,cl,cs"; ALL="L,S,C,X,Ku,K,Ka,W"
step() {  # step <name> <command...>
    local name=$1; shift
    [ -f $R/$name.done ] && return
    echo "$(date +%H:%M) start $name" >> $R/night.log
    "$@" > $R/$name.jsonl 2> $R/$name.err && touch $R/$name.done
    echo "$(date +%H:%M) end $name ($(wc -l < $R/$name.jsonl) lines)" >> $R/night.log
}
step fw4      $PY exp_fw.py  $ALL
step tc2      $PY exp_tc.py  $ALL $T
step cov2     $PY exp_cov.py $ALL $T
step m_lima   env STUDY_MICRO=LIMA $PY exp_m.py   $ALL $T
step cov_lima env STUDY_MICRO=LIMA $PY exp_cov.py $ALL $T
echo "$(date +%H:%M) all done" >> $R/night.log
