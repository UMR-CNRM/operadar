#!/bin/bash
# Figures of the AROME case recomputed with the current operadar code (David2026PhD tables, sensitivity threshold,
# no false echoes at 0 dBZ), from the model fields stored in the output of 20 Nov. (the .fa file needs epygram
# with the FA library). Run from the repository root:
#     bash docs/reports/figures_2026-11_outputs/make_figures_rerun.sh BAND [PYTHON]
set -e
B=$1; PY=${2:-python}
ROOT=$(pwd)
OUT=$ROOT/docs/reports/figures_2026-11_outputs/rerun_${B}
OLD=$ROOT/operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
RUN=$ROOT/operadarFiles/AROME/rerun_2026-11-23
F=$RUN/dpolvar_Arome_ICE3_${B}band_20250831_1800.nc
[ -f $F ] || $PY docs/reports/figures_2026-11_outputs/rerun_from_output.py $OLD $B $RUN/
cd $ROOT/plot_tools
$PY plot_output_check.py $F --out $OUT/check
$PY compute_cfad.py $F --label "cores (Zmax >= 40 dBZ)" --mask core --out $OUT/stats/core.nc
$PY compute_cfad.py $F --label "envelopes (Zmax >= 25 dBZ)" --mask envelope --out $OUT/stats/envelope.nc
$PY compute_cfad.py $F --label "all points" --mask all --out $OUT/stats/all.nc
$PY plot_cfad.py $OUT/stats/core.nc $OUT/stats/envelope.nc $OUT/stats/all.nc \
    --colors tab:red tab:orange tab:blue --figs summary --hmax 13 --name masks --out $OUT/cfad
$PY plot_cfad.py $OUT/stats/core.nc --figs freq --hmax 13 --name core --out $OUT/cfad
$PY plot_cross_section.py $F --points 2.3 45.1 2.8 45.9 --vars Zh Zdr Kdp Ah --contents --hmax 13 --out $OUT/sections
if [ "$B" = "C" ]; then
  $PY plot_compare_outputs.py $OLD $F --labels "20 Nov. run" "rerun 23 Nov." --vars Zh Zdr Kdp Ah --out $OUT/compare
else
  C=$RUN/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
  $PY plot_compare_outputs.py $C $F --labels "C band" "$B band" --vars Zh Zdr Kdp Ah --zh-min -10 --out $OUT/compare
  # CFADs of both bands together, in the cores of the C band (Zh rarely reaches 40 dBZ at Ka/W band)
  $PY compute_cfad.py $F --cells-from $C --label "$B band, C band cores" --mask core --out $OUT/stats/core_Ccells.nc
  $PY plot_cfad.py $OUT/stats/core_Ccells.nc --figs freq --hmax 13 --name core_Ccells --out $OUT/cfad
  $PY plot_cfad.py $ROOT/docs/reports/figures_2026-11_outputs/rerun_C/stats/core.nc $OUT/stats/core_Ccells.nc \
      --labels "C band" "$B band" --colors tab:blue tab:red --figs radar freq --diff --hmax 13 --name C_vs_${B}_core --out $OUT/cfad
fi
echo "done $B $(date '+%F %T')"
