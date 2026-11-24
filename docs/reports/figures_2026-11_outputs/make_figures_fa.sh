#!/bin/bash
# Figures of section 6-8 of 2026-11_operadar_output_analysis.md: AROME case recomputed from the .fa file with the
# current code (C and W band, with attenuation) and MesoNH case (3BL13, cyclone near La Reunion).
# The operadar runs (environment with epygram + FA library) are:
#     python -m operadar "historic.arome.franmg-01km30+0018:00.fa" _verif_fa_C     (and _verif_fa_W)
#     python -m operadar 3BL13.1.MXX25.004.nc _verif_mesonh_C   (subDomain 53-60 E, 24-17 S)
# Run from the repository root:   bash docs/reports/figures_2026-11_outputs/make_figures_fa.sh [PYTHON] [PART]
# PART: arome, mesonh or all (default)
set -e
PY=${1:-python}; PART=${2:-all}
ROOT=$(pwd)
FIG=$ROOT/docs/reports/figures_2026-11_outputs
OLD=$ROOT/operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
FA=$ROOT/operadarFiles/AROME/fa_2026-11-24
RERUN=$ROOT/operadarFiles/AROME/rerun_2026-11-23
C=$FA/dpolvar_Arome_ICE3_Cband_20250831_1800.nc
W=$FA/dpolvar_Arome_ICE3_Wband_20250831_1800.nc
MNH_OLD=$ROOT/modelFiles/MesoNH/dpolvar_MesoNH_ICE3_Cband_20240115_0400.nc
MNH=$ROOT/operadarFiles/MesoNH/2026-11-24/dpolvar_MesoNH_ICE3_Cband_20240115_0400.nc
cd $ROOT/plot_tools

if [ "$PART" = all ] || [ "$PART" = arome ]; then
  # ---- AROME C band from the .fa file
  O=$FIG/fa_C
  $PY plot_output_check.py $C --out $O/check
  $PY compute_cfad.py $C --label "cores (Zmax >= 40 dBZ)" --mask core --vars Zh Zh_att Zdr Kdp --out $O/stats/core.nc
  $PY compute_cfad.py $C --label "envelopes (Zmax >= 25 dBZ)" --mask envelope --out $O/stats/envelope.nc
  $PY compute_cfad.py $C --label "all points" --mask all --out $O/stats/all.nc
  $PY plot_cfad.py $O/stats/core.nc $O/stats/envelope.nc $O/stats/all.nc \
      --colors tab:red tab:orange tab:blue --figs summary --hmax 13 --name masks --out $O/cfad
  $PY plot_cfad.py $O/stats/core.nc --figs freq --hmax 13 --name core --out $O/cfad
  $PY plot_cross_section.py $C --points 2.3 45.1 2.8 45.9 --vars Zh Zh_att Zdr Kdp Ah --contents --hmax 13 --out $O/sections
  $PY plot_compare_outputs.py $OLD $C --labels "20 Nov. run" "fa 24 Nov." --vars Zh Zdr Kdp Ah Zh_att --out $O/compare
  # the rerun from the stored fields (sections 6-7 before) against the real .fa run
  $PY plot_compare_outputs.py $RERUN/dpolvar_Arome_ICE3_Cband_20250831_1800.nc $C \
      --labels "rerun from output" "fa run" --vars Zh Zdr Kdp Ah --out $O/compare_rerun

  # ---- AROME W band from the .fa file
  O=$FIG/fa_W
  $PY plot_output_check.py $W --out $O/check
  $PY compute_cfad.py $W --cells-from $C --label "W band, C band cores" --mask core --vars Zh Zh_att Zdr Kdp \
      --zh-min -50 --out $O/stats/core_Ccells.nc
  $PY compute_cfad.py $C --label "C band" --mask core --vars Zh Zh_att Zdr Kdp --zh-min -50 --out $O/stats/C_core.nc
  $PY plot_cfad.py $O/stats/C_core.nc $O/stats/core_Ccells.nc --labels "C band" "W band" \
      --colors tab:blue tab:red --figs radar freq --diff --hmax 13 --name C_vs_W_core --out $O/cfad
  $PY plot_cross_section.py $W --points 2.3 45.1 2.8 45.9 --vars Zh Zh_att Zdr Kdp Ah --hmax 13 --out $O/sections
  $PY plot_compare_outputs.py $C $W --labels "C band" "W band" --vars Zh Zdr Kdp Ah Zh_att --zh-min -10 --out $O/compare
  $PY plot_compare_outputs.py $RERUN/dpolvar_Arome_ICE3_Wband_20250831_1800.nc $W \
      --labels "rerun from output" "fa run" --vars Zh Zdr Kdp Ah --out $O/compare_rerun
fi

if [ "$PART" = all ] || [ "$PART" = mesonh ]; then
  # ---- MesoNH, output of 20 Nov. on the whole domain (1201 x 1201 x 82)
  O=$FIG/mesonh_old
  $PY plot_output_check.py $MNH_OLD --out $O/check
  $PY compute_cfad.py $MNH_OLD --label "cores, whole domain (20 Nov.)" --mask core --out $O/stats/core.nc
  # ---- MesoNH recomputed with the current code on the subdomain of the cyclone
  O=$FIG/mesonh
  $PY plot_output_check.py $MNH --out $O/check
  $PY compute_cfad.py $MNH --label "cores (Zmax >= 40 dBZ)" --mask core --vars Zh Zh_att Zdr Kdp --out $O/stats/core.nc
  $PY compute_cfad.py $MNH --label "envelopes (Zmax >= 25 dBZ)" --mask envelope --out $O/stats/envelope.nc
  $PY compute_cfad.py $MNH --label "all points" --mask all --out $O/stats/all.nc
  $PY plot_cfad.py $O/stats/core.nc $O/stats/envelope.nc $O/stats/all.nc \
      --colors tab:red tab:orange tab:blue --figs summary --hmax 16 --name masks --out $O/cfad
  $PY plot_cfad.py $O/stats/core.nc --figs freq --hmax 16 --name core --out $O/cfad
  # west-east section through the maximum of Zh (56.96 E, 20.07 S)
  $PY plot_cross_section.py $MNH --points 55.5 -20.07 58.5 -20.07 --vars Zh Zh_att Zdr Kdp Ah --contents --hmax 16 \
      --out $O/sections
  $PY plot_compare_outputs.py $MNH_OLD $MNH --labels "20 Nov. run" "rerun 24 Nov." --vars Zh Zdr Kdp --out $O/compare
fi
echo "done $PART $(date '+%F %T')"
