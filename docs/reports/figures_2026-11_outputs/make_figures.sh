#!/bin/bash
# Figures of the report 2026-11_operadar_output_analysis.md, made with the plot tools of plot_tools/.
# Run from the repository root, in the operadar Python environment:
#     bash docs/reports/figures_2026-11_outputs/make_figures.sh [PYTHON]
set -e
PY=${1:-python}
ROOT=$(pwd)
OUT=$ROOT/docs/reports/figures_2026-11_outputs
F=$ROOT/operadarFiles/AROME/dpolvar_Arome_ICE3_Cband_20250831_1800.nc        # run of 20 Nov. 2026
F_OLD=$ROOT/modelFiles/AROME/ncfiles/dpolvar_Arome_ICE3_Cband_20250831_1800.nc # run of 20 Apr. 2026
cd $ROOT/plot_tools

# 1. Quick look and checks
$PY plot_output_check.py $F --out $OUT/check

# 2. CFADs: cores, envelopes, all points
$PY compute_cfad.py $F --label "cores (Zmax >= 40 dBZ)" --mask core --out $OUT/stats/core.nc
$PY compute_cfad.py $F --label "envelopes (Zmax >= 25 dBZ)" --mask envelope --out $OUT/stats/envelope.nc
$PY compute_cfad.py $F --label "all points" --mask all --out $OUT/stats/all.nc
$PY plot_cfad.py $OUT/stats/core.nc $OUT/stats/envelope.nc $OUT/stats/all.nc \
    --colors tab:red tab:orange tab:blue --figs summary --hmax 13 --name masks --out $OUT/cfad
$PY plot_cfad.py $OUT/stats/core.nc --figs freq --hmax 13 --name core --out $OUT/cfad
$PY plot_cfad.py $OUT/stats/core.nc --figs contents --contents-log --hmax 13 --name core_log --out $OUT/cfad

# 3. Cross-sections
$PY plot_cross_section.py $F --hmax 13 --out $OUT/sections
$PY plot_cross_section.py $F --points 2.3 45.1 2.8 45.9 --vars Zh Zdr Kdp Ah --contents --hmax 13 --out $OUT/sections

# 4. Comparison with the run of April 2026 (same model file)
$PY plot_compare_outputs.py $F_OLD $F --labels "April 2026" "November 2026" --out $OUT/compare

# 5. HTML version of the report (figures embedded; needs the markdown package)
cd $ROOT
$PY docs/reports/figures_2026-11_outputs/make_html.py docs/reports/2026-11_operadar_output_analysis.md
