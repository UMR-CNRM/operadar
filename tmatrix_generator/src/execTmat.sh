#!/bin/bash

function show_help {
    echo "----------------------------------------"
    echo " Executable to create the lookup tables "
    echo "----------------------------------------"
    echo ""
    echo "Usage: $0 -hydro HYDRO -af ARfunc -av ARvalue -c CANTING -dsty DSTYfunc -riming RIMING -diel DIELfunc"
    echo ""
    echo "  -hydro HYDRO   : rr, ii, gg, ss, tt, wg, hh, wh "
    echo "  -af ARfunc     : AUds, CNST, BR02, RYdg, RYwg"
    echo "  -av ARvalue    : any float value."
    echo "  -c CANTING     : any float value."
    echo "  -dsty DSTYfunc : BR07, RHOX"
    echo "  -riming RIMING : any float value >= 1 (1=unrimed)"
    echo "  -diel DIELfunc : Liebe91, RY19dry, LBwetgr, MGwMA08"
    echo ""
}

if [ $# -lt 4 ]; then
    show_help
    exit 1
fi

HYDRO=$1
ARfunc=$2
ARvalue=$3
CANTING=$4
DSTYfunc=$5
RIMING=$6
DIELfunc=$7
BAND="C"

tables_folder="../tables/"$HYDRO"/"$ARfunc"_AR"$(printf "%.1f" $ARvalue)\
"_CANT"$(printf "%.1f" $CANTING)"_"$DSTYfunc"_Frim"$(printf "%.1f" $RIMING)"_"$DIELfunc"/"

mkdir -p $tables_folder

tables_param="../param/TmatParam_"$BAND$HYDRO

if [ $CANTING == 0 ]; then
    sed -i "s/^canting:.*/canting:nocant/" $tables_param
else
    sed -i "s/^canting:.*/canting:analyt/" $tables_param
fi

sed -i "s/^ARfunc:.*/ARfunc:"$ARfunc"/" $tables_param
sed -i "s/^ARcnst:.*/ARcnst:"$(printf "%.1f" $ARvalue)"/" $tables_param
sed -i "s/^sigbeta:.*/sigbeta:"$(printf "%.1f" $CANTING)"/" $tables_param
sed -i "s/^DSTYfunc:.*/DSTYfunc:"$DSTYfunc"/" $tables_param
sed -i "s/^DIEL:.*/DIEL:"$DIELfunc"/" $tables_param
sed -i "s/^Frim:.*/Frim:"$(printf "%.1f" $RIMING)"/" $tables_param

Tmat

mv "../tables/"$HYDRO"/TmatCoefDiff_"$BAND$HYDRO $tables_folder"/TmatCoefDiff_"$BAND$HYDRO
mv "../tables/"$HYDRO"/TmatResu_"$BAND$HYDRO $tables_folder"/TmatResu_"$BAND$HYDRO

