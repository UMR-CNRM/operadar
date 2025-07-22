#!/bin/bash

# Variables initialization
MODE=""
NEW_CONF=""
BAND=""
HYDRO=""
ARF=""
ARV=""
CANTING=""
DSTY=""
RIMING=""
DIEL=""

# List of (fixed) parameters
HYDRO_LIST=("cs" "cl" "rr" "ii" "ss" "ws" "gg" "wg" "hh" "wh")
BAND_LIST=("C" "S" "X" "K" "W")
ARfunc_LIST=("AUds" "CNST" "BR02" "RYdg" "RYwg")
DSTYfunc_LIST=("BR07" "RHOX" "LS15" "ZA05")
DIELfunc_LIST=("Liebe91" "RY19dry" "LBwetgr" "MGwMA08")
MICRO_LIST=("ICE3" "LIMA")

# Errors storage
MISSING_FILES=()

# Paths
PARAM_FOLDER="./tables_generator/param"
TABLE_FOLDER="./tables_generator/tables"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMAT_DIR="$(dirname "$(realpath ./tables_generator/src/Tmat)")"
TMATINT_DIR="$(dirname "$(realpath ./tables_generator/src/TmatInt)")"


# Help function
usage() {
    echo " "
    echo "How to use :"
    echo "  DEFAULT MODE : $0 --default --band <val>"
    echo "  NEWCONF MODE : $0 --newConf <confName> --band <val>"
    echo "  EDIT MODE    : $0 --hydro <val> --band <val> [--arf <val>] [--arv <val>] [--canting <val>] [--dsty <val>] [--riming <val>] [--diel <val>]"
    echo " "
    echo "Accepted values for :"
    echo "  --band    (radar band)          : C, K, S, W, X"
    echo "  --hydro   (hydrometeor type)    : rr, ii, gg, ss, cc, wg, hh, wh "
    echo "  --arf     (axis ratio function) : AUds, CNST, BR02, RYdg, RYwg"
    echo "  --arv     (axis ratio value)    : any float value."
    echo "  --canting (canting angle)       : any float value."
    echo "  --dsty    (density function)    : BR07, RHOX, LS15, ZA05"
    echo "  --riming  (fraction of riming)  : any float value >= 1 (1=unrimed)"
    echo "  --diel    (dielectric function) : Liebe91, RY19dry, LBwetgr, MGwMA08"
    echo "Further details are available in the README.md "
    echo " "
    exit 1
}

# Checking validity
valid_band() {
    [[ "${BAND_LIST[*]}" =~ "$1" ]]
}
valid_hydro() {
    [[ "${HYDRO_LIST[*]}" =~ "$1" ]]
}
valid_arf(){
    [[ "${ARfunc_LIST[*]}" =~ "$1" ]]
}
valid_dstyf(){
    [[ "${DSTYfunc_LIST[*]}" =~ "$1" ]]
}
valid_dielf(){
    [[ "${DIELfunc_LIST[*]}" =~ "$1" ]]
}

# Reading arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --default)
            [[ -n "$MODE" ]] && { echo "/!\ Error: mutually exclusive options."; usage; exit 1; }
            MODE="default"
            ;;
        --newConf)
            [[ -n "$MODE" ]] && { echo "/!\ Error: mutually exclusive options."; usage; exit 1; }
            MODE="newConf"
            shift
            [[ -z "$1" || "$1" == --* ]] && { echo "/!\ Error: --newConf expect a folder name."; usage; exit 1; }
            NEW_CONF="$1"
            ;;
        --hydro)
            [[ -n "$MODE" && "$MODE" != "edit" ]] && { echo "/!\ Error: mutually exclusive options."; usage; exit 1; }
            MODE="edit"
            shift
            [[ -z "$1" || "$1" == --* ]] && { echo "/!\ Error: --hydro expect a value among : ${HYDRO_LIST[*]}"; exit 1; }
            valid_hydro "$1" || { echo "/!\ Error: wrong value for --hydro. Expect one of the following : ${HYDRO_LIST[*]}"; exit 1; }
            HYDRO="$1"
            ;;
        --arf)
            shift
            valid_arf "$1" || { echo "/!\ Error: invalid axis ratio function. Expect one of the following : ${ARfunc_LIST[*]}"; exit 1; }
            ARF="$1"
            ;;
        --arv) shift; ARV="$1" ;;
        --canting) shift; CANTING="$1" ;;
        --dsty)
            shift
            valid_dstyf "$1" || { echo "/!\ Error: invalid density function. Expect one of the following : ${DSTYfunc_LIST[*]}"; exit 1; }
            DSTY="$1"
            ;;
        --riming) shift; RIMING="$1" ;; # RAJOUTER TEST VALEUR <1 CAR =1 PAS DE RIMING, >1 RIMING
        --diel)
            shift
            valid_dielf "$1" || { echo "/!\ Error: invalid dielectric function. Expect one of the following : ${DIELfunc_LIST[*]}"; exit 1; }
            DIEL="$1"
            ;;
        --band)
            shift
            valid_band "$1" || { echo "/!\ Error: invalid band. Expect one of the following : ${BAND_LIST[*]}"; exit 1; }
            BAND="$1"
            ;;
        *)
            echo "/!\ Unknown option: $1"
            usage
            exit 1;
            ;;
    esac
    shift
done

# Checking arguments
if [[ -z "$MODE" ]]; then
    echo "/!\ Error: wrong use of the tables_generator. Can only be :"
    usage
fi

if [[ -z "$BAND" ]]; then
    echo "/!\ Error: --band argument is mandatory."
    usage
fi


generate_tables() {
    local output_subfolder="$1"

    if [[ "$output_subfolder" == "default" ]]; then
        echo "=========================================="
        echo "               DEFAULT MODE               "
        echo "=========================================="
        echo -e "Tables will be generated for all hydrometeor types with the values given by the tables_generator/param/TmatParam_*_default files."
        echo "Tables will be stored under ${TABLE_FOLDER}/${output_subfolder}/"
        echo "Progression of the table's generation for each hydrometeor is displayed under ./logs/{radarBand}_{hydrometeor}.log"
        echo "/!\ Table generation is time-consuming and can take several hours."
        echo ""
    else
        echo "=========================================="
        echo "               NewConf MODE               "
        echo "=========================================="
        echo -e "Tables will be generated for all hydrometeor types with the values given by the tables_generator/param/TmatParam_* files."
        echo "Tables will be stored under ${TABLE_FOLDER}/${output_subfolder}/"
        echo "Progression of the table's generation for each hydrometeor is displayed under ./logs/{radarBand}_{hydrometeor}.log"
        echo "/!\ Table generation is time-consuming and can take several hours."
        echo ""
    fi

    for H in "${HYDRO_LIST[@]}"; do

        echo -e "\n====== START OF THE PROGRAM FOR ${H} ======"
        
        if [[ "$output_subfolder" == "default" ]]; then
            PARAM_FILE="${PARAM_FOLDER}/TmatParam_${BAND}${H}_default"
        else
            PARAM_FILE="${PARAM_FOLDER}/TmatParam_${BAND}${H}"
            cp ${PARAM_FILE} "${PARAM_FILE}_${output_subfolder}"
        fi

        for MICRO in "${MICRO_LIST[@]}"; do
            MOMENT=$( cat "${PARAM_FILE}"| grep $MICRO |cut -d : -f2 )
            OUT_FILE="${TABLE_FOLDER}/${output_subfolder}/TmatCoefInt_${MICRO}_${MOMENT}_${BAND}${H}"

            if [[ -f "$OUT_FILE" ]]; then
                echo "$OUT_FILE already exist."
                
            else
                mkdir -p "${TABLE_FOLDER}/${output_subfolder}"
                mkdir -p "${TABLE_FOLDER}/${H}"
                if [[ -f "$PARAM_FILE" ]]; then
                    DIAMETER_TABLE="${TABLE_FOLDER}/${H}/TmatCoefDiff_${BAND}${H}"
                    if [[ ! -f "$DIAMETER_TABLE" ]]; then
                        cp "$PARAM_FILE" "${PARAM_FOLDER}/tmp_config"
                        echo "Launching the creation of the tables for a range of diameters."
                        # Temporary mv into Tmat directory to execute the f77 Tmat binary
                        pushd "$SCRIPT_DIR/tables_generator/src" > /dev/null
                        ./Tmat
                        if [[ $? -ne 0 ]]; then
                            echo "Error: Table creation failed for $H."                        
                        fi
                        popd > /dev/null
                        # Back to the Launching directory
                        if [ -f "$DIAMETER_TABLE" ]; then
                        cp "${DIAMETER_TABLE}" "${TABLE_FOLDER}/${output_subfolder}/TmatCoefDiff_${BAND}${H}"
                        fi               
                    else
                        echo "Table for the range of diameters already exists."
                    fi

                    INTEGRATED_TABLE="${TABLE_FOLDER}/${H}/TmatCoefInt_${MICRO}_${BAND}${H}"
                    echo "Integrating over the ${H} PSD for ${MICRO} microphysics (${MOMENT})"
                    if "$TMATINT_DIR/TmatInt" "$TMATINT_DIR" "$H" "$BAND" "$MICRO" "$MOMENT"; then
                        mv "$INTEGRATED_TABLE" "$OUT_FILE"
                        echo "Tables generated for $H with $MICRO microphysics."
                    else
                        echo "Error: Failed integration for $H with $MICRO microphysics."
                    fi
                else
                    echo "Missing or unknown file: $PARAM_FILE"
                fi
            fi
        done

        echo -e "\n====== END OF THE PROGRAM FOR ${H} ======"
                
    done

    echo " "
}


# Mode 1 : default
if [[ "$MODE" == "default" ]] ; then
    generate_tables "default"

# Mode 2 : newConf
elif [[ "$MODE" == "newConf" ]] ; then
    generate_tables "$NEW_CONF"
fi

# Mode 3 : edit/modification personnalisée
if [[ "$MODE" == "edit" ]]; then
    [[ -z "$HYDRO" ]] && { echo "Erreur : --hydro est requis en mode modification."; usage; }

    FILE="toto_${BAND}.txt"
    [[ ! -f "$FILE" ]] && { echo "Fichier $FILE introuvable."; exit 1; }

    echo "Génération de config.txt à partir de $FILE avec les paramètres modifiés..."

    cp "$FILE" config.txt

    # Fonction de remplacement dans le fichier
    update_param() {
        local key="$1"
        local value="$2"
        if grep -q "^$key=" config.txt; then
            sed -i "s/^$key=.*/$key=$value/" config.txt
        else
            echo "$key=$value" >> config.txt
        fi
    }

    update_param "hydro" "$HYDRO"
    [[ -n "$ARF" ]] && update_param "arf" "$ARF"
    [[ -n "$ARV" ]] && update_param "arv" "$ARV"
    [[ -n "$CANTING" ]] && update_param "canting" "$CANTING"
    [[ -n "$DSTY" ]] && update_param "dsty" "$DSTY"
    [[ -n "$RIMING" ]] && update_param "riming" "$RIMING"
    [[ -n "$DIEL" ]] && update_param "diel" "$DIEL"

    echo "config.txt généré avec modifications personnalisées."
    exit 0
fi
