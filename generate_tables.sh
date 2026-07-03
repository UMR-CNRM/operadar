#!/bin/bash
#====================================================================
#  Tmatrix tables generator 
# @authors: Cloé David and Clotilde Augros
#====================================================================

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

# Parameters list
HYDRO_LIST=("cs" "cl" "rr" "ii" "ss" "gg" "hh" "wg") #"wh" "ws"
BAND_LIST=('K') #("L" "S" "C" "X" "Ku" "K" "Ka" "W") #
ARfunc_LIST=("AUds" "CNST" "BR02" "RYdg" "RYwg")
DSTYfunc_LIST=("BR07" "RHOX" "LS15" "ZA05" "TH08")
DIELfunc_LIST=("Liebe91" "RY19dry" "LBwetgr" "MGwMA08")
MICRO_LIST=("THOM") #("ICE3" "ICJW" "LIMA" "LIMC" "THOM")

DEFAULT_MODES=("default" "David2025AMT" "David2026PhD" "vertical" "thomfull" "thompsd" "thomdens" "thomfall" "thomtest")

# Errors storage
MISSING_FILES=()

# Paths
PARAM_FOLDER="./tables_generator/param"
TABLE_FOLDER="./tables_generator/tables"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMAT_DIR="$(dirname "$(realpath ./tables_generator/src/Tmat)")"
TMATINT_DIR="$(dirname "$(realpath ./tables_generator/src/TmatInt)")"

# Validation functions
validate_value() {
    local -n array=$1
    [[ " ${array[*]} " =~ " ${2} " ]] && return 0
    return 1
}

valid_band() { validate_value BAND_LIST "$1"; }
valid_hydro() { validate_value HYDRO_LIST "$1"; }
valid_arf() { validate_value ARF_LIST "$1"; }
valid_dstyf() { validate_value DSTYfunc_LIST "$1"; }
valid_dielf() { validate_value DIELfunc_LIST "$1"; }

# Help function
usage() {
    cat <<EOF
Usage:
  DEFAULT MODE:     $0 --default --band <value>
  Specialized MODE: $0 --David2025AMT/--David2026PhD/--vertical/--thomfull/--thompsd/--thomdens/--thomfall/--thomtest --band <value>
  NEWCONF MODE:     $0 --newConf <folder_name> --band <value>
  EDIT MODE:        $0 --hydro <value> --band <value> [options]

Options:
  --arf:     ${ARF_LIST[@]}
  --arv:     Any float
  --canting: Any float
  --dsty:    ${DSTYfunc_LIST[@]}
  --riming:  Any float >= 1
  --diel:    ${DIELfunc_LIST[@]}
EOF
    exit 1
}

# Argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --default|--David2025AMT|--David2026PhD|--vertical|--thomfull|--thompsd|--thomdens|--thomfall|--thomtest)
            [[ -n "$MODE" ]] && { echo "/!\ Error: mutually exclusive options."; usage; exit 1; }
            MODE="$1"
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

# Validation checks
[[ -z "$MODE" || -z "$BAND" || "$MODE" == "edit" && -z "$HYDRO" ]] && usage

if [[ -z "$BAND" ]]; then
    echo "/!\ Error: --band argument is mandatory."
    usage
fi

# Core table generation function
generate_tables() {
    local output_subfolder="$1"
    local is_special_mode=false

    # Check if special mode
    if [[ " ${DEFAULT_MODES[@]} " =~ " ${output_subfolder} " ]]; then
        is_special_mode=true
    fi

    echo "====================================================================="
    echo "                ${mode_type} MODE - ${output_subfolder}              "
    echo "====================================================================="
    echo -e "Generating tables for all hydrometeor types ${mode_desc}"
    echo "Results will be stored in: ${TABLE_FOLDER}/${output_subfolder}/"
    echo ""
    echo "/!\ Table generation may take several hours..."

    # Création des répertoires
    mkdir -p "${TABLE_FOLDER}/${output_subfolder}"

    for H in "${HYDRO_LIST[@]}"; do
        echo -e "\n================ PROCESSING ${H} ============================="

        # Handle specialized modes (files in param/ with extensions)
        if [[ $is_special_mode == true ]]; then
            PARAM_FILE="${PARAM_FOLDER}/TmatParam_${BAND}${H}_${output_subfolder}"

            if [[ ! -f "$PARAM_FILE" ]]; then
                echo "Error: Parameter file not found: $PARAM_FILE"
                continue
            fi

            # Generate diameter tables if needed
            if [[ ! -f "${TABLE_FOLDER}/${output_subfolder}/TmatCoefDiff_${BAND}${H}" ]]; then
                cp "$PARAM_FILE" "${PARAM_FOLDER}/tmp_config"
                pushd "$SCRIPT_DIR/tables_generator/src" > /dev/null
                ./Tmat
                popd > /dev/null

                if [[ -f "${TABLE_FOLDER}/${H}/TmatCoefDiff_${BAND}${H}" ]]; then
                    mv "${TABLE_FOLDER}/${H}/TmatCoefDiff_${BAND}${H}" \
                       "${TABLE_FOLDER}/${output_subfolder}/TmatCoefDiff_${BAND}${H}"
                fi
            fi
        
        # Handle custom modes (files in param/ without extension)
        else
            PARAM_FILE="${PARAM_FOLDER}/TmatParam_${BAND}${H}"
            OUTPUT_PARAM_FILE="${PARAM_FILE}_${output_subfolder}"

            if [[ ! -f "$PARAM_FILE" ]]; then
                echo "Error: Base parameter file not found: $PARAM_FILE"
                continue
            fi

            # Create copy with extension for this custom config
            if [[ ! -f "$OUTPUT_PARAM_FILE" ]]; then
                cp "$PARAM_FILE" "$OUTPUT_PARAM_FILE"
            fi

            PARAM_FILE="$OUTPUT_PARAM_FILE"
        fi

        # Generate integration tables for all microphysics
        for MICRO in "${MICRO_LIST[@]}"; do
            MOMENT=$(grep "^${MICRO}" "$PARAM_FILE" | cut -d: -f2)
            OUT_FILE="${TABLE_FOLDER}/${output_subfolder}/TmatCoefInt_${MICRO}_${MOMENT}_${BAND}${H}"

            if [[ ! -f "$OUT_FILE" ]]; then
                if "$TMATINT_DIR/TmatInt" "$TMATINT_DIR" "$output_subfolder" "$H" "$BAND" "$MICRO" "$MOMENT"; then
                    echo "Successfully generated tables for ${H} with ${MICRO} microphysics"
                else
                    echo "Failed to generate tables for ${H} with ${MICRO}"
                fi
            else
                echo "Using existing table: $OUT_FILE"
            fi
        done

        echo -e "\n================== END OF THE PROGRAM FOR ${H} =================="
                
    done

    echo " "
}

# Main execution
case "$MODE" in
    --default|--David2025AMT|--David2026PhD|--vertical|--thomfull|--thompsd|--thomdens|--thomfall|--thomtest)
        generate_tables "${MODE#--}"  # removes -- to get the name only
        ;;
    --newConf)
        generate_tables "$NEW_CONF"  # Personnalized config
        ;;
    *)
        # Edit mode - create custom config
        [[ -z "$HYDRO" ]] && { echo "Error: --hydro is required in edit mode"; usage; }

        if [[ ! -f "toto_${BAND}.txt" ]]; then
            echo "Error: Configuration file 'toto_${BAND}.txt' not found"
            exit 1
        fi

        echo "Creating custom config.txt..."
        cp "toto_${BAND}.txt" "config.txt"

        # Update parameters
        local config_file="config.txt"
        {
            echo "[EDITED]"
            echo "hydro=${HYDRO}"
            [[ -n "$ARF" ]] && echo "arf=${ARF}"
            [[ -n "$ARV" ]] && echo "arv=${ARV}"
            [[ -n "$CANTING" ]] && echo "canting=${CANTING}"
            [[ -n "$DSTY" ]] && echo "dsty=${DSTY}"
            [[ -n "$RIMING" ]] && echo "riming=${RIMING}"
            [[ -n "$DIEL" ]] && echo "diel=${DIEL}"
        } | sed -i '/^[[:space:]]*$/d' "$config_file"

        echo "config.txt created successfully with custom parameters"
        exit 0
        ;;
esac
