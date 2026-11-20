#!/bin/bash
#====================================================================
#  Tmatrix tables generator 
# @authors: Cloé David and Clotilde Augros
#====================================================================

# Variables initialization
MODE=""
NEW_CONF=""
BAND=""

# Parameters list
HYDRO_LIST=("cs" "cl" "rr" "ii" "ss" "gg" "hh" "wg") #"wh" "ws"
BAND_LIST=("L" "S" "C" "X" "Ku" "K" "Ka" "W") #
MICRO_LIST=("ICE3") #("ICE3" "ICJW" "LIMA" "LIMC")

DEFAULT_MODES=("default" "David2025AMT" "David2026PhD" "vertical")


# Paths (relative to this script, so it can be launched from any directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARAM_FOLDER="${SCRIPT_DIR}/tables_generator/param"
TABLE_FOLDER="${SCRIPT_DIR}/tables_generator/tables"
SRC_DIR="${SCRIPT_DIR}/tables_generator/src"

# Validation functions
validate_value() {
    local -n array=$1
    [[ " ${array[*]} " =~ " ${2} " ]] && return 0
    return 1
}

valid_band() { validate_value BAND_LIST "$1"; }

# Help function
usage() {
    cat <<EOF
Usage:
  DEFAULT MODE:     $0 --default --band <value>
  Specialized MODE: $0 --David2025AMT/--David2026PhD/--vertical --band <value>
  NEWCONF MODE:     $0 --newConf <folder_name> --band <value>

Bands: ${BAND_LIST[@]}

DEFAULT and Specialized modes use param/TmatParam_<band><type>_<mode>.
NEWCONF mode copies param/TmatParam_<band><type> to TmatParam_<band><type>_<folder_name>
(edit the base files, or the copies, to change the configuration) and writes the
tables in tables/<folder_name>/.
EOF
    exit 1
}

# Argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --default|--David2025AMT|--David2026PhD|--vertical)
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
[[ -z "$MODE" || -z "$BAND" ]] && usage


# Core table generation function
FAILURES=0

generate_tables() {
    local output_subfolder="$1"
    local is_special_mode=false

    # Check if special mode
    if [[ " ${DEFAULT_MODES[@]} " =~ " ${output_subfolder} " ]]; then
        is_special_mode=true
        mode_type="PREDEFINED"
        mode_desc="with the predefined parameter files TmatParam_${BAND}<type>_${output_subfolder}"
    else
        mode_type="NEW CONFIGURATION"
        mode_desc="from the base parameter files TmatParam_${BAND}<type>"
    fi

    echo "====================================================================="
    echo "                ${mode_type} MODE - ${output_subfolder}              "
    echo "====================================================================="
    echo -e "Generating tables for all hydrometeor types ${mode_desc}"
    echo "Results will be stored in: ${TABLE_FOLDER}/${output_subfolder}/"
    echo ""
    echo "/!\ Table generation may take several hours..."

    mkdir -p "${TABLE_FOLDER}/${output_subfolder}"

    for H in "${HYDRO_LIST[@]}"; do
        echo -e "\n================ PROCESSING ${H} ============================="

        if [[ $is_special_mode == true ]]; then
            # Specialized modes: files in param/ with the mode name as extension
            PARAM_FILE="${PARAM_FOLDER}/TmatParam_${BAND}${H}_${output_subfolder}"
            if [[ ! -f "$PARAM_FILE" ]]; then
                echo "Error: Parameter file not found: $PARAM_FILE"
                FAILURES=$((FAILURES+1))
                continue
            fi
        else
            # New configuration: copy of the base file (without extension), kept for traceability
            local base_file="${PARAM_FOLDER}/TmatParam_${BAND}${H}"
            PARAM_FILE="${base_file}_${output_subfolder}"
            if [[ ! -f "$PARAM_FILE" ]]; then
                if [[ ! -f "$base_file" ]]; then
                    echo "Error: Base parameter file not found: $base_file"
                    FAILURES=$((FAILURES+1))
                    continue
                fi
                cp "$base_file" "$PARAM_FILE"
            fi
        fi

        # Scattering coefficients as a function of the diameter (Tmat), if not already there
        local coefdiff="${TABLE_FOLDER}/${output_subfolder}/TmatCoefDiff_${BAND}${H}"
        if [[ ! -f "$coefdiff" ]]; then
            cp "$PARAM_FILE" "${PARAM_FOLDER}/tmp_config"
            mkdir -p "${TABLE_FOLDER}/${H}"   # Tmat writes into tables/<type>/
            if ( cd "$SRC_DIR" && ./Tmat ) && [[ -f "${TABLE_FOLDER}/${H}/TmatCoefDiff_${BAND}${H}" ]]; then
                mv "${TABLE_FOLDER}/${H}/TmatCoefDiff_${BAND}${H}" "$coefdiff"
            else
                echo "Failed to compute the scattering coefficients for ${H} (Tmat)"
                FAILURES=$((FAILURES+1))
                continue
            fi
        else
            echo "Using existing table: $coefdiff"
        fi

        # Integration over the PSD for each microphysics scheme (TmatInt)
        for MICRO in "${MICRO_LIST[@]}"; do
            MOMENT=$(grep "^${MICRO}" "$PARAM_FILE" | cut -d: -f2)
            OUT_FILE="${TABLE_FOLDER}/${output_subfolder}/TmatCoefInt_${MICRO}_${MOMENT}_${BAND}${H}"

            if [[ ! -f "$OUT_FILE" ]]; then
                if "$SRC_DIR/TmatInt" "$SRC_DIR" "$output_subfolder" "$H" "$BAND" "$MICRO" "$MOMENT"; then
                    echo "Successfully generated tables for ${H} with ${MICRO} microphysics"
                else
                    echo "Failed to generate tables for ${H} with ${MICRO}"
                    FAILURES=$((FAILURES+1))
                fi
            else
                echo "Using existing table: $OUT_FILE"
            fi
        done
    done

    echo ""
    if [[ $FAILURES -gt 0 ]]; then
        echo "/!\ ${FAILURES} step(s) failed, see the messages above."
        return 1
    fi
    echo "All tables are in ${TABLE_FOLDER}/${output_subfolder}/"
}

# Main execution
case "$MODE" in
    --default|--David2025AMT|--David2026PhD|--vertical)
        generate_tables "${MODE#--}" || exit 1  # removes -- to get the name only
        ;;
    newConf)
        generate_tables "$NEW_CONF" || exit 1  # Personnalized config
        ;;
esac
