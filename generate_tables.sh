#!/bin/bash
#====================================================================
#  Tmatrix tables generator 
# @authors: Cloé David and Clotilde Augros
#====================================================================

# Variables initialization
MODE=""
NEW_CONF=""
BAND=""
JOBS=1   # number of Tmat processes run in parallel (option -j)

# Parameters list
HYDRO_LIST=("cs" "cl" "rr" "ii" "ss" "gg" "hh" "wg") #"wh" "ws"
BAND_LIST=("L" "S" "C" "X" "Ku" "K" "Ka" "W") #
# Types skipped for some bands (not needed for now). Remove an entry to compute those tables again.
declare -A SKIP_TYPES=(["Ku"]="hh wh" ["K"]="hh wh" ["Ka"]="hh wh" ["W"]="hh wh")
MICRO_LIST=("ICE3" "LIMA") #("ICE3" "ICJW" "LIMA" "LIMC")

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

Option:
  -j, --jobs <N>    run N Tmat and TmatInt processes in parallel (each computes a
                    block of temperatures), e.g. -j \$(nproc). Default: 1.

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
        -j|--jobs)
            shift
            [[ "$1" =~ ^[1-9][0-9]*$ ]] || { echo "/!\ Error: -j/--jobs expects a positive integer."; exit 1; }
            JOBS="$1"
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


# Run Tmat for one parameter file, split over $JOBS processes (blocks of temperatures).
# Each process writes its own files; they are then concatenated in temperature order
# (Tmat writes the header for the full range, so only the first block keeps it).
# Usage: run_tmat <param file> <output TmatCoefDiff> <output TmatResu> <log file>
run_tmat() {
    local param="$1" out="$2" resu="$3" log="$4"
    local ntc
    ntc=$(awk -F: '/^Tcmin:/{a=$2} /^Tcmax:/{b=$2} /^Tcstep:/{c=$2} END{printf "%d", (b-a)/c+1.5}' "$param")
    local njobs=$(( JOBS < ntc ? JOBS : ntc ))
    local chunk=$(( (ntc + njobs - 1) / njobs ))
    local tmp
    tmp=$(mktemp -d "${out}.tmp.XXXXXX")
    local pids=() j first last status=0
    njobs=$(( (ntc + chunk - 1) / chunk ))
    echo "Tmat: ${ntc} temperatures in ${njobs} process(es), log in ${log}"
    for (( j=0; j<njobs; j++ )); do
        first=$(( j * chunk ))
        last=$(( first + chunk - 1 ))
        (( last > ntc - 1 )) && last=$(( ntc - 1 ))
        ( cd "$SRC_DIR" && ./Tmat "$param" "$tmp/coef_$j" "$tmp/resu_$j" "$first" "$last" ) > "$tmp/log_$j" 2>&1 &
        pids+=($!)
    done
    for (( j=0; j<njobs; j++ )); do
        wait "${pids[$j]}" || status=1
    done
    : > "$log"
    for (( j=0; j<njobs; j++ )); do
        echo "===== Tmat block $j" >> "$log"
        cat "$tmp/log_$j" >> "$log"
    done
    if [[ $status -eq 0 ]]; then
        { cat "$tmp/coef_0"; for (( j=1; j<njobs; j++ )); do tail -n +4 "$tmp/coef_$j"; done; } > "$out"
        { cat "$tmp/resu_0"; for (( j=1; j<njobs; j++ )); do tail -n +2 "$tmp/resu_$j"; done; } > "$resu"
    else
        grep -h -iE "error|convergence|terminated" "$tmp"/log_* | sort -u | head -5
    fi
    rm -rf "$tmp"
    return $status
}

# Run TmatInt for one type and scheme, split over $JOBS processes (blocks of temperatures).
# Usage: run_tmatint <output subfolder> <type> <micro> <moment> <output TmatCoefInt>
run_tmatint() {
    local repo="$1" h="$2" micro="$3" moment="$4" out="$5"
    local coefdiff="${TABLE_FOLDER}/${repo}/TmatCoefDiff_${BAND}${h}"
    local ntc
    ntc=$(sed -n 2p "$coefdiff" | awk '{printf "%d", ($14-$12)/$13+1.5}')
    local njobs=$(( JOBS < ntc ? JOBS : ntc ))
    local chunk=$(( (ntc + njobs - 1) / njobs ))
    njobs=$(( (ntc + chunk - 1) / chunk ))
    local tmp
    tmp=$(mktemp -d "${out}.tmp.XXXXXX")
    local pids=() j first last status=0
    for (( j=0; j<njobs; j++ )); do
        first=$(( j * chunk ))
        last=$(( first + chunk - 1 ))
        (( last > ntc - 1 )) && last=$(( ntc - 1 ))
        "$SRC_DIR/TmatInt" "$SRC_DIR" "$repo" "$h" "$BAND" "$micro" "$moment" "$first" "$last" "$tmp/int_$j" \
            > "$tmp/log_$j" 2>&1 &
        pids+=($!)
    done
    for (( j=0; j<njobs; j++ )); do
        wait "${pids[$j]}" || status=1
    done
    if [[ $status -eq 0 ]]; then
        { cat "$tmp/int_0"; for (( j=1; j<njobs; j++ )); do tail -n +4 "$tmp/int_$j"; done; } > "$out"
    else
        cat "$tmp"/log_* | grep -iE "error|not found" | sort -u | head -5
    fi
    rm -rf "$tmp"
    return $status
}

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
        if [[ " ${SKIP_TYPES[$BAND]:-} " == *" $H "* ]]; then
            echo -e "\n================ SKIPPING ${H} (not computed for ${BAND} band, see SKIP_TYPES) ==="
            continue
        fi
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
            if run_tmat "$PARAM_FILE" "$coefdiff.part" \
                        "${TABLE_FOLDER}/${output_subfolder}/TmatResu_${BAND}${H}" \
                        "${TABLE_FOLDER}/${output_subfolder}/Tmat_${BAND}${H}.log"; then
                mv "$coefdiff.part" "$coefdiff"
            else
                rm -f "$coefdiff.part"
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
                if run_tmatint "$output_subfolder" "$H" "$MICRO" "$MOMENT" "$OUT_FILE"; then
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
