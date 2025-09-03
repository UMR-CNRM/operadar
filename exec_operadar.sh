#!/bin/bash

function show_help {
    echo "----------------------------------------------------------------"
    echo " This radar forward operator is developped at the CNRM, France. "
    echo "----------------------------------------------------------------"
    echo ""
    echo "Usage: $0 -f FILENAME -c CONFIG [--verbose]"
    echo ""
    echo "  -f FILENAME : Only the filename. Please use the config file to provide the path to access the file."
    echo "  -c CONFIG   : Before running the code, you need to create a configuration file in ./configFiles/ based on the template provided."
    echo "  --verbose   : Optional, show more details if activated (default value: False)"
    echo "  --append    : Optional, append the computed dual-polarimetric fields into the input file (default value: False)"
    echo "                /!\ If this option is activated, no netcdf output file will be saved. Only working for AROME .fa files"
    echo ""
}

if [ $# -lt 4 ]; then
    show_help
    exit 1
fi

FILENAME=""
CONFIG=""
VERBOSE=""
APPEND=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--filename)
            FILENAME="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE="--verbose"
            shift
            ;;
        --append)
            APPEND="--append"
            shift
            ;;
        *)
            show_help
            exit 1
            ;;
    esac
done

# Checking if FILENAME and CONFIG are provided
if [ -z "$FILENAME" ] || [ -z "$CONFIG" ]; then
    echo "/!\ Missing arguments /!\ "
    show_help
    exit 1
fi


python3 -u -m operadar $FILENAME $CONFIG $VERBOSE $APPEND
