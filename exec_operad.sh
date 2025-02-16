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
}

if [ $# -lt 4 ]; then
    show_help
    exit 1
fi

FILENAME=""
CONFIG=""
VERBOSE=""

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


cp ./configFiles/$CONFIG operadar/operadar_conf.py

python3 -u operadar/operadar.py $FILENAME $VERBOSE