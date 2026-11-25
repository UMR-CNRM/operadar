#!/bin/bash
# Generation of all the tables David2026PhD all bands
# Existing tables are kept (generate_tables.sh skips them).
#cd "$(dirname "$0")/../.."
bandlist="L S C X Ku K Ka W"
for b in $bandlist; do
    echo "$b"
done

# Run function
run() {
    echo "=== $(date '+%F %T') start $*"
    ./generate_tables.sh "$@" -j 8 > "tables_generator/tables/gen_${1#--}_$3.log" 2>&1;
    local st=$?
    echo "=== $(date '+%F %T') end $* (exit $st)";
    }

# Execution for each band
for b in $bandlist; do
    run --David2026PhD --band $b
done

echo "=== $(date '+%F %T') ALL DONE"
