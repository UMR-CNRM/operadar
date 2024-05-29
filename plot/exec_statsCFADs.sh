#!/bin/bash
[ "$1" = "-h" -o "$1" = "--help" ] && echo "
  Before running exec_statsCFADs.sh code, you need to:
  * select the simulation and directory options in stats_for_CFADs.py
  _________________________________________
  Examples :
   >>> ./exec_statsCFADs.sh
  _________________________________________  
  STEPS :
     a) Creates or check the existance of the directory where all nohupFiles.txt will be saved (nohupOUT)
     b) Run stats_for_CFADs.py and save the prompt output into ./nohupOUT/nohupFile.txt
     c) Return in ./process_id_historic.txt the process id (for killing purpose)
  
  Note : process_id_historic.txt keep in memory all process id and the date/time the program was launched.
         You may want to manually erase some old lines or delete the file to clear process_id_historic.txt
" && exit

if [ ! -d ./nohupOUT ]; then
  mkdir nohupOUT
fi

currentTime=`date +"%F %R"`

# --- For tests -----
#python3 -i operad.py $1 $2 $3 $4

# --- For multiple times running ---
nohup python3 -u stats_for_CFADs.py >nohupOUT/out_statCFAD.txt 2> nohupOUT/err_statCFAD.txt &

echo ${currentTime} $! >> process_id_historic.txt

echo program output redirected to ./nohupOUT/out_statCFAD.txt
echo corresponding process id written in ./process_id_historic.txt
