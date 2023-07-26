#!/bin/bash

[ "$1" = "-h" -o "$1" = "--help" ] && echo "
  Before running exec_operad.sh code, you need to modify/check common_settings.py
  This function takes 3 arguments in this order :
      1 - Arome or MesoNH
      2 - date into the yyyymmdd format
      3 - microphysics scheme name in capital letter 
  _________________________________________
  Example :
   >>> ./exec_operad.sh AROME 20220522 ICE4
  _________________________________________  
  STEPS :
     a) Creates or check the existance of the directory where all nohupFiles.txt will be saved (nohupOUT)
     b) Copy the right configuration file thanks to arguments 1 and 3
     c) Run operad.py and save the prompt output into ./nohupOUT/nohupFile.txt
     d) Return in ./process_id_historic.txt the process id (for killing purpose)
  
  Note : process_id_historic.txt keep in memory all process id + the function arguments and the date/time the program was launched.
         You may want to manually erase some old lines or delete the file to clear process_id_historic.txt
" && exit

if [ ! -d ./nohupOUT ]; then
  mkdir nohupOUT
fi

cp ./configFiles/conf_$1_$3.py operad_conf.py

currentTime=`date +"%D %R"`

nohup python3 operad.py $1 $2 $3 1>nohupOUT/output_$2_$3.txt 2>&1 &

echo ${currentTime} $! $1 $2 $3 >> process_id_historic.txt
