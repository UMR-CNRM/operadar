#!/bin/bash

[ "$1" = "-h" -o "$1" = "--help" ] && echo "
  Before running exec_operad.sh code, you need to:
  * create a case study file in ./study_cases/
  * create a configuration file in ./configFiles/
  
  This script takes 4 arguments in this order :
      1 - Arome or MesoNH
      2 - date into the yyyymmdd format or "all"
      3 - microphysics scheme name in capital letter (ICE3, ICE4, LIMA, LIMAAG)
      4 - Config file specifying directories and forward operator options 
  _________________________________________
  Examples :
   >>> ./exec_operad.sh MesoNH 20220818 ICE3 conf_MesoNH_ICE3_CORSEbe.py
   >>> ./exec_operad.sh MesoNH 20220818 LIMA conf_MesoNH_LIMA_CORSEbe.py
   >>> ./exec_operad.sh Arome 20220818 ICE3 conf_Arome_ICE3_CORSEbe.py
   >>> ./exec_operad.sh MesoNH 20220818 LIMAAG conf_MesoNH_LIMAAG_CORSEbe.py
  _________________________________________  
  STEPS :
     a) Creates or check the existance of the directory where all nohupFiles.txt will be saved (nohupOUT)
     b) Copy the right configuration file thanks to arguments 4
     c) Run operad.py and save the prompt output into ./nohupOUT/nohupFile.txt
     d) Return in ./process_id_historic.txt the process id (for killing purpose)
  
  Note : process_id_historic.txt keep in memory all process id + the function arguments and the date/time the program was launched.
         You may want to manually erase some old lines or delete the file to clear process_id_historic.txt
" && exit

if [ ! -d ./nohupOUT ]; then
  mkdir nohupOUT
fi

cp ./configFiles/$4 operad_conf.py

currentTime=`date +"%F %R"`

# --- For tests -----
#python3 -i operad.py $1 $2 $3

# --- For multiple times running ---
nohup python3 -u operad.py $1 $2 $3 >nohupOUT/out_$1_$2_$3.txt 2> nohupOUT/err_$1_$2_$3.txt &

echo ${currentTime} $! $1 $2 $3 >> process_id_historic.txt

echo program output redirected to ./nohupOUT/out_$1_$2_$3.txt
echo corresponding process id written in ./process_id_historic.txt
