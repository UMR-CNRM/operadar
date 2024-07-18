#!/bin/bash
[ "$1" = "-h" -o "$1" = "--help" ] && echo "
  Before running compute_CFADs.sh code, you need to:
  * select the simulation and directory options in compute_CFADs.py
  _________________________________________
  Examples :
   >>> ./exec_computeCFAD.sh Arome_oper
   >>> ./exec_computeCFAD.sh MesoNH_ICE3
   >>> ./exec_computeCFAD.sh MesoNH_LIMA
   >>> ./exec_computeCFAD.sh MesoNH_LIMAAG
   >>> ./exec_computeCFAD.sh obs
   
  _________________________________________  
  STEPS :
     a) Creates or check the existance of the directory where all nohupFiles.txt will be saved (nohupOUT)
     b) Run compute_CFADs.py and save the prompt output into ./nohupOUT/nohupFile.txt
     c) Return in ./process_id_historic.txt the process id (for killing purpose)
  
  Note : process_id_historic.txt keep in memory all process id and the date/time the program was launched.
         You may want to manually erase some old lines or delete the file to clear process_id_historic.txt
" && exit

if [ ! -d ./nohupOUT ]; then
  mkdir nohupOUT
fi

currentTime=`date +"%F %R"`

# --- For tests -----
#python3 -i compute_CFADs.py $1


# --- For multiple times running ---
nohup python3 -u compute_CFADs.py $1 >nohupOUT/out_statCFAD_$1.txt 2> nohupOUT/err_statCFAD_$1.txt &

echo ${currentTime} $! >> process_id_historic_cfad.txt

echo program output redirected to ./nohupOUT/out_statCFAD_$1.txt
echo corresponding process id written in ./process_id_historic_cfad.txt
