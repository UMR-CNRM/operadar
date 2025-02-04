#!/bin/bash

[ "$1" = "-h" -o "$1" = "--help" ] && echo "
Before running exec_operad.sh code, you need to create a configuration file in ./configFiles/ based on the template provided.

This script takes 3 arguments in this order :
  1 - the filename to read (+ additional tree structure if not provided in the configuration file)
  2 - date in the yyyymmdd format
  3 - name of the configuration file

----------
 Examples
----------
 * With input_file_dir  = "/cnrm/precip/users/davidcl/expeOLIVE/arome/3dvarfr/GO2V/20220818T1200P/forecast/" in config file :
   >>> ./exec_operad.sh historic.arome.franmg-01km30+0007:00.fa 20220818 conf_Arome_with_ICE3.py

 * With input_file_dir  = "/home/davidcl/Programmation/expeOLIVE/" in config file :
   >>> ./exec_operad.sh LIMA_simus_Arome/00Z/20220620/historic.arome.franmg-01km30+0003:00.fa 20220620 conf_Arome_with_LIMA_chaud.py
  
 * With input_file_dir  = "/home/cnrm_other/ge/mrmp/davidcl/MNH-V5-7-0/MY_RUN/" in config file :
   >>> ./exec_operad.sh LIMAH.1.SEG01.003.nc 20220818 conf_MesoNH_LIMA_CORSEbe.py
_____________________________________________________________________________________________________________________________________
" && exit


cp ./configFiles/$3 operadar/operad_conf.py

currentTime=`date +"%F_%R"`

python3 -u operadar/operad.py $1 $2