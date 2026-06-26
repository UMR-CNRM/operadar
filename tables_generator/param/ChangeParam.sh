#!/bin/bash

# Change param in TmatParam files

refband="C"
newband="L" #"S"
config="default" # improved
 
param1="band"
val1="C"
newval1="L" #"S"
param2="LAM"
val2="53.2"
newval2="190.4" #"106.2" "190.4"
for ficparam in `ls TmatParam_${refband}??_${config}`
do 
    echo $ficparam
    type=`echo -e $ficparam | tail -c11`
    newficparam=TmatParam_$newband$type
    cp $ficparam $newficparam
    echo "newficparam= "$newficparam
    echo "replace "$param1":"$val1 " by " $param1":"$newval1
    sed -i -e "s/$param1":"$val1/$param1":"$newval1/g" $newficparam
    sed -i -e "s/$param2":"$val2/$param2":"$newval2/g" $newficparam    
done
