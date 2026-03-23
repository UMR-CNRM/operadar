#!/bin/bash

# Change param in TmatParam files

newband="L" #"S"
param1="band"
val1="C"
newval1="L" #"S"
param2="LAM"
val2="53.2"
newval2="190.4" #"106.2" "190.4"
for ficparam in `ls TmatParam_C??`
do 
    type=`echo -e $ficparam | tail -c3`
    newficparam=TmatParam_$newband$type
    cp $ficparam $newficparam
    echo "newficparam= "$newficparam
    echo "replace "$param1":"$val1 " by " $param1":"$newval1
    sed -i -e "s/$param1":"$val1/$param1":"$newval1/g" $newficparam
    sed -i -e "s/$param2":"$val2/$param2":"$newval2/g" $newficparam    
done


#band=X
# param="expDstep"
# val="0.02"
# newval="0.005"
# for ficparam in `ls TmatParam_?w? TmatParam_?hh`
# do
#     type=`echo -e $ficparam | tail -c3`
#     echo "ficparam= "$ficparam
#     echo "replace "$param":"$val " by " $param":"$newval
#     sed -i -e "s/$param":"$val/$param":"$newval/g" $ficparam
    
# done
