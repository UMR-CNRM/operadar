#!/bin/bash

# Change param in TmatParam files

newband="K"
param="LAMmin"
val="3.15"
newval="8.40"
for ficparam in `ls TmatParam_K??`
do 
    type=`echo -e $ficparam | tail -c3`
    newficparam=TmatParam_$newband$type
    cp $ficparam $newficparam
    echo "newficparam= "$newficparam
    echo "replace "$param":"$val " by " $param":"$newval
    sed -i -e "s/$param":"$val/$param":"$newval/g" $newficparam
    
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
