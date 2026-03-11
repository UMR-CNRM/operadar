#!/bin/bash

# Change param in TmatParam files

newband="L"
param="LAM"
val="53.2"
newval="190.4"
for ficparam in `ls TmatParam_C??`
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
