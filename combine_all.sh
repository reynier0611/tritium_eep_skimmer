#!/bin/bash

# Combines all the skimmed files for a given kinematic and nucleus into one root file

for A in 3H 3He 1H 2H dum
do
    for kin in mid fast slow fast2 mid2
    do 
	if [ -e skim/${A}_${kin} ]
	then

	    if [ `ls skim/${A}_${kin} | wc -l` -eq "0" ]
	    then
		echo "$A $kin has no runs to combine."
	    else 
		echo "Working on $A $kin ..."
		./combiner skim/${A}_${kin}/*.root
		mv combine_out.root combiner_out/${A}_${kin}.root
	    fi
	fi
    done
done
