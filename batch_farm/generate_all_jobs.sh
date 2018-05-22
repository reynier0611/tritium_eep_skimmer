#!/bin/bash

INFILEDIR=/work/halla/triton/Rey/skimmer/infiles/
for kin in mid fast slow fast2 mid2
do
    for A in 1H 2H 3H 3He dum
    do

        # Test if a runlist exists
        runlist=${INFILEDIR}${A}_${kin}_kin.dat

        if [ -e $runlist ]
        then
            ./genjob ${A}_${kin}_kin.dat
        fi

    done
done
