#!/bin/bash

INFILEDIR=/work/halla/triton/Rey/skimmer/optics_skimmer/infiles/
for kin in mid fast slow fast2 mid2
do
    for A in Opt Csing
    do

        # Test if a runlist exists
        runlist=${INFILEDIR}${A}_${kin}_kin.dat

        if [ -e $runlist ]
        then
            ./genjob ${A}_${kin}_kin.dat
        fi

    done
done
