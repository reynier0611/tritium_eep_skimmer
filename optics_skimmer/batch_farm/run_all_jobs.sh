#!/bin/bash

for kin in fast mid slow fast2 mid2
do
    for A in Opt Csing
    do
	swif run -workflow optics_skim_eep_${A}_${kin}
    done
done
