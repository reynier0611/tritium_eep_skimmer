Welcome to skimmer.

We know it looks messy, but there is a lot of useful stuff.

    Executable programs:
	   skimmer	       skims 1trackL & 1trackR events from multiple
	   		           files and stores them in one
	   combiner            combines the output of multiple skims

    Shell scripts
	   run_skims.sh        runs the skimmer on all unskimmed-files
	   combine_all.sh      combines skim output for all settings
	   conbine_last.sh     combines skim output for the last shift
	   		           as well as all the slow kin data
				   up until the last shift

    Data directories
    	   infiles:            lists of files in each run setting
	   skim:               soft link to the directory on /chafs1
                                   where the skim files are stored
           combiner_out:       directory where combined rootfiles go


Novice-level instructions

   1) Update the run lists in ./infiles/ with the latest replayed runs
   2) Run  ./run_skims.sh
   3) Run  ./combine_all.sh
   4) Your combined rootfiles are now in ./combiner_out


More detailed instructions to follow
	... they never came
