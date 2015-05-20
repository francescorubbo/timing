#!/bin/bash

psis="5"
mus="10"
pixelSizes="500"
profiles="1"
timeModes="0"

path=files/*/

for mu in $mus ; do
    for psi in $psis ; do
	for px in $pixelSizes ; do
	    for tm in $timeModes ; do
		for prof in $profiles ; do
		    if [ $prof -eq 0 ]
		    then
			hadd -f -k ./Gauss_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root ${path}Sample_mu-${mu}_psi-${psi}_nevents-*_prof-0_px-${px}_time-${tm}_job-*
			if [ $? -gt 0 ]
			then
			    exit
			fi
		    fi
		    
		    if [ $prof -eq 1 ]
		    then
			hadd -f -k ./Square_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root ${path}Sample_mu-${mu}_psi-${psi}_nevents-*_prof-1_px-${px}_time-${tm}_job-*
			if [ $? -gt 0 ]
			then
			    exit
			fi
		    fi
		done
	    done
	done
    done
done