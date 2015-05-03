#!/bin/bash

psis="5"
mus="80"
pixelSizes="0 50"
profiles="0 1"
timeModes="0 1"

for mu in $mus ; do
    for psi in $psis ; do
	for px in $pixelSizes ; do
	    for tm in $timeModes ; do
		hadd -f Gauss_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root Sample_mu-${mu}_psi-${psi}_nevents-200_prof-0_px-${px}_time-${tm}_job-*
		if [ $? -gt 0 ]
		then
		    exit
		fi
		hadd -f Square_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root Sample_mu-${mu}_psi-${psi}_nevents-200_prof-1_px-${px}_time-${tm}_job-*
		if [ $? -gt 0 ]
		then
                    exit
                fi
	    done
	done
    done
done