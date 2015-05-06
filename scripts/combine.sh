#!/bin/bash

psis="5"
mus="140 200"
pixelSizes="0 200"
profiles="0 1"
timeModes="0 1"

path=files/*/

for mu in $mus ; do
    for psi in $psis ; do
	for px in $pixelSizes ; do
	    for tm in $timeModes ; do
		hadd -f -k ./Gauss_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root ${path}Sample_mu-${mu}_psi-${psi}_nevents-200_prof-0_px-${px}_time-${tm}_job-*
		if [ $? -gt 0 ]
		then
		    exit
		fi
		hadd -f -k ./Square_mu-${mu}_psi-${psi}_px-${px}_time-${tm}.root ${path}Sample_mu-${mu}_psi-${psi}_nevents-200_prof-1_px-${px}_time-${tm}_job-*
		if [ $? -gt 0 ]
		then
                    exit
                fi
	    done
	done
    done
done