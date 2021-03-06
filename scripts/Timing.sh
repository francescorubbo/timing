#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/u/at/rubbo/nfs/Software/pythia8183/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /u/at/rubbo/nfs/Software/root_v5.34.17/bin/thisroot.sh
}

setup_fastjet() {
	export FASTJETLOCATION=/u/at/rubbo/nfs/Software/fastjet-3.0.3/fastjet-install/
	#export FASTJETLOCATION=/u/at/pnef/Work/Code/TrackBasedGrooming/fastjet-3.0.3/fastjet-install/
    #export FASTJETLOCATION=/u/at/pnef/Work/Code/fastjet-install/
    export LD_LIBRARY_PATH=${FASTJETPATH}lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/usr/include
    export BOOSTLIBLOCATION=/usr/lib64
    export LD_LIBRARY_PATH=/u/at/rubbo/local/lib/:$LD_LIBRARY_PATH
}

env | grep -q "ROOT"
if [ $? -eq 1 ]
then
    setup_ROOT
fi

env | grep -q "PYTHIA"
if [ $? -eq 1 ]
then
    setup_PYTHIA
fi

env | grep -q "FAST"
if [ $? -eq 1 ]
then
    setup_fastjet
fi

env | grep -q "BOOST"
if [ $? -eq 1 ]
then
    setup_boost
fi

Timing $@