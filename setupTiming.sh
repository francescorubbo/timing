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
    export LD_LIBRARY_PATH=${FASTJETLOCATION}lib/:$LD_LIBRARY_PATH
}

setup_lhapdf() {
    export LHAPDFLOCATION=/u/at/rubbo/nfs/Software/lhapdf-5.9.1/install/
    export LD_LIBRARY_PATH=${LHAPDFLOCATION}lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/usr/include
    export BOOSTLIBLOCATION=/usr/lib64
    export LD_LIBRARY_PATH=/u/at/rubbo/local/lib/:$LD_LIBRARY_PATH
}

setup_ROOT
setup_PYTHIA
setup_lhapdf
setup_fastjet
setup_boost
