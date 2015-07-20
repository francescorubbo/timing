#!/usr/bin/env python
from pylab import *
import numpy
from timingFunctions import *

pfilename='plots/tf_densitybased.pdf'
filename='Square_mu-150_psi-5_px-500_time-0.root'

SetupATLAS()

#truth,timefracs=timeFractionDensityBased(filename,rcut=0.1,cut=0.05,ptCentering=True)
truth,timefracs=timeFractionTruthBranch(filename,rcut=0.1,cut=0.05,ptCentering=True)
n,b,p = hist(timefracs[truth==1],bins=30,alpha=0.5)
hist(timefracs[truth==0],bins=b,alpha=0.5)

savefig(pfilename)
