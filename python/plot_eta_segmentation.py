#!/usr/bin/env python

from timingFunctions import *

SetupATLAS()

scale=1e-3
binsize=numpy.array([0.25,0.5,1.0,2.0,5.0])
binsize*=scale
etaMin=2.5
etaMax=4.3
for i in range(0,len(binsize)):
    sigma=binsize[i]
    x,y=etaBinsize(sigma,etaMin=etaMin,etaMax=etaMax)
    color=colors[i]
    plot(x,y,label="$"+str(sigma*1e3).format("%e")+"$ mm Pixels",color=colors[i])

xlim(etaMin,etaMax)
xlabel('$\eta$')
ylabel('$\Delta\eta$ of Tracker Pixel')
yscale('log')
ylim(3e-4,1e-1)
legend(loc='upper left',ncol=1,frameon=False)
savefig("plots/tracker_resolution")
