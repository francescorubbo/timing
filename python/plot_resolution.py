#!/usr/bin/env python

import numpy
from pylab import *
from timingFunctions import *

SetupATLAS()

sigma=30.0
layers=[2,3,4,5,8,10]
effs=numpy.arange(0.5,1.01,0.01)
colors=['black','green','red','blue','cyan','magenta']
for i in range(0,len(layers)):
    color=colors[i]
    res=resolutionArray(sigma,layers[i],effs)
    plot(effs,res,label=str(layers[i]),color=color)
    idealres=resolution(sigma,layers[i],1)
    plot([0.5,1.0],[idealres,idealres],'--',color=color)

title("Tracker Timing Resolution v Layers")
xlim(0.5,1.0)
ylim(5,30)
ylabel('Timing Resolution [ps]')
xlabel('Hit Efficiency')
legend(loc='upper right',ncol=3,frameon=False,title="Layers")
savefig("plots/tracker_timing_resolution")
