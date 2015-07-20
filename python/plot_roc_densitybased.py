#!/usr/bin/env python
from pylab import *
import numpy
from timingFunctions import *

pfilename='plots/roc_densitybased.pdf'
files=['Square_mu-150_psi-5_px-500_time-0.root']
colors=['black','blue','red','green','orange','magenta','cyan']
smears=[0.01,0.02,0.03]

SetupATLAS()
clf() 
plot([0,1],[0,1],':',label="No Discrimination",color='black')

filename=files[0]
cuts=[0.1]
rcuts=[0.1]
linestyles=['-','--',':']
for i in range(0,len(rcuts)):
    for j in range(0,len(cuts)):
        print(filename)
        color=colors[j]
        rcut=rcuts[i]
        cut=cuts[j]
        truth,times=timeFractionDensityBased(filename,rcut=rcut,cut=cut)
        plotROC(truth,times,color=color,linestyle=linestyles[i])

for i in range(0,len(rcuts)): 
    plot([],[],linestyles[i],color='grey',label="$R_{cut}="+str(rcuts[i])+"$")
for i in range(0,len(cuts)): 
    plot([],[],'-',color=colors[i],label="$c="+str(cuts[i])+"$")

xlabel('Efficiency')
ylabel('Fake Rate')
xlim(0,1)
ylim(0,1)
legend(loc='upper left',frameon=False)
savefig(pfilename)
