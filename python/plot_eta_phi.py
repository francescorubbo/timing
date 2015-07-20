#!/usr/bin/env python
from timingFunctions import *

SetupATLAS()

files=['Timing.root','Timing_withB.root']
labels=['no magnetic field','B=2T']
filename="event_display"

figure(figsize=(12,12))
event=1
for i in range(0,len(files)):
    subplot(2,2,i+1)
    plotEvent(files[i],event)
    title(labels[i])

savefig('plots/event_display_pileup.png')
