#!/usr/bin/env python

import numpy
from math import *
from pylab import *
import rootpy
from rootpy.plotting.style import get_style, set_style

rootpy.log.basic_config_colorized()

#use latex for text
rc('text', usetex=True)
# set the style
style = get_style('ATLAS')
style.SetEndErrorSize(3)
set_style(style)

def y(eta):
    return 3.5/numpy.sinh(eta)

def eta_y(y):
    return numpy.arcsinh(3.5/y)

binsize=numpy.array([250.,500.0,1000.0])
binsize*=1e-6

print(y(2.5))
print(eta_y(y(2.5)))

etaMin=2.5
etaMax=4.3

xmin=y(etaMax)
xmax=y(etaMin)

for i in range(0,len(binsize)):
    sigma=binsize[i]
    print(sigma)
    x=numpy.arange(xmin,xmax-sigma,sigma)
    de=eta_y(x)-eta_y(x+sigma)
    xe=eta_y(x)
    plot(xe,de,label="$"+str(sigma*1e6).format("%e")+" \mu m$ Pixels")

title("$\eta$ Binsize of Tracker Pixels")
xlim(etaMin,etaMax)
xlabel('$\eta$')
ylabel('$\Delta\eta$ of Tracker Pixel')
legend(loc='upper left',ncol=1,frameon=False,title="Tracker Segmentation")
savefig("plots/tracker_resolution")
