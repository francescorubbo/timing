import rootpy
from rootpy.plotting.style import get_style, set_style
from root_numpy import root2rec
from pylab import *
from findhsutils import *
from numpy import random
import numpy as np

filename = 'Square_mu-150_psi-5_px-500_time-0.root'
timeres = 0.01
leaves=['jtruth','jeta','jphi','j0cleta','j0clphi','j0cltime','j0cltruth']

array = root2rec(filename,'tree',leaves)
trutharrays=array['jtruth']

timediff = []

for evt in range(len(trutharrays)):

    jeta=array['jeta'][evt]
    if len(jeta)<1: continue
    jphi=array['jphi'][evt]
    times=array['j0cltime'][evt]
    etas=array['j0cleta'][evt]
    phis=array['j0clphi'][evt]
    truth=array['j0cltruth'][evt]

#    if timeres>0.:
#        times = times+random.normal(0.,timeres,len(times))

    times = getcenter(times,etas-jeta[0],phis-jphi[0],0.1)
    truth = getcenter(truth,etas-jeta[0],phis-jphi[0],0.1)

    times = times[(times>-990)]
    truth = truth[(times>-990)]

    #n,bins,patches = hist(times,bins=200)
    #hist(ths,bins=bins)

    if len(times[truth==0])>0:
        ths = times[truth==0][0]
        tmeas = gettime(times,timeres)
        timediff.append(ths-tmeas)
        
hist(timediff,bins=30)
print np.array(timediff).mean()
print np.array(timediff).std()

savefig('test.png')

