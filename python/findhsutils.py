from sklearn.neighbors.kde import KernelDensity
from numpy import array,newaxis,linspace,exp,sqrt
import matplotlib.pyplot as plt

def gettime(times,width=0.05):
#    X = times[times>0.]
    X = times
    X_probe = linspace(min(X),max(X),200)
    X = X[:,newaxis]
    X_plot = X_probe[:,newaxis]
    
    kde = KernelDensity(kernel='gaussian', bandwidth=width).fit(X)
    log_dens = kde.score_samples(X_plot)

#    plt.fill(X_plot[:, 0], exp(log_dens))

    return X_probe[log_dens==max(log_dens)].mean()

def getcenter(var,detas,dphis,cut=0.1):

    drs = sqrt(detas*detas+dphis*dphis)
    return var[drs<cut]
