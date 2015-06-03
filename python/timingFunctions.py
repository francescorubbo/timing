import rootpy
from rootpy.plotting.style import get_style, set_style
from root_numpy import root2rec
from pylab import *
import numpy
from numpy import sin,cos,arctan2,sinh,arcsinh,pi
from itertools import chain
from matplotlib.patches import Circle

def SetupATLAS():
    rootpy.log.basic_config_colorized()
    #use latex for text
    rc('text', usetex=True)
    # set the style
    style = get_style('ATLAS')
    style.SetEndErrorSize(3)
    set_style(style)

def htype(x):
    if(x > 0.5):
        return 1
    else:
        return 0

def etaPhi(x,y,z):
    phi=arctan2(y,x)
    if(phi < 0):
        phi=phi+2*pi
    eta=arcsinh(z*cos(phi)/x)
    return eta,phi

def xy(eta,phi,z):
    h=z/sinh(eta)
    return h*cos(phi),h*sin(phi)

def pixelNum(filename):
    leaves=['j0clpixelNum']
    array = root2rec(filename,'tree',leaves)
    pnum=array['j0clpixelNum']

    pchain = chain.from_iterable(pnum)
    pnums=numpy.array(list(pchain))
    gpts=numpy.where(pnums > -1)
    return pnums[gpts]

ParticleMass={11:0.51,
              22:0.0,
              130:497.0,
              211:140.0,
              321:494.0,
              2112:940.0,
              2212:938.0}

ParticleName={11:'$e^{\pm}$',
              22:'$\gamma$',
              130:'$K_{L0}$',
              211:'$\pi^{\pm}$',
              321:'$K^{\pm}$',
              2112:'$N(\overline{N})$',
              2212:'$p^{\pm}$'}

def Gamma(pdgid,pt,eta):
    mass=ParticleMass[abs(pdgid)]/1.0e3
    if(mass != 0):
        ptot2=pow(pt*cosh(eta),2)
        mass2=pow(mass,2)
        return sqrt(1+ptot2/mass2)
    else:
        return 1.0e3

def GammaPtPzID(filename,truthOnly=False):
    leaves=['j0clpdgid','j0clpt','j0cleta','j0cltruth']
    array = root2rec(filename,'tree',leaves)
    ids=array['j0clpdgid']
    etas=array['j0cleta']
    pts=array['j0clpt']
    truth=array['j0cltruth']

    gammas=list()
    pTout=list()
    pz=list()
    idsOut=list()
    limit=min([len(ids),100])
    for i in range(0,limit):
        for j in range(0,len(ids[i])):
            if((not truthOnly) or (truth[i][j] == 0)):
                gammas.append(Gamma(ids[i][j],pts[i][j],etas[i][j]))
                pTout.append(pts[i][j])
                pz.append(pts[i][j]*sinh(etas[i][j]))
                idsOut.append(abs(ids[i][j]))
    return numpy.array(gammas),numpy.array(pTout),numpy.array(pz),numpy.array(idsOut)

def jetMultiplicity(filename):
    leaves=['jtruth','jpt','truejpt','jeta','jphi','truejeta','truejphi']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']

    truthnum=list()
    pileupnum=list()
    for i in range(0,len(trutharrays)):
        truth=trutharrays[i]
        tnum=(truth > 0.5).sum()
        if(tnum > 0):
            truthnum.append(tnum)
        pnum=(truth <= 0.5).sum()
        if(pnum > 0):
            pileupnum.append(pnum)

    return numpy.array(truthnum),numpy.array(pileupnum)

def ptCorrection(filename,ptMin=0.0,ptMax=1000.0):
    leaves=['jtruth','jpt','truejpt','jeta','jphi','truejeta','truejphi']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']
    allpt=array['jpt']
    alleta=array['jeta']
    allphi=array['jphi']
    trueeta=array['truejeta']
    truephi=array['truejphi']
    truept=array['truejpt']

    correction=list()
    for i in range(0,len(trutharrays)):
        if((len(trutharrays[i]) > 0) and (len(truept[i]) > 0)):
            frac=trutharrays[i][0]
            if(frac > 0.5):
                eta=alleta[i][0]
                phi=allphi[i][0]
                for j in range(0,len(trueeta[i])):
                    teta=trueeta[i][j]
                    tphi=truephi[i][j]
                    dist=sqrt(pow(eta-teta,2)+pow(phi-tphi,2))
                    if(dist <= 0.4):
                        pt=allpt[i][0]
                        if((pt > ptMin) and (pt < ptMax)):
                            correction.append(allpt[i][0]-truept[i][j])
                        j=len(trueeta[i])

    print(len(correction))
    return numpy.array(correction)

def ptFractionTruthBranch(filename,keepAll=False,minEta=0,maxEta=10,minPt=0,maxPt=100000,absTime=True):
    leaves=['jtruth','jtime','jeta','jpt']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']
    alltimes=array['jtime']
    alleta=array['jeta']
    allpt=array['jpt']

    times=numpy.zeros(len(alltimes))
    for j in range(0,len(alltimes)):
        if(len(alltimes[j]) > 0):
            if(absTime):
                times[j]=abs(alltimes[j][0])
            else:
                times[j]=alltimes[j][0]
        else:
            times[j]=999

    jetTruth=list()
    jetTimes=list()
    for j in range(0,len(trutharrays)):
        if((len(trutharrays[j]) > 0) and (times[j] < 100)):
            frac=trutharrays[j][0]
            eta=abs(alleta[j][0])
            pt=allpt[j][0]
            if((eta >= minEta) and (eta <= maxEta) and (pt >= minPt) and (pt <= maxPt)):
                jetTruth.append(htype(frac))
                jetTimes.append(times[j])
            elif(keepAll):
                jetTruth.append(-1)
                jetTimes.append(times[j])
        elif (keepAll):
            jetTruth.append(-1)
            jetTimes.append(times[j])

    return numpy.array(jetTruth),numpy.array(jetTimes)

def timeFractionTruthBranch(filename,keepAll=False,minEta=0,maxEta=10,minPt=0,maxPt=1000,rcut=0.3,cut=0.075,ptCentering=False):
    leaves=['jtruth','jtime','jeta','jphi','jpt','j0cleta','j0clphi','j0cltime','j0clpixelID','j0clpt']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']
    pixarrays=array['j0clpixelID']
    alltimes=array['jtime']
    alleta=array['jeta']
    allphi=array['jphi']
    allpt=array['jpt']

    etas=array['j0cleta']
    phis=array['j0clphi']
    pts=array['j0clpt']
    timearrays=array['j0cltime']

    jetTruth=list()
    timefrac=list()
    for j in range(0,len(trutharrays)):
        times=abs(timearrays[j])
        if(len(times) > 0):
            frac=trutharrays[j][0]
            eta=abs(alleta[j][0])
            pt=allpt[j][0]
            if((eta >= minEta) and (eta <= maxEta) and (pt >= minPt) and (pt <= maxPt)):
                jetTruth.append(htype(frac))
                pnums=pixarrays[j]
                if(ptCentering):
                    pti=numpy.argmax(pts[j])
                    deta=etas[j]-etas[j][pti]
                    dphi=phis[j]-phis[j][pti]
                else:
                    deta=etas[j]-alleta[j][0]
                    dphi=phis[j]-allphi[j][0]
                dist=sqrt(pow(deta,2)+pow(dphi,2))
                gpts=numpy.where((times < cut) & (dist < rcut) & (pnums != 0))
                dpts=numpy.where(dist < rcut)
                if(len(dpts[0]) > 0):
                    timefrac.append(1-float(len(gpts[0]))/float(len(dpts[0])))
                else:
                    timefrac.append(-1)
            elif(keepAll):
                jetTruth.append(-1)
                timefrac.append(-1)
        elif (keepAll):
            jetTruth.append(-1)
            timefrac.append(-1)

    return numpy.array(jetTruth),numpy.array(timefrac)

def chargeParticleDistance(filename,minEta=0,maxEta=10,minPt=0,maxPt=1000):
    leaves=['jtruth','jtime','jeta','jphi','jpt','j0cleta','j0clphi','j0clpixelID','j0clpt','j0cltime']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']
    alltimes=array['jtime']
    pixarrays=array['j0clpixelID']
    alleta=array['jeta']
    allphi=array['jphi']
    allpt=array['jpt']

    etas=array['j0cleta']
    phis=array['j0clphi']
    pts=array['j0clpt']

    timearrays=array['j0cltime']

    times=numpy.zeros(len(alltimes))
    for j in range(0,len(alltimes)):
        if(len(alltimes[j]) > 0):
            times[j]=abs(alltimes[j][0])
        else:
            times[j]=999

    jetTruth=list()
    jetTimes=list()
    dists=list()
    for j in range(0,len(trutharrays)):
        if((len(trutharrays[j]) > 0) and (times[j] < 100)):
            frac=trutharrays[j][0]
            eta=abs(alleta[j][0])
            pt=allpt[j][0]
            if((eta >= minEta) and (eta <= maxEta) and (pt >= minPt) and (pt <= maxPt)):
                jetTruth.append(htype(frac))
                pnums=pixarrays[j]
                pti=numpy.argmax(pts[j])
                deta=etas[j]-etas[j][pti]
                dphi=phis[j]-phis[j][pti]
                dist=sqrt(pow(deta,2)+pow(dphi,2))
                gpts=numpy.where(pnums != 0)
                pixdists=dist[gpts]
                pixtimes=abs(timearrays[j][gpts])
                if(len(pixdists) > 0):
                    mini=numpy.argmin(pixdists)
                    dists.append(pixdists[mini])
                    jetTimes.append(pixtimes[mini])
                else:
                    dists.append(-1)
                    jetTimes.append(-1)

    return numpy.array(jetTruth),numpy.array(dists),numpy.array(jetTimes)

def ptFractionTruthBranchFrac(filename,keepAll=False,minEta=0,maxEta=10,minPt=0,maxPt=100000):
    leaves=['jtruth','jtime','jeta','jpt']
    array = root2rec(filename,'tree',leaves)

    trutharrays=array['jtruth']
    alltimes=array['jtime']
    alleta=array['jeta']
    allpt=array['jpt']

    times=numpy.zeros(len(alltimes))
    for j in range(0,len(alltimes)):
        if(len(alltimes[j]) > 0):
            times[j]=abs(alltimes[j][0])
        else:
            times[j]=999

    jetTruth=list()
    jetTimes=list()
    jetFrac=list()
    for j in range(0,len(trutharrays)):
        if((len(trutharrays[j]) > 0) and (times[j] < 100)):
            frac=trutharrays[j][0]
            eta=abs(alleta[j][0])
            pt=allpt[j][0]
            if((eta >= minEta) and (eta <= maxEta) and (pt >= minPt) and (pt <= maxPt)):
                jetFrac.append(frac)
                jetTruth.append(htype(frac))
                jetTimes.append(times[j])
            elif(keepAll):
                jetFrac.append(frac)
                jetTruth.append(-1)
                jetTimes.append(times[j])
        elif (keepAll):
            jetFrac.append(-1)
            jetTruth.append(-1)
            jetTimes.append(times[j])

    return numpy.array(jetTruth),numpy.array(jetTimes),numpy.array(jetFrac)

def ptFractionTruth(filename,keepAll=False):
    leaves=['j0clpt','j0cltruth','jtime','jpt']
    array = root2rec(filename,'tree',leaves)
    
    trutharrays=array['j0cltruth']
    ptarrays=array['j0clpt']
    alltimes=array['jtime']
    jetpt=array['jpt']

    times=numpy.zeros(len(alltimes))
    for j in range(0,len(alltimes)):
        if(len(alltimes[j]) > 0):
            times[j]=abs(alltimes[j][0])
        else:
            times[j]=999

    jetTruth=list()
    jetTimes=list()
    for j in range(0,len(trutharrays)):
        if((len(jetpt[j]) > 0) and (times[j] < 100)):
            truth=trutharrays[j]
            pt=ptarrays[j]
            jpt=jetpt[j][0]
            frac=numpy.sum(((1.0-truth)*pt)/jpt)
            jetTruth.append(htype(frac))
            jetTimes.append(times[j])
        elif (keepAll):
            jetTruth.append(-1)
            jetTimes.append(times[j])

    return numpy.array(jetTruth),numpy.array(jetTimes)

def jetDistanceTruth(filename,keepAll=False):
    leaves=['jtime','jeta','jphi','truejeta','truejphi']
    array = root2rec(filename,'tree',leaves)

    jetEta=array['jeta']
    jetPhi=array['jphi']
    trueEta=array['truejeta']
    truePhi=array['truejphi']
    alltimes=array['jtime']

    times=numpy.zeros(len(alltimes))
    for j in range(0,len(alltimes)):
        if(len(alltimes[j]) > 0):
            times[j]=abs(alltimes[j][0])
        else:
            times[j]=999

    jetTruth=list()
    jetTimes=list()
    dists=list()
    for j in range(0,len(jetEta)):
        if((len(jetEta[j]) > 0) and (len(trueEta[j]) > 0) and (times[j] < 100)):
            mindist=10
            for i in range(0,len(truePhi[j])):
                dphi=abs(jetPhi[j][0]-truePhi[j][i])
                if(dphi > 3.14159):
                    dphi=dphi-2*(3.14159)
                dist=sqrt(pow(jetEta[j][0]-trueEta[j][i],2)+pow(dphi,2))
                if(dist < mindist):
                    mindist=dist
            dist=mindist
            if (dist < 0.3):
                jetTruth.append(1)
                jetTimes.append(times[j])
            elif (dist > 0.6):
                jetTruth.append(0)
                jetTimes.append(times[j])
            elif (keepAll):
                jetTruth.append(-1)
                jetTimes.append(times[j])
            if (dist > 1.0):
                dist=1.0
            dists.append(dist)
        elif (keepAll):
            jetTruth.append(-1)
            jetTimes.append(times[j])
            dists.append(-1)

    return numpy.array(jetTruth),numpy.array(jetTimes),numpy.array(dists)

def highestPtPixelTimes(filename,PixelSize,z):
    leaves=['j0cltime','j0clpt','j0cleta','j0clphi']
    array = root2rec(filename,'tree',leaves)

    Aeta=array['j0cleta']
    Aphi=array['j0clphi']
    Apt=array['j0clpt']
    Atimes=array['j0cltime']

    jetTimes=list()
    for i in range(0,len(Apt)):
        eta=Aeta[i]
        phi=Aphi[i]
        zs=z*eta/numpy.abs(eta)
        x,y=xy(eta,phi,zs)
        pt=Apt[i]
        times=Atimes[i]
        gpts=numpy.where(times > -100)
        if(len(gpts[0]) > 0):
            x=x[gpts]
            y=y[gpts]
            pt=pt[gpts]
            times=times[gpts]
            if(len(pt) > 0):
                mpt=numpy.argmax(pt)
                minx=floor(x[mpt]/PixelSize)
                maxx=ceil(x[mpt]/PixelSize)
                miny=floor(y[mpt]/PixelSize)
                maxy=ceil(y[mpt]/PixelSize)
                selpts=numpy.where((y > miny) & (y <= maxy) & (x>minx) & (x<=maxx))
                if(len(selpts[0]) > 0):
                    stimes=times[selpts]
                    mstimes=numpy.mean(stimes)
                    tstd=numpy.std(stimes)
                    rpts=numpy.where((stimes < mstimes+tstd*2) & (stimes > mstimes-tstd*2))
                    bestTimes=stimes[rpts]
                    if(len(bestTimes) > 0):
                        jetTimes.append(numpy.mean(bestTimes))
    print(jetTimes)
    return numpy.array(jetTimes)

def plotROC(truth,variable,color='black',linestyle='-',label=''):

    if(len(truth) != len(variable)):
        print "Warning: truth and variable arrays have different lengths"

    vMax=max(variable)
    dv=vMax/1000

    truthvars=variable[numpy.where(truth == 1)]
    fakevars=variable[numpy.where(truth == 0)]
    truthTotal=float(len(truthvars))
    fakeTotal=float(len(fakevars))

    print len(truthvars),len(fakevars),len(variable)

    cuts=numpy.arange(1e-10,vMax,dv)
    efficiency=list()
    fakerate=list()
    efficiency.append(0)
    fakerate.append(0)
    for cut in cuts:
        tnum=float(len(numpy.where(truthvars < cut)[0]))
        fnum=float(len(numpy.where(fakevars < cut)[0]))
        efficiency.append(tnum/truthTotal)
        fakerate.append(fnum/fakeTotal)
    efficiency.append(1)
    fakerate.append(1)

    if(label != ''):
        plot(efficiency,fakerate,linestyle,label=label,color=color)
    else:
        plot(efficiency,fakerate,linestyle,color=color)

def plotROCSmear(truth,variable,std,color='black',linestyle='-',label=''):
    vs=variable+numpy.random.normal(0,std,len(variable))
    plotROC(truth,vs,color=color,linestyle=linestyle,label=label)

def resolution(sigma,layers,efficiency):
    hits=numpy.arange(1,layers+1)
    probability=numpy.power(efficiency,hits)*numpy.power((1-efficiency),layers-hits)
    resolutions=probability*sigma/numpy.sqrt(hits)
    return resolutions.sum()/probability.sum()

def resolutionArray(sigma,layers,efficiencies):
    results=list()
    for i in range(0,len(efficiencies)):
        results.append(resolution(sigma,layers,efficiencies[i]))
    return numpy.array(results)

def bins (dat): return int(ceil(pow(len(dat),0.333)))
def binSize (dat): return 2.4*std(dat)/pow(len(dat),0.333)
def nbins (dat): return int(ceil((max(dat)-min(dat))/binSize(dat)))
def nbins2 (dat): return int(ceil((max(dat)-min(dat))/binSize(dat)/4))
def contour2d (x,y,color='grey',colors='black',showImage=False):
    scatter(x,y,marker=".",s=1,alpha=0.25,color=color)
    H, xedges, yedges = histogram2d(y, x, bins=(nbins2(y), nbins2(x)),normed=True)
    H=H*binSize(x)*binSize(y)*16
    scale=numpy.max(H)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    levels = (0.9*scale,0.75*scale,0.5*scale,0.25*scale,0.1*scale)
    cset = contour(H, levels, origin="lower",linewidths=(2.0,1.75, 1.5, 1.25, 1.0),colors=colors,extent=extent)
    clabel(cset, inline=1, fontsize=10, fmt="%0.3f")
    for c in cset.collections:
        c.set_linestyle("solid")
    if(showImage):
        imshow(H, interpolation='nearest', origin='low',extent=extent,cmap=cm.binary)

colors=['black','blue','red','green','orange','magenta','cyan','brown']

def y(eta):
    return 3.5/numpy.sinh(eta)

def eta_y(y):
    return numpy.arcsinh(3.5/y)

def etaBinsize(sigma,etaMin=2.5,etaMax=4.3):
    xmin=y(etaMax)
    xmax=y(etaMin)
    x=numpy.arange(xmin,xmax-sigma,sigma)
    binsizeArray=eta_y(x)-eta_y(x+sigma)
    etaArray=eta_y(x)
    return etaArray,binsizeArray

def plotEvent(filename,eventnum=0):
    leaves=['j0cleta','j0clphi','j0cltruth','j0cltime','j0clpixelID','jeta','jphi','jtruth','truejeta','truejphi']

    print(filename)
    array = root2rec(filename,'tree',leaves)

    etas=array['j0cleta']
    phis=array['j0clphi']
    truths=array['j0cltruth']
    times=array['j0cltime']
    ids=array['j0clpixelID']
    jetTruth=array['jtruth']

    frac=array['jtruth'][eventnum][0]
    jetEta=array['jeta'][eventnum]
    jetPhi=array['jphi'][eventnum]
    tjetEta=array['truejeta'][eventnum]
    tjetPhi=array['truejphi'][eventnum]
    
    eta=etas[eventnum]
    phi=phis[eventnum]
    truth=truths[eventnum]
    zerotimes=abs(times[eventnum])
    ID=ids[eventnum].astype(float64)
    
    pixpts=numpy.where((zerotimes < 100) & (ID != 0))
    print(len(pixpts[0]))
    etapix=eta[pixpts]
    phipix=phi[pixpts]
    
    inds=numpy.where((zerotimes < 100) & (ID == 0))
    zerotimes=zerotimes[inds]
    truth=truth[inds]
    eta=eta[inds]
    phi=phi[inds]
    ID=ID[inds]
    
    hpts=numpy.where(truth == 0)
    ppts=numpy.where(truth == 1)
    
    for j in range(0,len(jetEta)):
        xy=(jetEta[j],jetPhi[j])
        c=Circle(xy,radius=0.4,facecolor='green',alpha=0.1)
        gca().add_patch(c)
    for j in range(0,len(tjetEta)):
        xy=(tjetEta[j],tjetPhi[j])
        c=Circle(xy,radius=0.4,facecolor='blue',alpha=0.1)
        gca().add_patch(c)
    scatter(jetEta,jetPhi,color='green',label="Reco Jet")
    scatter(tjetEta,tjetPhi,color='blue',label="Truth Jet")
    scatter(eta[ppts],phi[ppts],marker='.',color='black',label='pileup')
    scatter(eta[hpts],phi[hpts],marker='.',color='red',label='hard scatter')    
    if(len(pixpts[0]) > 0):
        scatter(etapix,phipix,marker='o',alpha=0.25,s=30,color='grey',label='rebinned')

    R=0.4
    limit=R+0.3
    if(frac >= 0.5):
        jtype="Hard Scatter"
    else:
        jtype="Pileup"
    text(jetEta[0]-limit+0.1,jetPhi[0]-limit+0.1,jtype,fontsize=16)
    xlabel('$\eta$')
    ylabel('$\phi$')
    xlim(jetEta[0]-limit,jetEta[0]+limit)
    ylim(jetPhi[0]-limit,jetPhi[0]+limit+0.2)
    legend(loc='upper left',numpoints=1,frameon=False,ncol=2)
