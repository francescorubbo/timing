import rootpy
from rootpy.plotting.style import get_style, set_style
from root_numpy import root2rec
from pylab import *
import numpy
from numpy import sin,cos,arctan2,sinh,arcsinh,pi

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
    pnumOut=list()
    for i in range(0,len(pnum)):
        for j in range(0,len(pnum[i])):
            if(pnum[i][j] > -1):
                pnumOut.append(pnum[i][j])

    return numpy.array(pnumOut)

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

def ptFractionTruthBranch(filename,keepAll=False,minEta=0,maxEta=10,minPt=0,maxPt=1000):
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

def ptFractionTruthBranchFrac(filename,keepAll=False,minEta=0,maxEta=10,minPt=0,maxPt=1000):
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
