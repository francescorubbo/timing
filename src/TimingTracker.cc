#include "TimingTracker.h"

coordinate xy_to_EtaPhi(double x, double y, double radius){
  double phi = atan(y/x);
  return make_pair(asinh(radius*cos(phi)/x),phi);
}

coordinate EtaPhi_to_xy(double eta, double phi, double radius){
  double height=radius/sinh(eta);
  return make_pair(height*cos(phi),height*sin(phi))
}

TrackerPixel::TrackerPixel(double xMin, double yMin, double radius, double pixelSize){
  coordinate ll = xy_to_EtaPhi(xMin,yMin);
  coordinate rr = xy_to_EtaPhi(xMin+pixelSize,yMin+pixelSize);

  //eta and phi are average of top right and lower left corner values
  _eta=(ll.first+rr.first)/2;
  _phi=(ll.second+ll.second)/2;

  int width = static_cast<int>(radius/pixelSize);
  int xnum=static_cast<int>(xMin/pixelSize)+width;
  int ynum=static_cast<int>(yMin/pixelSize)+width;
  pixelID = xnum+ynum*2*width;
}

void TrackerPixel::detect(fastjet::PseudoJet &p){
  particles.push_back(p);
}

JetVector& TrackerPixel::getParticles(){

  double pt=0;
  double time=0;
  double frac=0;
  double num = static_cast<double>(particles.size());
  for(auto itr = particles.begin(); itr != particles.end(); ++itr){
    pt+=itr->pt();
    TimingInfo info=itr->user_info();
    time+=info.time();
    frac+=info.pileup() ? 1.0 : 0.0;
  }

  if(num > 1){
    time /= num;
    frac /= num;
  }

  //need to change type of TimingInfo to accomodate fraction
  bool tparticle=false;
  if(frac > 0.5)
    tparticle=true;
  
  fastjet::PseudoJet p();
  p.reset_PtYPhiM(pt, _eta, _phi);
  p.set_user_info(new TimingInfo(pixelID,num,0,tparticle,time)); 

  detParticles.push_back(p);

  return detParticles;
}

pixelCoordinate getPixel(double eta, double phi){
  coordinate xy = EtaPhi_to_xy(eta,phi,_radius);
  return make_pair(static_cast<int>(floor((xy.first)/_pixelSize)),static_cast<int>(floor((xy.second)/_pixelSize)));
}

TimingTracker::TimingTracker(double pixelSize, double radius){
  if(pixelSize < 1e-6){
    cerr << "Invalid pixel size " << pixelSize << endl;
    exit(30);
  }
  _pixelSize=pixelSize;
  _radius=radius;
}

void TimingTracker::DetectedParticles(JetVector &truthParticles, JetVector &detectedParticles){

  pixels.reset();

  //fill tracker
  pixelCoordinate pi;
  for(auto particle = truthParticles.begin(); particle != truthParticles.end(); ++particle){
    pi=getPixel(particle->eta(),particle->phi());
    if(pixels.count(pi) == 0){
      double xMin = static_cast<double>(pi.first)*_pixelSize;
      double yMix = static_cast<double>(pi.second)*_pixelSize;
      pixels[pi].reset(new TrackerPixel(xMin,yMin,_radius,_pixelSize));
    }
    pixels[pi]->detect(*particle);
  }

  detectedParticles.reset();

  for(auto pixel = pixels.begin(); pixel != pixels.end(); ++pixels){
    JetVector particles=pixel->getParticles();
    for(auto particle = particles.begin(); particle != particles.end(); ++particle)
      detectedParticles.push_back(*particle);
  }
  
}
