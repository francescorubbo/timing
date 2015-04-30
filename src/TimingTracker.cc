#include "TimingTracker.h"

coordinate xy_to_EtaPhi(double x, double y, double radius){
  double phi = atan(y/x);
  return make_pair(asinh(radius*cos(phi)/x),phi);
}

coordinate EtaPhi_to_xy(double eta, double phi, double radius){
  double height=radius/sinh(eta);
  return make_pair(height*cos(phi),height*sin(phi));
}

TrackerPixel::TrackerPixel(double xMin, double yMin, double radius, double pixelSize){
  coordinate ll = xy_to_EtaPhi(xMin,yMin,radius);
  coordinate rr = xy_to_EtaPhi(xMin+pixelSize,yMin+pixelSize,radius);

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

  detParticles.clear();

  double pt,time;
  int id;
  bool pileup;

  int snum=0;
  for(auto itr = particles.begin(); itr != particles.end(); ++itr){
    
    if(snum == 0)
      pt=itr->pt();
    else
      pt=1e-100;

    time=itr->user_info<TimingInfo>().time();
    pileup=itr->user_info<TimingInfo>().pileup();
    id=itr->user_info<TimingInfo>().pdg_id();

    fastjet::PseudoJet p;
    p.reset_PtYPhiM(pt, _eta, _phi);
    p.set_user_info(new TimingInfo(id,snum,pixelID,pileup,time));
    detParticles.push_back(p);
    
    snum++;
  }

  return detParticles;
}

pixelCoordinate TimingTracker::getPixel(double eta, double phi){
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

  pixels.clear();

  //fill tracker
  pixelCoordinate pi;
  for(auto particle = truthParticles.begin(); particle != truthParticles.end(); ++particle){
    pi=getPixel(particle->eta(),particle->phi());
    if(pixels.count(pi) == 0){
      double xMin = static_cast<double>(pi.first)*_pixelSize;
      double yMin = static_cast<double>(pi.second)*_pixelSize;
      pixels[pi].reset(new TrackerPixel(xMin,yMin,_radius,_pixelSize));
    }
    pixels[pi]->detect(*particle);
  }

  detectedParticles.clear();

  for(auto pixel = pixels.begin(); pixel != pixels.end(); ++pixel){
    JetVector particles=pixel->second->getParticles();
    for(auto particle = particles.begin(); particle != particles.end(); ++particle)
      detectedParticles.push_back(*particle);
  }
  
}
