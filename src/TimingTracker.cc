#include "TimingTracker.h"

coordinate xy_to_EtaPhi(double x, double y, double radius){
  double phi = atan2(y,x);
  phi = (phi < 0) ? phi+2*PI : phi;
  double eta=asinh(radius*cos(phi)/x);
  eta = (eta < 0) ? -eta : eta;
  return make_pair(eta,phi);
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

  double width =floor(radius/pixelSize);
  double xnum =floor(xMin/pixelSize)+width;
  double ynum =floor(yMin/pixelSize)+width;
  pixelID_forward = xnum+ynum*width+1; //add 1 such that all IDs > 0
  pixelID_backward = -pixelID_forward;
}

void TrackerPixel::detect(fastjet::PseudoJet &p){
  particles.push_back(p);
}

void TrackerPixel::getParticles(JetVector &detParticles){

  static const double pt=1e-10;
  unsigned long fnum=0;
  unsigned long bnum=0;
  unsigned long snum=0;

  for(auto itr = particles.begin(); itr != particles.end(); ++itr){
    fastjet::PseudoJet p;
    double newEta=_eta;
    double pixelID=pixelID_forward;
    if(itr->eta() < 0){
      newEta*=-1;
      pixelID=pixelID_backward;
      snum=bnum;
      bnum++;
    }
    else{
      snum=fnum;
      fnum++;
    }
    p.reset_PtYPhiM(pt, newEta, _phi);
    p.set_user_info(new TimingInfo(itr->user_info<TimingInfo>().pdg_id(),
				   itr->user_info<TimingInfo>().pythia_id(),
				   itr->user_info<TimingInfo>().pv(),
				   itr->user_info<TimingInfo>().pileup(),
				   itr->user_info<TimingInfo>().time(),
				   itr->user_info<TimingInfo>().abstime(),
				   pixelID,
				   snum));
    detParticles.push_back(p);    
  }
}

pixelCoordinate TimingTracker::getPixel(double eta, double phi){
  coordinate xy = EtaPhi_to_xy(eta,phi,_radius);
  return make_pair(static_cast<int>(floor((xy.first)/_pixelSize)),static_cast<int>(floor((xy.second)/_pixelSize)));
}

TimingTracker::TimingTracker(double pixelSize, double radius){
  if(pixelSize < 1e-7){
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

  for(auto pixel = pixels.begin(); pixel != pixels.end(); ++pixel){
    pixel->second->getParticles(detectedParticles);
  }
  
}

void TimingTracker::AddDetectedParticles(JetVector &truthParticles){
  DetectedParticles(truthParticles,truthParticles);
}
