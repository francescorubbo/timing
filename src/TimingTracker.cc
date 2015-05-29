#include "TimingTracker.h"

coordinate xy_to_EtaPhi(double x, double y, double z){
  double phi = atan2(y,x);
  phi = (phi < 0) ? phi+2*PI : phi;
  double eta=asinh(z*cos(phi)/x);
  eta = (eta < 0) ? -eta : eta;
  return make_pair(eta,phi);
}

coordinate EtaPhi_to_xy(double eta, double phi, double z){
  double height=z/sinh(eta);
  return make_pair(height*cos(phi),height*sin(phi));
}

TrackerPixel::TrackerPixel(double xMin, double yMin, double radius, double zbase, double pixelSize){
  coordinate ll = xy_to_EtaPhi(xMin,yMin,zbase);
  coordinate rr = xy_to_EtaPhi(xMin+pixelSize,yMin+pixelSize,zbase);

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

  //ghost particle pt low enough not to affect Clustering, high enough 
  //to not cause error
  static const double pt=1e-10;
  unsigned long fnum=0;
  unsigned long bnum=0;
  unsigned long snum=0;

  double minTimes[] = {999,999};
  int iminTimes[] = {-1,-1};
  //loop through all detected particles
  for(int i=0;i<static_cast<int>(particles.size());i++){
    double time=particles[i].user_info<TimingInfo>().time();
    //if backwards
    if(particles[i].eta() < 0){
      bnum++; //count particle
      //decide if it was detected first
      if(time < minTimes[0]){
	minTimes[0]=time;
	iminTimes[0]=i;
      }
    }
    //if forwards
    else{
      fnum++; //count particle
      //decide if it was detected first
      if(time < minTimes[1]){
	minTimes[1]=time;
	iminTimes[1]=i;
      }
    }
  }

  //loop over first backwards and forwards particles
  for(int i=0;i<2;i++){
    //if we found a minimum time particle
    if(iminTimes[i] > -1){
      fastjet::PseudoJet p;
      double newEta=_eta;
      double pixelID=pixelID_forward;
      //if backwards
      if(particles[iminTimes[i]].eta() < 0){
	newEta*=-1;
	pixelID=pixelID_backward;
	snum=bnum;
      }
      //if forwards
      else{
	snum=fnum;
      }
      //create particle
      p.reset_PtYPhiM(pt, newEta, _phi);
      p.set_user_info(new TimingInfo(particles[iminTimes[i]].user_info<TimingInfo>().pdg_id(),
				     particles[iminTimes[i]].user_info<TimingInfo>().charge(),
				     particles[iminTimes[i]].user_info<TimingInfo>().pythia_id(),
				     particles[iminTimes[i]].user_info<TimingInfo>().pv(),
				     particles[iminTimes[i]].user_info<TimingInfo>().pileup(),
				     particles[iminTimes[i]].user_info<TimingInfo>().pt(),
				     particles[iminTimes[i]].user_info<TimingInfo>().time(),
				     particles[iminTimes[i]].user_info<TimingInfo>().abstime(),
				     pixelID,
				     snum));
      //return particle
      detParticles.push_back(p);    
    }
  }
}

pixelCoordinate TimingTracker::getPixel(double eta, double phi){
  coordinate xy = EtaPhi_to_xy(abs(eta),phi,_zbase); //pretend all on forward tracker to save memory
  return make_pair(static_cast<int>(floor((xy.first)/_pixelSize)),static_cast<int>(floor((xy.second)/_pixelSize)));
}

TimingTracker::TimingTracker(double pixelSize, double radius, double zbase, bool filterCharge){
  if(pixelSize < 1e-7){
    cerr << "Invalid pixel size " << pixelSize << endl;
    exit(30);
  }
  _pixelSize=pixelSize;
  _radius=radius;
  _zbase=zbase;
  _filterCharge=filterCharge;

  _filterP=false;
  _minP=0;
}

void TimingTracker::SetPThreshold(double minP){
  if(minP > 0){
    _filterP=true;
    _minP=minP;
  }
  else
    cout << "TimingTracker::SetPThreshold - Warning: Invalid Momentum" << endl;
}

void TimingTracker::DetectedParticles(JetVector &truthParticles, JetVector &detectedParticles){

  pixels.clear();

  //fill tracker
  pixelCoordinate pi;
  for(auto particle = truthParticles.begin(); particle != truthParticles.end(); ++particle){
    
    if(_filterCharge and (particle->user_info<TimingInfo>().charge()==0)) 
      continue; //use only charged particles
    if(_filterP and (sqrt(particle->modp2()) < _minP))
      continue; //use only higher |p| particles

    pi=getPixel(particle->eta(),particle->phi());
    if(pixels.count(pi) == 0){
      double xMin = static_cast<double>(pi.first)*_pixelSize;
      double yMin = static_cast<double>(pi.second)*_pixelSize;
      pixels[pi].reset(new TrackerPixel(xMin,yMin,_radius,_zbase,_pixelSize));
    }
    pixels[pi]->detect(*particle);
  }

  for(auto pixel = pixels.begin(); pixel!= pixels.end(); ++pixel){
    pixel->second->getParticles(detectedParticles);
  }
  
}

void TimingTracker::AddDetectedParticles(JetVector &truthParticles){
  DetectedParticles(truthParticles,truthParticles);
}
