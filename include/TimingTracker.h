//-*-c++-*-

#ifndef TIMINGTRACKER_H
#define TIMINGTRACKER_H

#include "Definitions.h"
#include "TimingInfo.h"

using namespace std;

typedef pair<int,int> pixelCoordinate;
typedef pair<double,double> coordinate;

coordinate xy_to_EtaPhi(double x, double y, double z);
coordinate EtaPhi_to_xy(double eta, double phi, double z);

class TrackerPixel{
 private:
  JetVector particles;
  JetVector detParticles;
  double _phi;
  double _eta;
  double pixelID_forward;
  double pixelID_backward;
 public:
  TrackerPixel(double xMin, double yMin, double radius, double zbase, double pixelSize);
  void detect(fastjet::PseudoJet &p);
  void getParticles(JetVector &detParticles);
};

class TimingTracker{
 private:
  bool _filterCharge;
  bool _filterP;
  double _minP;
  double _pixelSize;
  double _radius;
  double _zbase;
  map<pixelCoordinate,shared_ptr<TrackerPixel> > pixels;
  pixelCoordinate getPixel(double eta, double phi);
 public:
  TimingTracker(double pixelSize, double radius, double zbase, bool filterCharge=true);
  void SetPThreshold(double minP);
  void DetectedParticles(JetVector &truthParticles, JetVector &detectedParticles);
  void AddDetectedParticles(JetVector &truthParticles);
};

#endif
