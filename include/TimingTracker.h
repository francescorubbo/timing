//-*-c++-*-

#ifndef TIMINGTRACKER_H
#define TIMINGTRACKER_H

#include "Definitions.h"
#include "TimingInfo.h"

using namespace std;

typedef pair<int,int> pixelCoordinate;
typedef pair<double,double> coordinate;

coordinate xy_to_EtaPhi(double x, double y, double radius);
coordinate EtaPhi_to_xy(double eta, double phi, double radius);

class TrackerPixel{
 private:
  JetVector particles;
  JetVector detParticles;
  double _phi;
  double _eta;
  unsigned long pixelID;
 public:
  TrackerPixel(double xMin, double yMin, double radius, double pixelSize);
  void detect(fastjet::PseudoJet &p);
  JetVector& getParticles();
};

class TimingTracker{
 private:
  double _pixelSize;
  double _radius;
  map<pixelCoordinate,shared_ptr<TrackerPixel> > pixels;
  pixelCoordinate getPixel(double eta, double phi);
 public:
  TimingTracker(double pixelSize, double radius);
  void DetectedParticles(JetVector &truthParticles, JetVector &detectedParticles);
};

#endif
