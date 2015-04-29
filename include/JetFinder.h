//-*-c++-*-

#ifndef JETFINDER_H
#define JETFINDER_H

#include <map>
#include <vector>
#include <math.h>
#include <memory>

#include "Definitions.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

using namespace fastjet;

typedef vector<fastjet::PseudoJet> JetVector;

class JetFinder{
private:
  double _pixelSize;
  double _R;
  double _etaMin;
  double _etaMax;

  //these need to be ptrs because constructor must be called
  unique_ptr<JetDefinition> jetDef;
  unique_ptr<AreaDefinition> active_area;
  unique_ptr<GridMedianBackgroundEstimator> bge;
  unique_ptr<Selector> select_fwd;

 public:
  JetFinder(double etaMin, double etaMax, double pixelSize=1e-4, double R=0.4);
  void PixelSize(double pixelSize);
  void selectJets(JetVector particlesForJets, JetVector selectedJets);
  void selectSegmentedJets(JetVector particlesForJets, JetVector selectedJets);
};

#endif
