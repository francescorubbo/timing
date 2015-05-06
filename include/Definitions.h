//-*-c++-*-

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define _USE_MATH_DEFINES

#include <utility>
#include <vector>
#include <memory>
#include <random>
#include <set>
#include <math.h>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

enum distribution {gaussian, pseudoRectangular, crabKissingGaussian, crabKissingSquare};
enum timingMode {highestPT,centralParticle,mean};
enum smearMode {Off,Z,T,ZT};

const double LIGHTSPEED = 299792458.;
const double PI = M_PI;

typedef std::pair<double,double> CorrInfo; // (time,eta)
typedef std::vector<fastjet::PseudoJet> JetVector;

#endif
