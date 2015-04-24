//-*-c++-*-

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"
#include "Pythia8/Pythia.h"

#include "Definitions.h"

using namespace std;
namespace po = boost::program_options;

class Configuration{
private:
  po::variables_map vm;
  int getSeed(int seed);
public:
  Configuration(int argc, char* argv[]);
  void ConfigurePythia(Pythia8::Pythia* pu, Pythia8::Pythia* hs);
  void print();
  
  string outName;

  int    pileup;
  int    nEvents;
  int    fDebug;
  int    proc;
  int    seed;

  float  bunchsize;
  float  minEta;
  float  maxEta;  
  float  pThatmin;
  float  pThatmax;
  float  boson_mass;

  float  phi;
  float  psi;

  bool   useCK;
  smearMode HSmode;
  smearMode PUmode;
  distribution dtype;  
};

#endif
