//-*-c++-*-

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"
#include "Pythia8/Pythia.h"

using namespace std;
namespace po = boost::program_options;

void printBanner();
void printOptions(po::variables_map vm);

struct PythiaSettings{
  PythiaSettings(float m,float pmin,float pmax,float s){
    bosonMass=m;
    pthatmin=pmin;
    pthatmax=pmax;
    seed=s;
  }
  float bosonMass;
  float pthatmin;
  flaot pthatmax;
  float seed;
};

void configurePythia(Pythia8::Pythia* pu, Pythia8::Pythia* hs, int proc, PythiaSettings settings);

struct Configuration{
  Configuration(int argc, char* argv[]);
  
  string outName;
  int    pileup;
  float  bunchsize;
  float  minEta;
  int    nEvents;
  int    fDebug;
  float  pThatmin;
  float  pThatmax;
  float  boson_mass;
  int    proc;
  int    seed;
  smearMode HSmode;
  smearMode PUmode;
  bool   useCK;
  float  phi;
  float  psi;
  int    profile;
};

#endif
