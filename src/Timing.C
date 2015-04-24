#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "Pythia8/Pythia.h"
#include "TimingAnalysis.h"
#include "Configuration.h"

using namespace std;

int getSeed(int seed){                                                      
  if (seed > -1) return seed;
  int timeSeed = time(NULL);                                                                 
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}

int main(int argc, char* argv[]){

  Configuration q(argc, argv);

  //obtain random seed or preset seed
  q.seed = getSeed(q.seed);
  std::shared_ptr<Pythia8::Pythia> pythiaHS(new Pythia8::Pythia("../xmldoc",false));
  std::shared_ptr<Pythia8::Pythia> pythiaPU(new Pythia8::Pythia("../xmldoc",false));
  PythiaSettings settings(q.boson_mass,q.pThatmin,q.pThatmax,q.seed);
  ConfigurePythia(pythiaHS.get(),pythiaPU.get(),q.proc,settings);
    
  // TimingAnalysis
  TimingAnalysis analysis(pythiaHS.get(),pythiaPU.get(),q.bunchsize,q.PUmode,q.HSmode,q.fDebug);
  analysis.SetOutName(q.outName);
  //seeds vertex generator, different seed than pythia
  analysis.Initialize(q.dtype,2*q.seed,q.phi,q.psi);
  
  std::cout << "Number of Pileup Events: " << q.pileup << std::endl;
  
  // Event loop
  cout << "Progress:" << endl;
  for (Int_t iev = 0; iev < q.nEvents; iev++) {
    if (iev%20==0)
      cout << "\tCurrent: " << iev << endl;
    analysis.AnalyzeEvent(iev, q.pileup, q.minEta);
  }
  
  cout << "Timing Analysis Complete!" << endl;
  
  return 0;
}

