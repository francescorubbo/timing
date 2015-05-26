#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "Definitions.h"
#include "TimingAnalysis.h"
#include "Configuration.h"

using namespace std;

int main(int argc, char* argv[]){

  //read options, configure Pythia
  Configuration settings(argc, argv);
  std::shared_ptr<Pythia8::Pythia> pythiaHS(new Pythia8::Pythia("../xmldoc",false));
  std::shared_ptr<Pythia8::Pythia> pythiaPU(new Pythia8::Pythia("../xmldoc",false));
  settings.ConfigurePythiaSignal(pythiaHS.get());
  settings.ConfigurePythiaPileup(pythiaPU.get());
    
  // TimingAnalysis
  TimingAnalysis analysis(pythiaHS.get(),pythiaPU.get(), settings);
  analysis.Initialize(settings.minEta,settings.maxEta,settings.dtype,2*settings.seed);
  
  // Event loop
  cout << "Progress:" << endl;
  for (Int_t iev = 0; iev < settings.nEvents; iev++) {
    if (iev%10==0)
      cout << "\tCurrent: " << iev << endl;
    analysis.AnalyzeEvent(iev, settings.pileup);
  }  
  cout << "Timing Analysis Complete!" << endl;
  
  return 0;
}

