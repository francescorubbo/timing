#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "TString.h"
#include "TSystem.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "TimingAnalysis.h"

#include "boost/program_options.hpp"

using namespace std;
namespace po = boost::program_options;

void printBanner();
void printOptions(po::variables_map vm);

int getSeed(int seed){                                                                                                                                                                     
  if (seed > -1) return seed;                                                                                                                                                           
  int timeSeed = time(NULL);                                                                 
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}

int main(int argc, char* argv[]){
    // arguments 
    string outName   = "Timing.root";
    int    pileup    = 0;
    float  bunchsize = 0.075;
    float  minEta    = 2.5;
    int    nEvents   = 1;
    int    fDebug    = 1;
    float  pThatmin  =100;
    float  pThatmax  =500;
    float  boson_mass=1500;
    int    proc      = 4;
    int    seed      =-1;

    bool   randZ     =true;
    bool   randT     =false;
    bool   smearHS   =false;
    bool   useCK     =false;
    float  phi       =0;
    float  psi       =0;

    po::options_description gen_desc("Allowed options");
    gen_desc.add_options()
      ("help", "produce help message")
      ("Debug",     po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile",   po::value<string>(&outName)->default_value("Timing.root"), "output file name")
      ("Seed",      po::value<int>(&seed)->default_value(-1), "Seed. -1 means random seed");

    po::options_description sim_flag("Simulation Flags");
    sim_flag.add_options()
      ("VaryZ",     "Vary only Z Vertex of Pileup")
      ("VaryT",     "Vary only Vertex Time of Pileup")
      ("VaryZT",    "Vary both Z and Time of Vertex of Pileup")
      ("SmearHS",   "Smear Hard Scatter Vertex in Time")
      ("ForceCK",   "Force Crab-Kissing PDF even if Phi=Psi=0");

    po::options_description sim_desc("Simulation Settings");
    sim_desc.add_options()
      ("NEvents",   po::value<int>(&nEvents)->default_value(1) ,    "Number of Events ")
      ("Pileup",    po::value<int>(&pileup)->default_value(0), "Number of Additional Interactions")
      ("BunchSize", po::value<float>(&bunchsize)->default_value(0.075), "Size of Proton Bunches")
      ("Phi",     po::value<float>(&phi)->default_value(0), "Phi Parameter, Crab-Kissing PDF")
      ("Psi",     po::value<float>(&psi)->default_value(0), "Psi Parameter, Crab-Kissing PDF")
      ("MinEta",    po::value<float>(&minEta)->default_value(2.5), "Minimum Pseudorapidity for Particles")
      ("Proc",      po::value<int>(&proc)->default_value(4), "Process:\n - 1: Z'T->ttbar\n - 2: W'->WZ+lept\n - 3: W'->WZ+had\n - 4: QCD")
      ("pThatMin",  po::value<float>(&pThatmin)->default_value(100), "pThatMin for QCD")
      ("pThatMax",  po::value<float>(&pThatmax)->default_value(500), "pThatMax for QCD")
      ("BosonMass", po::value<float>(&boson_mass)->default_value(1500), "Z' or W' mass in GeV")
      ;

    po::options_description desc;
    desc.add(gen_desc).add(sim_flag).add(sim_desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }
    else{
      printBanner();
      printOptions(vm);
    }

    cout << "\t";
    if (vm.count("SmearHS")>0){
      cout <<"Smearing Hard-Scatter Timing" << endl;
      smearHS=true;
    }
    else
      cout <<"No Hard-Scatter Smearing" << endl;
    
    cout << "\t";
    if (vm.count("VaryZ")>0){
      cout <<"Varying Z of Pileup vertex" << endl;
      randZ=true;
      randT=false;
    }
    else if (vm.count("VaryT")>0){
      cout <<"Varying T of Pileup vertex" <<endl;
      randZ=false;
      randT=true;
    }
    else if ((vm.count("VaryZT")>0) or ((vm.count("VaryZ")>0) and vm.count("VaryT"))){
      cout <<"Varying Z and T of Pileup vertex" <<endl;
      randZ=true;
      randT=true;
    }
    else{
      cout <<"No Pileup Smearing" << endl;
      randZ=false;
      randT=false;
    }

    if ((vm.count("ForceCK")>0) or (phi != 0) or (psi != 0)){
      cout << "\tUsing Crab-Kissing PDF" << endl;
      useCK=true;
    }
    else
      cout << "\tUsing Gaussian PDF" << endl;

    cout << endl;
    
    //seed 
    seed = getSeed(seed);
    
    // Configure and initialize pythia
    std::unique_ptr<Pythia8::Pythia> pythia8(new Pythia8::Pythia("../xmldoc",false));
    pythia8->readString("Print:quiet=on");
    pythia8->readString("Random:setSeed = on"); 
    std::stringstream ss; 
    ss << "Random:seed = " << seed;
    if(fDebug)
      cout << ss.str() << endl; 
    pythia8->readString(ss.str());

   if(proc ==1){
     std::stringstream bosonmass_str; bosonmass_str<< "32:m0=" << boson_mass ;
     pythia8->readString(bosonmass_str.str());
     pythia8->readString("NewGaugeBoson:ffbar2gmZZprime= on");
     pythia8->readString("Zprime:gmZmode=3");
     pythia8->readString("32:onMode = off");
     pythia8->readString("32:onIfAny = 6");
     pythia8->readString("24:onMode = off");
     pythia8->readString("24:onIfAny = 1 2 3 4");
     pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line! 
   }else if(proc ==2){
      std::stringstream bosonmass_str; bosonmass_str<< "34:m0=" << boson_mass ;
      pythia8->readString(bosonmass_str.str());
      pythia8->readString("NewGaugeBoson:ffbar2Wprime = on");
      pythia8->readString("Wprime:coup2WZ=1");
      pythia8->readString("34:onMode = off");
      pythia8->readString("34:onIfAny = 23 24");
      pythia8->readString("24:onMode = off");
      pythia8->readString("24:onIfAny = 1 2 3 4");
      pythia8->readString("23:onMode = off");
      pythia8->readString("23:onIfAny = 12");
      pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 3){
      std::stringstream bosonmass_str; bosonmass_str<< "34:m0=" << boson_mass ;
      pythia8->readString(bosonmass_str.str());
      pythia8->readString("NewGaugeBoson:ffbar2Wprime = on");
      pythia8->readString("Wprime:coup2WZ=1");
      pythia8->readString("34:onMode = off");
      pythia8->readString("34:onIfAny = 23 24");
      pythia8->readString("24:onMode = off");
      pythia8->readString("24:onIfAny = 11 12");
      pythia8->readString("23:onMode = off");
      pythia8->readString("23:onIfAny = 1 2 3 4 5");
      pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 4){ 
      pythia8->readString("HardQCD:all = on");
      std::stringstream ptHatMin;
      std::stringstream ptHatMax;
      ptHatMin << "PhaseSpace:pTHatMin  =" << pThatmin;
      ptHatMax << "PhaseSpace:pTHatMax  =" << pThatmax;
      pythia8->readString(ptHatMin.str());
      pythia8->readString(ptHatMax.str());
      pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else{ throw std::invalid_argument("received invalid 'process'");}

   //Setup the pileup
   std::unique_ptr<Pythia8::Pythia> pythia_MB(new Pythia8::Pythia("../xmldoc",false));
   pythia_MB->readString("Random:setSeed = on");   
   ss.clear(); 
   ss.str(""); 
   ss << "Random:seed = " << seed+1; 
   if(fDebug)
     cout << ss.str() << endl; 
   pythia_MB->readString(ss.str());
   pythia_MB->readString("Print:quiet=on");
   pythia_MB->readString("SoftQCD:nonDiffractive = on");
   pythia_MB->readString("HardQCD:all = off");
   pythia_MB->readString("PhaseSpace:pTHatMin  = .1");
   pythia_MB->readString("PhaseSpace:pTHatMax  = 20000");
   pythia_MB->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */);

   // TimingAnalysis
   TimingAnalysis analysis(bunchsize,randZ,randT,smearHS);
   analysis.SetOutName(outName);

   //determine which PDF should be used
   if(useCK){
     analysis.Initialize(distribution::crabKissing,2*seed,phi,psi);
   }
   else
     analysis.Initialize(distribution::gaussian,2*seed); //seeds vertex generator, different seed than pythia
   
   analysis.Debug(fDebug);

   std::cout << "Number of Pileup Events: " << pileup << std::endl;

   // Event loop
   cout << "Progress:" << endl;
   for (Int_t iev = 0; iev < nEvents; iev++) {
     if (iev%20==0)
       cout << "\tCurrent: " << iev << endl;
     analysis.AnalyzeEvent(iev, pythia8.get(), pythia_MB.get(), pileup, minEta);
   }

   cout << "Timing Analysis Complete!" << endl;

   return 0;
}

void printBanner(){
  cout << endl << "=================================================================" << endl;
  cout << "=                        Timing Analysis                        =" << endl;
  cout << "=================================================================" << endl << endl;
}

void printOptions(po::variables_map vm){
  cout << "Settings:" << endl;
  for (po::variables_map::const_iterator itr=vm.begin();itr != vm.end();++itr){
    printf("%15s\t",itr->first.c_str());
    if((itr->first == "VaryT") or (itr->first == "VaryZ") or (itr->first == "VaryZT") or (itr->first == "ForceCK") or (itr->first == "SmearHS")){
      cout << endl;
      continue;
    }

    try { 
      cout << "= " << itr->second.as<double>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<float>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<int>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<std::string>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
  }
}
