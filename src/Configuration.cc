#include "Configuration.h"

Configuration::Configuration(int argc, char* argv[]){

    HSmode =smearMode::ZT;
    PUmode =smearMode::ZT;
    useCK     =false;
    filterCharge=true;
    int profile;
    int timing;

    po::options_description gen_desc("Allowed options");
    gen_desc.add_options()
      ("help", "produce help message")
      ("Debug",     po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile",   po::value<string>(&outName)->default_value("Timing.root"), "output file name")
      ("Seed",      po::value<int>(&seed)->default_value(-1), "Seed. -1 means random seed");

    po::options_description sim_flag("Simulation Flags");
    sim_flag.add_options()
      ("FixZT",     "Fix both Z and Time of Vertex of Pileup")
      ("VaryZ",     "Vary only Z Vertex of Pileup")
      ("VaryT",     "Vary only Vertex Time of Pileup")
      ("VaryZT",    "Vary both Z and Time of Vertex of Pileup")
      ("noSmearHS", "Smear Hard Scatter Vertex in Time")
      ("SmearHST",  "Smear Hard Scatter Vertex in Time")
      ("SmearHSZ",  "Smear Hard Scatter Vertex in Z (correcting time)")
      ("SmearHSZT", "Smear Hard Scatter Vertex in Time and Z (correcting time)")
      ("ForceCK",   "Force Crab-Kissing PDF even if Phi=Psi=0")
      ("KeepNeutral","Keep neutrals in timing tracker");

    po::options_description sim_desc("Simulation Settings");
    sim_desc.add_options()
      ("NEvents",   po::value<int>(&nEvents)->default_value(1) ,    "Number of Events ")
      ("Pileup",    po::value<int>(&pileup)->default_value(80), "Number of Additional Interactions")
      ("BunchSize", po::value<float>(&bunchsize)->default_value(0.075), "Size of Proton Bunches")
      ("Profile",   po::value<int>(&profile)->default_value(0), "Bunch Profile Type:\n - 0: Gaussian\n - 1: PseudoRectangular")
      ("JetTiming", po::value<int>(&timing)->default_value(0), "Bunch Profile Type:\n - 0: Highest-pT\n - 1: Central Particle\n - 2: Mean Particle Time")
      ("Phi",       po::value<float>(&phi)->default_value(0), "Phi Parameter, Crab-Kissing PDF")
      ("Psi",       po::value<float>(&psi)->default_value(0), "Psi Parameter, Crab-Kissing PDF")
      ("MinEta",    po::value<float>(&minEta)->default_value(2.5), "Minimum Pseudorapidity for Particles")
      ("MaxEta",    po::value<float>(&maxEta)->default_value(4.3), "Minimum Pseudorapidity for Particles")
      ("PixelSize", po::value<float>(&pixelSize)->default_value(0), "Pixel Size for Segmentation (in microns)")
      ("Proc",      po::value<int>(&proc)->default_value(4), "Process:\n - 1: Z'T->ttbar\n - 2: W'->WZ+lept\n - 3: W'->WZ+had\n - 4: QCD")
      ("pThatMin",  po::value<float>(&pThatmin)->default_value(100), "pThatMin for QCD")
      ("pThatMax",  po::value<float>(&pThatmax)->default_value(500), "pThatMax for QCD")
      ("BosonMass", po::value<float>(&boson_mass)->default_value(1500), "Z' or W' mass in GeV");
    
    po::options_description desc;
    desc.add(gen_desc).add(sim_flag).add(sim_desc);
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
	exit(0);
    }

    print();

    cout << "\t";
    if (vm.count("SmearHST")>0){
      cout <<"Smearing Hard-Scatter Timing" << endl;
      HSmode=smearMode::T;
    }
    else if (vm.count("SmearHSZ")>0){
      cout <<"Smearing Hard-Scatter Z (corrected time)" << endl;
      HSmode=smearMode::Z;
    }
    else if(vm.count("SmearHSZT")>0){
      cout <<"Smearing Hard-Scatter Timing and Z (corrected time)" << endl;
      HSmode=smearMode::ZT;
    }
    else if(vm.count("noSmearHS")>0){
      cout <<"No Hard-Scatter Smearing" << endl;
      HSmode=smearMode::Off;
    }
    else{
      cout <<"Smearing Hard-Scatter Timing and Z (corrected time)" << endl;
      HSmode=smearMode::ZT;
    }

    cout << "\t";
    if (vm.count("VaryZ")>0){
      cout <<"Varying Z of Pileup vertex" << endl;
      PUmode=smearMode::Z;
    }
    else if (vm.count("VaryT")>0){
      cout <<"Varying T of Pileup vertex" <<endl;
      PUmode=smearMode::T;
    }
    else if (vm.count("VaryZT")>0){
      cout <<"Varying Z and T of Pileup vertex" <<endl;
      PUmode=smearMode::ZT;
    }
    else if (vm.count("FixZT")>0){
      cout <<"No Pileup Smearing" << endl;
      PUmode=smearMode::Off;
    }
    else{
      cout <<"Varying Z and T of Pileup vertex" <<endl;
      PUmode=smearMode::ZT;
    }

    cout << "\t";
    if(vm.count("KeepNeutral")){
      cout << "Keeping Neutral Particles in Tracker" << endl;
      filterCharge=false;
    }
    else{
      cout << "Discarding Neutral Particles in Tracker" << endl;
      filterCharge=true;
    }

    if((profile < 0) or (profile > 1)){
      cerr << "ERROR: Invalid profile \"" << profile << "\", choices are Gaussian (0) and PseudoRectangular (1)" << endl;
      exit(2);
    }
    
    dtype=distribution::gaussian;
    cout << "\t";
    if ((vm.count("ForceCK")>0) or (phi != 0) or (psi != 0)){
      useCK=true;
      if(profile == 0){
	dtype=distribution::crabKissingGaussian;
	cout << "Using Gaussian Crab-Kissing Bunch Profile" << endl;
      }
      else{
	dtype=distribution::crabKissingSquare;
	cout << "Using PseudoRectangular Crab-Kissing Bunch Profile" << endl;
      }
    }
    else if (profile == 1){
      dtype=distribution::pseudoRectangular;
      cout << "Using PseudoRectangular Bunch Profile" << endl;
    }
    else
      cout << "Using Gaussian Bunch Profile" << endl;

    seed=getSeed(seed);
    
    if(pixelSize > 1e-1){
      pixelSize*=1e-6;
      segmentation=true;
      cout << "\tUsing Finite Segmentation" << endl;
    }
    else{
      segmentation=false;
      cout << "\tUsing Infinite Segmentation" << endl;
    }

    switch(timing){
    case timingMode::highestPT:
      cout << "\tUsing time of highest pT particle" << endl;
      timemode=highestPT;
      break;
    case timingMode::centralParticle:
      cout << "\tUsing time of central particle" << endl;
      timemode=centralParticle;
      break;
    case timingMode::mean:
      cout << "\tUsing mean of particle times" << endl;
      timemode=mean;
      break;
    default:
      cerr << "Invalid timing mode " << timing << endl;
      exit(1);
    }
    
    cout << endl;
}

void Configuration::print(){

  cout << endl << "=================================================================" << endl;
  cout << "=                        Timing Analysis                        =" << endl;
  cout << "=================================================================" << endl << endl;

  cout << "Settings:" << endl;
  for (po::variables_map::const_iterator itr=vm.begin();itr != vm.end();++itr){
    printf("%15s\t",itr->first.c_str());
    
    if(
       (itr->first == "VaryT") 
       or (itr->first == "VaryZ") 
       or (itr->first == "VaryZT") 
       or (itr->first == "ForceCK") 
       or (itr->first == "SmearHSZ")
       or (itr->first == "SmearHST")
       or (itr->first == "SmearHSZT")
       ){
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

void Configuration::ConfigurePythiaSignal(Pythia8::Pythia* hs){

    hs->readString("Print:quiet=on");
    hs->readString("Random:setSeed = on"); 
    std::stringstream ss; 
    ss << "Random:seed = " << seed;
    hs->readString(ss.str());

   if(proc ==1){
     std::stringstream bosonmass_str; 
     bosonmass_str<< "32:m0=" << boson_mass;
     hs->readString(bosonmass_str.str());
     hs->readString("NewGaugeBoson:ffbar2gmZZprime= on");
     hs->readString("Zprime:gmZmode=3");
     hs->readString("32:onMode = off");
     hs->readString("32:onIfAny = 6");
     hs->readString("24:onMode = off");
     hs->readString("24:onIfAny = 1 2 3 4");
     hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line! 
   }else if(proc ==2){
      std::stringstream bosonmass_str; 
      bosonmass_str<< "34:m0=" << boson_mass;
      
      hs->readString(bosonmass_str.str());
      hs->readString("NewGaugeBoson:ffbar2Wprime = on");
      hs->readString("Wprime:coup2WZ=1");
      hs->readString("34:onMode = off");
      hs->readString("34:onIfAny = 23 24");
      hs->readString("24:onMode = off");
      hs->readString("24:onIfAny = 1 2 3 4");
      hs->readString("23:onMode = off");
      hs->readString("23:onIfAny = 12");
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 3){
      std::stringstream bosonmass_str; 
      bosonmass_str<< "34:m0=" << boson_mass;
      hs->readString(bosonmass_str.str());
      hs->readString("NewGaugeBoson:ffbar2Wprime = on");
      hs->readString("Wprime:coup2WZ=1");
      hs->readString("34:onMode = off");
      hs->readString("34:onIfAny = 23 24");
      hs->readString("24:onMode = off");
      hs->readString("24:onIfAny = 11 12");
      hs->readString("23:onMode = off");
      hs->readString("23:onIfAny = 1 2 3 4 5");
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }else if(proc == 4){ 
      hs->readString("HardQCD:all = on");
      std::stringstream ptHatMin;
      std::stringstream ptHatMax;
      ptHatMin << "PhaseSpace:pTHatMin  =" << pThatmin;
      ptHatMax << "PhaseSpace:pTHatMax  =" << pThatmax;
      hs->readString(ptHatMin.str());
      hs->readString(ptHatMax.str());
      hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
   }
   else{ 
     throw std::invalid_argument("received invalid 'process'");
   }
   return;
}

void Configuration::ConfigurePythiaPileup(Pythia8::Pythia* pu){
   //Setup the pileup
   pu->readString("Random:setSeed = on");   
   std::stringstream ss;
   ss.clear(); 
   ss.str(""); 
   ss << "Random:seed = " << seed+1; 
   pu->readString(ss.str());
   pu->readString("Print:quiet=on");
   pu->readString("SoftQCD:nonDiffractive = on");
   pu->readString("HardQCD:all = off");
   pu->readString("PhaseSpace:pTHatMin  = .1");
   pu->readString("PhaseSpace:pTHatMax  = 20000");
   pu->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */);
   return;
}

int Configuration::getSeed(int seed){                                                      
  if (seed > -1) return seed;
  int timeSeed = time(NULL);                                                                 
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}
