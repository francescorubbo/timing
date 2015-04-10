#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TF2.h"

#include "TimingAnalysis.h"
#include "TimingInfo.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include "Pythia8/Pythia.h"
#include "TROOT.h"

using namespace std;

const double LIGHTSPEED = 299792458.;
const double PI  =3.141592653589793238463;

double sgn(double val){
  if(val < 0)
    return -1;
  else if (val > 0)
    return 1;
  else
    return 0;
}

// Constructor 
TimingAnalysis::TimingAnalysis(float bunchsize_, bool randomZ_, bool randomT_, bool smear_,bool Debug){

    fDebug=Debug;

    if(fDebug) 
      cout << "TimingAnalysis::TimingAnalysis Start " << endl;
    ftest = 0;
    fOutName = "test.root";

    bunchsize = bunchsize_;
    randomZ=randomZ_;
    randomT=randomT_;
    smear=smear_;

    if(fDebug) 
      cout << "TimingAnalysis::TimingAnalysis End " << endl;

    //suppress fastjet banner
    fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);
}

// Destructor 
TimingAnalysis::~TimingAnalysis(){
  if(tT != NULL){
    tT->Write();
    tF->Close();
    delete tF;
    delete tT;
  }

  if(jpt != NULL){
    delete jpt;
    delete jphi;
    delete jeta;
    delete jtime;
    delete j0clpt;
    delete j0clphi;
    delete j0cleta;
    delete j0cltime;
    delete j0cltruth;
    delete truejpt;
    delete truejphi;
    delete truejeta;
    delete truejtime;
  }
}

// Begin method
void TimingAnalysis::Initialize(distribution dtype, int seed){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("tree", "Event Tree for Timing");
   rnd.reset(new TimingDistribution(bunchsize,seed));
   _dtype=dtype;

   // for shit you want to do by hand
   DeclareBranches();

   jpt = new timingBranch();  
   jphi = new timingBranch();  
   jeta = new timingBranch();  
   jtime = new timingBranch();

   j0clpt = new timingBranch();  
   j0clphi = new timingBranch();  
   j0cleta = new timingBranch();  
   j0cltime = new timingBranch();  
   j0cltruth = new timingBranch();
   
   truejpt = new timingBranch();  
   truejphi = new timingBranch();  
   truejeta = new timingBranch();  
   truejtime = new timingBranch();  

   ResetBranches();
   
   return;
}

// Analyze
void TimingAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV,
				  float minEta){
    if(fDebug) 
      cout << "TimingAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) 
      cout << "TimingAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    JetVector particlesForJets;
    JetVector particlesForJets_np;

    //Pileup Loop

    fTNPV = NPV;
    std::pair<double,double> randomVariates;  
    randomVariates=rnd->get(_dtype);
    fzvtxspread = randomVariates.first;
    ftvtxspread = randomVariates.second;
    
    //Loop over Pileup Events
    for (int iPU = 0; iPU <= NPV; ++iPU) {
      
      //determine random vertex position in z-t space                                                                                                                                                     
      randomVariates=rnd->get(_dtype);
      double zvtx = 0;
      double tvtx = 0;
      if(randomZ)
	zvtx = randomVariates.first;
      if(randomT)
	tvtx = randomVariates.second;
      
      //Loop over pileup particles
      for (int i = 0; i < pythia_MB->event.size(); ++i) {

	//skipping Leptons
        if (!pythia_MB->event[i].isFinal()    ) continue;
        if (fabs(pythia_MB->event[i].id())==12) continue;
        if (fabs(pythia_MB->event[i].id())==14) continue;
        if (fabs(pythia_MB->event[i].id())==13) continue;
        if (fabs(pythia_MB->event[i].id())==16) continue;
	
	//Instantiate new pseudojet
	PseudoJet p(pythia_MB->event[i].px(), pythia_MB->event[i].py(), pythia_MB->event[i].pz(),pythia_MB->event[i].e() ); 
	
	//extract event information
	double eta = p.rapidity();
	double sinheta = sinh(eta);
	double cosheta = cosh(eta);

	//calculate eta from displacement (minEta pos)
	double radius = 1.2; // barrel radius=1.2 meter
	double zbase = radius*sinh(minEta)*sgn(eta); //should eta be minEta?
	double corrEta = asinh(zbase*sinheta/(zbase-zvtx));
	if(fabs(corrEta)>5.0) continue;

	//calculate time measured relative to if event was at 0
	double dist = (zbase-zvtx)*cosheta/sinheta;
	double time = fabs(dist)/LIGHTSPEED + tvtx; //plus random time
	double refdist = zbase*cosh(corrEta)/sinh(corrEta);
	double reftime = fabs(refdist)/LIGHTSPEED;
	double corrtime = (time-reftime)*1e9;
	if(fabs(corrEta)<minEta) 
	  corrtime = -999.;

	p.reset_PtYPhiM(p.pt(), corrEta, p.phi(), 0.);
	
	p.set_user_info(new TimingInfo(pythia_MB->event[i].id(),i,iPU,true,corrtime)); 
	particlesForJets.push_back(p); 
      }
      if (!pythia_MB->next()) continue;
    }

    //determine random vertex position in z-t space                                                                                                                             

    double tvtx=0.0;
    if(smear)
      tvtx=rnd->get(_dtype).second;
    
    // Particle loop -----------------------------------------------------------
    for (int ip=0; ip<pythia8->event.size(); ++ip){
      
      if (!pythia8->event[ip].isFinal() )      continue;
      //if (fabs(pythia8->event[ip].id())  ==11) continue;
      if (fabs(pythia8->event[ip].id())  ==12) continue;
      if (fabs(pythia8->event[ip].id())  ==13) continue;
      if (fabs(pythia8->event[ip].id())  ==14) continue;
      if (fabs(pythia8->event[ip].id())  ==16) continue;

      fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() ); 
      double eta = p.rapidity();
      if (fabs(eta)>5.0) continue;
      double corrtime = tvtx*1e9;
      if (fabs(eta)<minEta) corrtime = -999.;
      p.reset_PtYPhiM(p.pt(), eta, p.phi(), 0.);
      p.set_user_info(new TimingInfo(pythia8->event[ip].id(),ip,0, false,corrtime)); //0 for the primary vertex. 
      
      particlesForJets.push_back(p);
      particlesForJets_np.push_back(p);

    } // end particle loop -----------------------------------------------

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
    fastjet::AreaDefinition active_area(fastjet::active_area);
    fastjet::ClusterSequenceArea clustSeq(particlesForJets, jetDef, active_area);

    fastjet::GridMedianBackgroundEstimator bge(4.5,0.6);
    bge.set_particles(particlesForJets);
    fastjet::Subtractor subtractor(&bge);

    JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
    JetVector subtractedJets = subtractor(inclusiveJets);
    Selector select_fwd = SelectorAbsRapRange(minEta,4.5);
    JetVector selectedJets = select_fwd(subtractedJets);

    FillTree(selectedJets);

    fastjet::ClusterSequenceArea clustSeqTruth(particlesForJets_np, jetDef, active_area);
    JetVector truthJets = sorted_by_pt(clustSeqTruth.inclusive_jets(10.));
    JetVector selectedTruthJets = select_fwd(truthJets);
    FillTruthTree(selectedTruthJets);

    tT->Fill();

    if(fDebug) 
      cout << "TimingAnalysis::AnalyzeEvent End " << endl;

    return;
}

// worker function to actually perform an analysis
void TimingAnalysis::FillTree(JetVector jets){  

  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    jpt->push_back(jets[ijet].pt());
    jeta->push_back(jets[ijet].eta());
    jphi->push_back(jets[ijet].phi());
    jtime->push_back(ComputeTime(jets[ijet]));
  }

  if(jets.size()>0)
    for (unsigned int icl=0; icl<jets[0].constituents().size(); icl++){    
      j0clpt->push_back(jets[0].constituents()[icl].pt());
      j0clphi->push_back(jets[0].constituents()[icl].phi());
      j0cleta->push_back(jets[0].constituents()[icl].eta());
      j0cltime->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().time());
      j0cltruth->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().pileup() ? 1.0 : 0.0);
    }  
}

void TimingAnalysis::FillTruthTree(JetVector jets){  
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    truejpt->push_back(jets[ijet].pt());
    truejeta->push_back(jets[ijet].eta());
    truejphi->push_back(jets[ijet].phi());
    truejtime->push_back(ComputeTime(jets[ijet]));
  }
}

double TimingAnalysis::ComputeTime(fastjet::PseudoJet jet){
  double time = 0.;
  float maxpt = 0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    // time+=jet.constituents()[i].user_info<TimingInfo>().time();
    float pt = jet.constituents()[i].pt();
    if(pt>maxpt){
      maxpt = pt;
      time = jet.constituents()[i].user_info<TimingInfo>().time();
    }//endif
  }//endloop
  return time;
}

// declare branches
void TimingAnalysis::DeclareBranches(){
   // Event Properties 
  gROOT->ProcessLine("#include <vector>");
  
  tT->Branch("EventNumber",&fTEventNumber,"EventNumber/I");
  tT->Branch("NPV",&fTNPV,"NPV/I");
  tT->Branch("zvtxspread",&fzvtxspread,"zvtxspread/F");
  tT->Branch("tvtxspread",&ftvtxspread,"tvtxspread/F"); 
  tT->Branch("jpt","std::vector<float>",&jpt);
  tT->Branch("jphi","std::vector<float>",&jphi);
  tT->Branch("jeta","std::vector<float>",&jeta);
  tT->Branch("jtime","std::vector<float>",&jtime);
  tT->Branch("j0clpt","std::vector<float>",&j0clpt);
  tT->Branch("j0clphi","std::vector<float>",&j0clphi);
  tT->Branch("j0cleta","std::vector<float>",&j0cleta);
  tT->Branch("j0cltime","std::vector<float>",&j0cltime);
  tT->Branch("j0cltruth","std::vector<float>",&j0cltruth);

  tT->Branch("truejpt", "std::vector<float>",&truejpt);
  tT->Branch("truejphi","std::vector<float>",&truejphi);
  tT->Branch("truejeta","std::vector<float>",&truejeta);
  tT->Branch("truejtime","std::vector<float>",&truejtime);
  
  return;
}

// resets vars
void TimingAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV = -1;
      fzvtxspread = -1;
      ftvtxspread = -1;

      jpt->clear();
      jphi->clear();
      jeta->clear();
      jtime->clear();
      j0clpt->clear();
      j0clphi->clear();
      j0cleta->clear();
      j0cltime->clear();
      j0cltruth->clear();
      truejpt->clear();
      truejphi->clear();
      truejeta->clear();
      truejtime->clear();

}

double TimingDistribution::probability(double zpos, double time, distribution dtype){
  switch(dtype){
  case gaussian:
    static double ampl = LIGHTSPEED/(PI*_bunchsize*_bunchsize);
    return ampl*exp(-( pow(zpos,2) + pow(LIGHTSPEED*time,2) ) / (pow(_bunchsize,2)));
  case crabKissing:
    return 1; //no implemented yet
  default:
    cerr << "Invalid RNG Distribution" << endl;
    exit(10);
  }
}

int TimingDistribution::randomSeed(){
  int timeSeed = time(NULL);                                                                                                                                                                          \
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729); 
}

TimingDistribution::TimingDistribution(float bunchsize, int seed){
  if(seed == -1){
    cout << "Timing Distribution Generating Random Seed" << endl;
    _seed=randomSeed();
  }
  else
    _seed=seed;
  
  rng.seed(_seed);  
  _bunchsize=bunchsize;
}

pair<double,double> TimingDistribution::get(distribution dtype){
  
  double maxprob = probability(0,0,dtype);
  double zpos;
  double time;
  
  while(true){
    zpos = uniform(-3.0*_bunchsize,3.0*_bunchsize);
    time = uniform(-3*_bunchsize/LIGHTSPEED,3*_bunchsize/LIGHTSPEED);
    if(probability(zpos,time,dtype) > maxprob*uniform())
      break;
  }
  cout << "Out: " << zpos << '\t' << time << endl;
  return std::make_pair(zpos,time);
}

double TimingDistribution::uniform(double min, double max){
  static unsigned int range_min=rng.min();
  static unsigned int range_max=rng.max();

  unsigned int rn=rng();
  double drn=static_cast<double>(rn-range_min)/static_cast<double>(range_max-range_min);
  if((max != 1) and (min != 0))
    return (max-min)*drn+min;
  else
    return drn;
}
