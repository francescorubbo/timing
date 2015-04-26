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

double sgn(double val){
  if(val < 0)
    return -1;
  else if (val > 0)
    return 1;
  else
    return 0;
}

// Constructor 
TimingAnalysis::TimingAnalysis(Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q){

  if((pythiaHS != NULL) and (pythiaPU != NULL)){
    _pythiaHS=pythiaHS;
    _pythiaPU=pythiaPU;
  }
  else{
    cerr << "Invalid Pythia pointer passed to TimingAnalysis" << endl;
    exit(1);
  }

  fDebug=q.fDebug;
  if(fDebug) 
    cout << "TimingAnalysis::TimingAnalysis Start " << endl;

  ftest = 0;
  fOutName = q.outName;
  
  Bunchsize(q.bunchsize);
  PileupMode(q.PUmode);
  SignalMode(q.HSmode);
  Psi(q.psi);
  Phi(q.phi);
  
  if(fDebug) 
    cout << "TimingAnalysis::TimingAnalysis End " << endl;
  
  //suppress fastjet banner
  fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);
}

void TimingAnalysis::PileupMode(smearMode PU){
  switch(PU){
  case Off:
    randomT=false;
    randomZ=false;
    break;
  case Z:
    randomT=false;
    randomZ=true;
    break;
  case T:
    randomT=true;
    randomZ=false;
    break;
  case ZT:
    randomT=true;
    randomZ=false;
  }
}

void TimingAnalysis::SignalMode(smearMode HS){
  switch(HS){
  case Off:
    smear=false;
    displace=false;
    break;
  case Z:
    smear=false;
    displace=true;
    break;
  case T:
    smear=true;
    displace=false;
    break;
  case ZT:
    smear=true;
    displace=false;
  }
}

void TimingAnalysis::Phi(double phi){
  this->phi= phi;
}

void TimingAnalysis::Psi(double psi){
  this->psi= psi;
}

// Destructor 
TimingAnalysis::~TimingAnalysis(){
  if(tT != NULL){
    tT->Write();
    tF->Close();
    delete tF;
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
   rnd.reset(new TimingDistribution(bunchsize,seed,phi,psi));
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
void TimingAnalysis::AnalyzeEvent(int ievt, int NPV, float minEta, float maxEta){

  if(fDebug) 
    cout << "TimingAnalysis::AnalyzeEvent Begin " << endl;
  
  // -------------------------
  if (!_pythiaHS->next()) return;
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
  
  double zhs=0.0;
  double ths=0.0;
  if(displace) //randomly distribute "ideally measured" hard-scatter vertex
    zhs=fzvtxspread;
  if(smear)
    ths=ftvtxspread;

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
    for (int i = 0; i < _pythiaPU->event.size(); ++i) {

      if(Ignore(_pythiaPU->event[i]))
	continue;
      
      //Instantiate new pseudojet
      PseudoJet p(_pythiaPU->event[i].px(), 
		  _pythiaPU->event[i].py(), 
		  _pythiaPU->event[i].pz(),
		  _pythiaPU->event[i].e() ); 
      
      //extract event information
      double eta = p.rapidity();

      //extract event information
      double eta = p.rapidity();
      double sinheta = sinh(eta);
      double cosheta = cosh(eta);
      
      //calculate eta from displacement (minEta pos)
      static const double radius = 1.2; // barrel radius=1.2 meter
      double zbase = radius*sinh(minEta)*sgn(eta); //displace due to new location of Hard-Scatter Vertex
      double corrEta = asinh(zbase*sinheta/(zbase-zvtx));
      if(fabs(corrEta)>maxEta) continue;
      
      //calculate time measured relative to if event was at 0
      double dist = (zbase-zvtx)*cosheta/sinheta;
      double time = fabs(dist)/LIGHTSPEED + tvtx; //plus random time
      
      double refdist = sqrt(pow(dist,2)+pow((zbase-zhs),2)-pow((zbase-zvtx),2))
      double reftime = fabs(refdist)/LIGHTSPEED;
      double corrtime = (time-reftime)*1e9;
      if(fabs(corrEta)<minEta) 
	corrtime = -999.;
      
      p.reset_PtYPhiM(p.pt(), corrEta, p.phi(), 0.);
      
      p.set_user_info(new TimingInfo(_pythiaPU->event[i].id(),i,iPU,true,corrtime)); 
      particlesForJets.push_back(p); 
    }
    if (!_pythiaPU->next()) continue;
  }
  
  // Particle loop -----------------------------------------------------------
  for (int ip=0; ip<_pythiaHS->event.size(); ++ip){
    
    if(Ignore(_pythiaHS->event[ip]))
      continue;
    
    fastjet::PseudoJet p(_pythiaHS->event[ip].px(), 
			 _pythiaHS->event[ip].py(), 
			 _pythiaHS->event[ip].pz(),
			 _pythiaHS->event[ip].e() ); 
    
    double eta = p.rapidity();

    if (fabs(eta)>maxEta) continue;
    double corrtime = ths*1e9;
    if (fabs(eta)<minEta) corrtime = -999.;
    p.reset_PtYPhiM(p.pt(), eta, p.phi(), 0.);
    p.set_user_info(new TimingInfo(_pythiaHS->event[ip].id(),ip,0, false,corrtime)); //0 for the primary vertex. 
    
    p.reset_PtYPhiM(p.pt(), corrEta, p.phi(), 0.);
    //0 for the primary vertex.
    p.set_user_info(new TimingInfo(_pythiaHS->event[ip].id(),ip,0, false,corrtime));  
    
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

bool TimingAnalysis::Ignore(Pythia8::Particle &p){
  if (!p.isFinal() )      
    return true;
  switch(abs(p.id())){
  case 12:
  case 13:
  case 14:
  case 16:
    return true;
  default:
    return false;
  }
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
  double ampl=0;
  
  switch(dtype){
  case gaussian:
    return _gauss_norm*exp(-( pow(zpos,2) + pow(LIGHTSPEED*time,2) ) / (pow(_bunchsize,2)));
  case crabKissingGaussian:
    ampl=_gauss_norm*_phi_nums[1]*_psi_nums[1];
    return ampl*exp(-((pow(zpos,2)*_phi_nums[0]) + (pow(LIGHTSPEED*time,2)*_psi_nums[0])) / (pow(_bunchsize,2)));
  case pseudoRectangular:
    return _square_norm*exp(-(4*PI*PI/pow((tgamma(0.25)*_bunchsize),4))*(pow(zpos,4)+pow(LIGHTSPEED*time,4)+6*pow(LIGHTSPEED*time*zpos,2)));
  case crabKissingSquare:
    ampl=_square_norm*_phi_nums[1]*_psi_nums[1];
    return _square_norm*exp(-(4*PI*PI/pow((tgamma(0.25)*_bunchsize),4))*(pow(zpos,4)+pow(LIGHTSPEED*time,4)+6*pow(LIGHTSPEED*time*zpos,2)))*exp(-(pow(zpos*_phi,2) + pow(LIGHTSPEED*time*_psi,2)) / (pow(_bunchsize,2)));
  default:
    cerr << "Invalid RNG Distribution" << endl;
    exit(10);
  }
}

int TimingDistribution::randomSeed(){
  int timeSeed = time(NULL);
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729); 
}

TimingDistribution::TimingDistribution(float bunchsize, int seed, double phi, double psi) : _bunchsize(bunchsize){
  if(seed == -1){
    cout << "Timing Distribution Generating Random Seed" << endl;
    _seed=randomSeed();
  }
  else
    _seed=seed;

  _gauss_norm=LIGHTSPEED/(PI*_bunchsize*_bunchsize);
  _square_norm=pow(2,3.5)*LIGHTSPEED*PI/(pow(tgamma(0.25),4)*pow(_bunchsize,2));
  
  rng.seed(_seed);  
  this->phi(phi);
  this->psi(psi);
}

void TimingDistribution::phi(double phi){
  _phi=phi;
  _phi_nums[0]=1+pow(phi,2);
  _phi_nums[1]=sqrt(_phi_nums[0]);
}

void TimingDistribution::psi(double psi){
  _psi=psi;
  _psi_nums[0]=1+pow(psi,2);
  _psi_nums[1]=sqrt(_psi_nums[0]);
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
