//-*-c++-*-

#ifndef  TimingAnalysis_H
#define  TimingAnalysis_H

#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "Configuration.h"
#include "Definitions.h"
#include "TimingInfo.h"
#include "TimingTracker.h"

using namespace std;
using namespace fastjet;

typedef vector<float> timingBranch;

class TimingDistribution{
 private:
  float _bunchsize;
  int _seed;
  mt19937 rng;  // mt19937 is a standard mersenne_twister_engine

  double _phi;
  double _psi;
  double _phi_nums[2];
  double _psi_nums[2];

  double _gauss_norm;
  double _square_norm;

  double probability(double zpos, double time, distribution dtype);
  int randomSeed();

 public:
  TimingDistribution(float bunchsize=0.075, int seed=-1, double phi=0.0, double psi=0.0);
  void phi(double phi);
  void psi(double psi);
  double psi(){return _psi;};
  double phi(){return _phi;};
  pair<double,double> get(distribution dtype=gaussian);
  double uniform(double min=0, double max=1);
};

class TimingAnalysis{
 private:
  int  ftest;
  bool  fDebug;
  string fOutName;
  
  TFile *tF;
  TTree *tT;
  unique_ptr<TimingDistribution> rnd;

  //these need to be ptrs because constructor must be called
  unique_ptr<JetDefinition> jetDef;
  unique_ptr<AreaDefinition> active_area;
  unique_ptr<GridMedianBackgroundEstimator> bge;
  unique_ptr<Selector> select_fwd;
  unique_ptr<TimingTracker> tracker;
  
  float bunchsize;
  float _minEta;
  float _maxEta;
  double _R;
  double _pixelSize;

  distribution _dtype;
  double psi;
  double phi;
  double timefractioncut;
  
  // Tree Vars ---------------------------------------
  int fTEventNumber;
  int fTNPV;
  float fzvtxspread;
  float ftvtxspread;
  
  timingBranch *jpt;
  timingBranch *jphi;
  timingBranch *jeta;
  timingBranch *jtime;
  timingBranch *jabstime;
  timingBranch *jtruth;
  timingBranch *j0clpt;
  timingBranch *j0clphi;
  timingBranch *j0cleta;
  timingBranch *j0cltime;
  timingBranch *j0clabstime;
  timingBranch *j0cltruth;
  timingBranch *j0clpixelID;
  timingBranch *j0clpixelNum;
  
  timingBranch *truejpt;
  timingBranch *truejphi;
  timingBranch *truejeta;
  timingBranch *truejtime;
  timingBranch *truejabstime;

  bool randomZ;
  bool randomT;
  bool smear;
  bool displace;
  bool segmentation;
  timingMode timeMode;

  Pythia8::Pythia *_pythiaHS;
  Pythia8::Pythia *_pythiaPU;
  
  void DeclareBranches();
  void ResetBranches();
  void FillTree(JetVector jets);
  void FillTruthTree(JetVector jets);
  double ComputeTime(PseudoJet jet, double &abstime);
  double TruthFrac(PseudoJet jet);
  double TimeFrac(PseudoJet jet);
  bool Ignore(Pythia8::Particle &p);

  //Jet selection functions
  void selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets);
  
 public:
  TimingAnalysis (Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q);
  ~TimingAnalysis ();
  
  void AnalyzeEvent(int iEvt, int NPV);
  void Initialize(float minEta, float maxEta, distribution dtype=gaussian,int seed=123);
  
  //settings (call before initialization)
  void Debug(int debug){fDebug = debug;}
  void Bunchsize(float bunchsize_){bunchsize=(bunchsize_>0)?bunchsize_:0;}
  void PileupMode(smearMode PU);
  void SignalMode(smearMode HS);
  void Phi(double phi);
  void Psi(double psi);
  void SetOutName(string outname){fOutName = outname;}
  
};

#endif

