//-*-c++-*-

#ifndef  TimingAnalysis_H
#define  TimingAnalysis_H

#include <vector>
#include <math.h>
#include <string>
#include <random>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "fastjet/PseudoJet.hh"  
#include "TFile.h"
#include "TTree.h"
#include "Pythia8/Pythia.h"
#include "TH2F.h"

#include "Configuration.h"
#include "Definitions.h"

using namespace std;
using namespace fastjet;

typedef vector<float> timingBranch;
typedef vector<fastjet::PseudoJet> JetVector;

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
  
  float bunchsize;
  distribution _dtype;
  double psi;
  double phi;
  
  // Tree Vars ---------------------------------------
  int fTEventNumber;
  int fTNPV;
  float fzvtxspread;
  float ftvtxspread;
  
  timingBranch *jpt;
  timingBranch *jphi;
  timingBranch *jeta;
  timingBranch *jtime;
  timingBranch *j0clpt;
  timingBranch *j0clphi;
  timingBranch *j0cleta;
  timingBranch *j0cltime;
  timingBranch *j0cltruth;
  
  timingBranch *truejpt;
  timingBranch *truejphi;
  timingBranch *truejeta;
  timingBranch *truejtime;
  
  bool randomZ;
  bool randomT;
  bool smear;
  bool displace;

  Pythia8::Pythia *_pythiaHS;
  Pythia8::Pythia *_pythiaPU;
  
  void DeclareBranches();
  void ResetBranches();
  void FillTree(JetVector jets);
  void FillTruthTree(JetVector jets);
  CorrInfo ComputeTime(float z, float eta, float mineta);
  double ComputeTime(PseudoJet jet);
  bool Ignore(Particle &p);

  //void BinParticles(float binScale, JetVector input, JetVector output);
  
 public:
  TimingAnalysis (Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q);
  ~TimingAnalysis ();
  
  void AnalyzeEvent(int iEvt, int NPV, float minEta, float maxEta);
  void Initialize(distribution dtype=gaussian,int seed=123);
  
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

