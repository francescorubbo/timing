#ifndef  TimingAnalysis_H
#define  TimingAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/PseudoJet.hh"  

#include "TFile.h"
#include "TTree.h"

#include "Pythia8/Pythia.h"

#include "TH2F.h"
#include "TRandom3.h"

using namespace std;
using namespace fastjet;

class TimingAnalysis{
    private:
        int  ftest;
        int  fDebug;
        string fOutName;

        TFile *tF;
        TTree *tT;
	TRandom3* rnd;

	float bunchsize;

        // Tree Vars ---------------------------------------
        int              fTEventNumber;
	int fTNPV;
	float fzvtxspread;

	std::vector<float> *jpt;
	std::vector<float> *jphi;
	std::vector<float> *jeta;
	std::vector<float> *jtime;
	std::vector<float> *j0clpt;
	std::vector<float> *j0clphi;
	std::vector<float> *j0cleta;
	std::vector<float> *j0cltime;

	std::vector<float> *truejpt;
	std::vector<float> *truejphi;
	std::vector<float> *truejeta;
	std::vector<float> *truejtime;

	bool randomZ;
	bool randomT;

    public:
        TimingAnalysis (float bunchsize_=0.075, bool randomZ_=true, bool randomT_=false);
        ~TimingAnalysis ();
        
        void Begin();
        void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8,  Pythia8::Pythia *pythia_MB, int NPV, float minEta);
        void End();
        void DeclareBranches();
        void ResetBranches();
        void Debug(int debug){
            fDebug = debug;
        }
        void SetOutName(string outname){
            fOutName = outname;
        }
       
       	void FillTree(vector<fastjet::PseudoJet> jets);
       	void FillTruthTree(vector<fastjet::PseudoJet> jets);

	std::pair<double,double> GetVtxZandT();
	double GetIPprob(double zpos, double time);
	double ComputeTime(fastjet::PseudoJet jet);

};

#endif

