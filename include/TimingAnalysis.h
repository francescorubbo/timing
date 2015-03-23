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

        // Tree Vars ---------------------------------------
        int              fTEventNumber;
	int fTNPV;
	float fzvtxspread;
	float bunchsize;

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

    public:
        TimingAnalysis ();
        ~TimingAnalysis ();
        
        void Begin();
        void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8,  Pythia8::Pythia *pythia_MB, int NPV, float zspread, float minEta);
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

	double GetIPprob(double zpos, double time);
	double ComputeTime(fastjet::PseudoJet jet);

};

#endif

