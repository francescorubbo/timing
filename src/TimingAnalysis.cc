#include "TimingAnalysis.h"

using namespace std;

double sgn(double val){
  if(val < 0)
    return -1;
  else if (val > 0)
    return 1;
  else
    return 0;
}

double distance(double eta, double phi,double truthEta, double truthPhi){
  double de=eta-truthEta;
  double dp=phi-truthPhi;
  return sqrt(pow(de,2)+pow(dp,2));
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

  if(q.segmentation){
    segmentation=true;
    _pixelSize=q.pixelSize;
    filterCharge=q.filterCharge;
    _minP=q.minP;
    filterP=q.filterP;
  }
  else
    segmentation=false;

  timeMode=q.timemode;
  simMagneticField=q.magfield;

  if(fDebug) 
    cout << "TimingAnalysis::TimingAnalysis End " << endl;
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
    randomZ=true;
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
    displace=true;
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
    delete jabstime;
    delete jtruth;
    delete j0clpt;
    delete j0clphi;
    delete j0cleta;
    delete j0cltime;
    delete j0clabstime;
    delete j0cltruth;
    delete j0clpixelID;
    delete j0clpixelNum;
    delete j0clpdgid;
    delete j1clpt;
    delete j1clphi;
    delete j1cleta;
    delete j1clabstime;
    delete j1cltruth;
    delete truejpt;
    delete truejphi;
    delete truejeta;
    delete truejtime;
    delete truejabstime;
  }
}

// Begin method
void TimingAnalysis::Initialize(float minEta, float maxEta, distribution dtype, int seed){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("tree", "Event Tree for Timing");
   rnd.reset(new TimingDistribution(bunchsize,seed,phi,psi));
   _dtype=dtype;

   _minEta= (minEta > 0) ? minEta : 0;
   if(maxEta <= minEta){
     cerr << "Invalid Eta Limits " << minEta << " -> " << maxEta << "Passed to TimingAnalysis::Initialize" << endl;
     exit(20);
   }
   _maxEta= maxEta;

   const double zbase=3.5;
   const double radius=zbase/sinh(_minEta);
   if(segmentation){
     tracker.reset(new TimingTracker(_pixelSize,radius,zbase,filterCharge));
     if(filterP)
       tracker->SetPThreshold(_minP);
   }
   
   const double R=0.4;
   const double grid_spacing(0.6);
   
   //suppress fastjet banner
   fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);
   
   jetDef.reset(new JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best));
   active_area.reset(new AreaDefinition(fastjet::active_area));
   bge.reset(new GridMedianBackgroundEstimator(_maxEta, grid_spacing));
   rescaling.reset(new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05)); //function for rapidity rescaling of rho
   bge->set_rescaling_class(rescaling.get());
   select_fwd.reset(new Selector(SelectorAbsRapRange(_minEta,_maxEta)));
 
   DeclareBranches();
   
   jpt = new timingBranch();  
   jphi = new timingBranch();  
   jeta = new timingBranch();  
   jtime = new timingBranch();
   jabstime = new timingBranch();
   jtruth = new timingBranch();

   j0clpt = new timingBranch();  
   j0clphi = new timingBranch();  
   j0cleta = new timingBranch();  
   j0cltime = new timingBranch();  
   j0clabstime = new timingBranch();
   j0cltruth = new timingBranch();

   j0clpixelID = new timingBranch();
   j0clpixelNum = new timingBranch();
   j0clpdgid = new timingBranch();

   j1clpt = new timingBranch();  
   j1clphi = new timingBranch();  
   j1cleta = new timingBranch();  
   j1clabstime = new timingBranch();
   j1cltruth = new timingBranch();
   
   truejpt = new timingBranch();  
   truejphi = new timingBranch();  
   truejeta = new timingBranch();  
   truejtime = new timingBranch();  
   truejabstime = new timingBranch();

   ResetBranches();
   
   return;
}

// Analyze
void TimingAnalysis::AnalyzeEvent(int ievt, int NPV){

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

  static const double z0 = 3.5; // FIX THIS! DEFINE Z0 ONLY IN ONE PLACE! 

  //Loop over Pileup Events
  for (int iPU = 0; iPU < NPV; ++iPU) {
    
    cout << "filling w/ pileup" << endl;

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

      //apply magnetic field
      
      //extract event information
      double eta = p.rapidity();
      double sinheta = sinh(eta);
      
      //calculate eta from displacement (minEta pos)
      double zbase = z0*sgn(eta); //displace due to new location of Hard-Scatter Vertex
      double dz = zbase-zvtx;
      double corrEta = asinh(zbase*sinheta/dz);
      double hsEta = asinh((zbase-zhs)*sinheta/dz);
      
      //calculate time measured relative to if event was at 0
      double betaz=_pythiaPU->event[i].pz()/_pythiaPU->event[i].e();
      if(betaz > 1.001)
	cout << "Error: Invalid Beta value!!!" << endl;
      double time = fabs(dz)/(betaz*LIGHTSPEED); //plus random time
      time+=tvtx;
      
      double reftime = fabs((zbase-zhs)/(LIGHTSPEED*sinh(hsEta)/cosh(hsEta)));
      double corrtime = (time-reftime)*1e9;
      double corrphi = p.phi();
      if((fabs(corrEta)< _minEta) or (fabs(corrEta)>_maxEta))
	time = corrtime = -999.;
      else if(simMagneticField){
	double charge = _pythiaPU->event[i].charge()*1.609e-19;
	double omega = charge*2./(_pythiaPU->event[i].e()*1.782661845e-27);
	corrphi = p.phi() + omega*time;
	if(abs(corrphi)>2*PI) corrphi = fmod(corrphi,2*PI);
      }
      p.reset_PtYPhiM(p.pt(), corrEta, corrphi);
      
      p.set_user_info(new TimingInfo(_pythiaPU->event[i].id(),_pythiaPU->event[i].charge(),
				     i,iPU,true,_pythiaPU->event[i].pT(),corrtime,time*1e9)); 
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
    double sinheta = sinh(eta);

    //calculate eta from displacement (minEta pos)
    double zbase = z0*sgn(eta); //displace due to new location of Hard-Scatter Vertex    
    double dz = zbase-zhs;
    double betaz=_pythiaHS->event[ip].pz()/_pythiaHS->event[ip].e();
    if(betaz > 1.001)
      cout << "Error: Invalid Beta value!!!" << endl;
    double time = fabs(dz)/(betaz*LIGHTSPEED) + ths; //plus random time  

    double corrEta = asinh(zbase*sinheta/dz);
    double reftime = fabs((dz)/(LIGHTSPEED*sinheta/cosh(eta)));
    double corrtime = (time-reftime)*1e9;

    double corrphi = p.phi();
    if ((fabs(corrEta)<_minEta) or (fabs(corrEta)>_maxEta))
      time = corrtime = -999.;
    else if(simMagneticField){
      double charge = _pythiaHS->event[ip].charge()*1.609e-19;
      double omega = charge*2./(_pythiaHS->event[ip].e()*1.782661845e-27);
      corrphi = p.phi() + omega*time;
      if(abs(corrphi)>2*PI) corrphi = fmod(corrphi,2*PI);
    }
    p.reset_PtYPhiM(p.pt(), corrEta, corrphi);
    //0 for the primary vertex.
    p.set_user_info(new TimingInfo(_pythiaHS->event[ip].id(),_pythiaHS->event[ip].charge(),
				   ip,0, false,_pythiaHS->event[ip].pT(),corrtime,time*1e9));  
    
    particlesForJets.push_back(p);
    particlesForJets_np.push_back(p);
    
  } // end particle loop -----------------------------------------------

  if(segmentation){
    tracker->AddDetectedParticles(particlesForJets);
  }

  JetVector selectedJets,selectedTruthJets;
  fastjet::ClusterSequenceArea clustSeq(particlesForJets, *jetDef, *active_area);
  selectJets(particlesForJets,clustSeq,selectedJets);

  fastjet::ClusterSequenceArea clustSeqTruth(particlesForJets_np, *jetDef, *active_area);
  selectedTruthJets = sorted_by_pt(clustSeqTruth.inclusive_jets(10.));
  
  FillTree(selectedJets,selectedTruthJets);
  FillTruthTree(selectedTruthJets);
  
  tT->Fill();
  
  if(fDebug) 
    cout << "TimingAnalysis::AnalyzeEvent End " << endl;
  
  return;
}

void TimingAnalysis::selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets){
  try{
    bge->set_particles(particlesForJets);

    fastjet::Subtractor subtractor(bge.get());    
    JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
    JetVector subtractedJets = subtractor(inclusiveJets);

    JetVector allSelectedJets;
    allSelectedJets.clear();
    allSelectedJets = (*select_fwd)(subtractedJets);
    
    //select jets with pt > 20
    selectedJets.clear();
    for( auto ijet = allSelectedJets.begin(); ijet != allSelectedJets.end(); ++ijet){
      if(ijet->pt() >= 20)
	selectedJets.push_back(*ijet);
    }
  }
  catch(...){
    cerr << "Fastjet error caught in selectJets" << endl;
    exit(20);
  }
}

// worker function to actually perform an analysis
void TimingAnalysis::FillTree(JetVector jets, JetVector TruthJets){

  double abstime;
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    jpt->push_back(jets[ijet].pt());
    jeta->push_back(jets[ijet].eta());
    jphi->push_back(jets[ijet].phi());
    jtime->push_back(ComputeTime(jets[ijet],abstime));
    jabstime->push_back(abstime);
    jtruth->push_back(TruthFrac(jets[ijet],TruthJets));
  }
  
  if(jets.size()>0)
    for (unsigned int icl=0; icl<jets[0].constituents().size(); icl++){    
      j0clpt->push_back(jets[0].constituents()[icl].pt());
      j0clphi->push_back(jets[0].constituents()[icl].phi());
      j0cleta->push_back(jets[0].constituents()[icl].eta());
      j0cltime->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().time());
      j0clabstime->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().abstime());
      j0cltruth->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().pileup() ? 0.0 : 1.0);
      j0clpixelID->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().pixel_id());
      j0clpixelNum->push_back(static_cast<double>(jets[0].constituents()[icl].user_info<TimingInfo>().pixel_num()));
      j0clpdgid->push_back(jets[0].constituents()[icl].user_info<TimingInfo>().pdg_id());
    }  
  if(jets.size()>1)
    for (unsigned int icl=0; icl<jets[1].constituents().size(); icl++){    
      j1clpt->push_back(jets[1].constituents()[icl].pt());
      j1clphi->push_back(jets[1].constituents()[icl].phi());
      j1cleta->push_back(jets[1].constituents()[icl].eta());
      j1clabstime->push_back(jets[1].constituents()[icl].user_info<TimingInfo>().abstime());
      j1cltruth->push_back(jets[1].constituents()[icl].user_info<TimingInfo>().pileup() ? 0.0 : 1.0);
    }  
}

void TimingAnalysis::FillTruthTree(JetVector jets){  
  double abstime;
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    truejpt->push_back(jets[ijet].pt());
    truejeta->push_back(jets[ijet].eta());
    truejphi->push_back(jets[ijet].phi());
    truejtime->push_back(ComputeTime(jets[ijet],abstime));
    truejabstime->push_back(abstime);
  }
}

double TimingAnalysis::ComputeTime(fastjet::PseudoJet jet, double &abstime){
  double time=0;
  float mindist = 1;
  float maxpt = 0;
  double pnum=0;

  abstime=0;
  double jetEta = jet.eta();
  double jetPhi = jet.phi();

  double pt,eta,phi,dist;
  double meanR=0.1;

  for (unsigned int i=0; i < jet.constituents().size(); i++){
    //if segmentation, only use timing from ghost particles    
    if(((not segmentation) or jet.constituents()[i].user_info<TimingInfo>().isGhost()) 
       and (abs(jet.constituents()[i].user_info<TimingInfo>().time()) < 100)){
      switch(timeMode){
      case highestPT:
	pt = jet.constituents()[i].user_info<TimingInfo>().pt();
	if(pt>maxpt){
	  maxpt = pt;
	  time = jet.constituents()[i].user_info<TimingInfo>().time();
	  abstime = jet.constituents()[i].user_info<TimingInfo>().abstime();
	}//endif
	break;
      case centralParticle:
	eta = jet.constituents()[i].eta();
	phi = jet.constituents()[i].phi();
	dist = sqrt(pow(eta-jetEta,2)+pow(phi-jetPhi,2));
	if(dist < mindist){
	  mindist = dist;
	  time = jet.constituents()[i].user_info<TimingInfo>().time();
	  abstime = jet.constituents()[i].user_info<TimingInfo>().abstime();
	}
	break;
      case mean:
	eta = jet.constituents()[i].eta();
        phi = jet.constituents()[i].phi();
        dist = sqrt(pow(eta-jetEta,2)+pow(phi-jetPhi,2));
        if(dist < meanR){
          time += jet.constituents()[i].user_info<TimingInfo>().time();
          abstime += jet.constituents()[i].user_info<TimingInfo>().abstime();
	  pnum++;
        }
	break;
      default:
	cerr << "ComputeTime called with invalid Timing Mode" << endl;
	return -999;
      }
    }
  }

  if(timeMode == mean){
    time/=pnum;
    abstime/=pnum;
  }

  return time;
}

double TimingAnalysis::TruthFrac(PseudoJet jet, JetVector truthJets){

  double maxFrac=0;

  //for each truth jet
  for (unsigned int tj = 0; tj < truthJets.size(); tj++){
    double ptTruthTot=0;
    double ptTot=0;
    //for each truth particle
    for (unsigned int ti=0; ti < truthJets[tj].constituents().size(); ti++){
      auto truthInfo = truthJets[tj].constituents()[ti].user_info<TimingInfo>();
      int truthID = truthInfo.pythia_id();
      ptTot += truthInfo.pt();
      
      //for each particle in the main jet
      for (unsigned int i=0; i < jet.constituents().size(); i++){
	auto pInfo = jet.constituents()[i].user_info<TimingInfo>();
	// if it is the identical particle, only use if we have timing info
	if ((not pInfo.pileup()) and (truthID == pInfo.pythia_id())){
	  //if using segmentation, only consider ghost particles
	  if(not pInfo.isGhost())
	    ptTruthTot += truthInfo.pt();
	}
      }
    }
    
    double frac=ptTruthTot/ptTot;
    if(frac > maxFrac)
      maxFrac=frac;
  }
  
  return maxFrac;
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
  tT->Branch("jabstime","std::vector<float>",&jabstime);
  tT->Branch("jtruth","std::vector<float>",&jtruth);

  tT->Branch("j0clpt","std::vector<float>",&j0clpt);
  tT->Branch("j0clphi","std::vector<float>",&j0clphi);
  tT->Branch("j0cleta","std::vector<float>",&j0cleta);
  tT->Branch("j0cltime","std::vector<float>",&j0cltime);
  tT->Branch("j0clabstime","std::vector<float>",&j0clabstime);
  tT->Branch("j0cltruth","std::vector<float>",&j0cltruth);
  tT->Branch("j0clpixelID","std::vector<float>",&j0clpixelID);
  tT->Branch("j0clpixelNum","std::vector<float>",&j0clpixelNum);
  tT->Branch("j0clpdgid","std::vector<float>",&j0clpdgid);

  tT->Branch("j1clpt","std::vector<float>",&j1clpt);
  tT->Branch("j1clphi","std::vector<float>",&j1clphi);
  tT->Branch("j1cleta","std::vector<float>",&j1cleta);
  tT->Branch("j1clabstime","std::vector<float>",&j1clabstime);
  tT->Branch("j1cltruth","std::vector<float>",&j1cltruth);

  tT->Branch("truejpt", "std::vector<float>",&truejpt);
  tT->Branch("truejphi","std::vector<float>",&truejphi);
  tT->Branch("truejeta","std::vector<float>",&truejeta);
  tT->Branch("truejtime","std::vector<float>",&truejtime);
  tT->Branch("truejabstime","std::vector<float>",&truejabstime);
  
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
      jabstime->clear();
      jtruth->clear();
      j0clpt->clear();
      j0clphi->clear();
      j0cleta->clear();
      j0cltime->clear();
      j0clabstime->clear();
      j0cltruth->clear();
      j0clpixelID->clear();
      j0clpixelNum->clear();
      j0clpdgid->clear();
      j1clpt->clear();
      j1clphi->clear();
      j1cleta->clear();
      j1clabstime->clear();
      j1cltruth->clear();
      truejpt->clear();
      truejphi->clear();
      truejeta->clear();
      truejtime->clear();
      truejabstime->clear();
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
