#include "TimingJetFinder.h"

TimingJetFinder::TimingJetFinder(double etaMin, double etaMax, double pixelSize, double R){
  
  _etaMin= (etaMin > 0) ? etaMin : 0;
  if(etaMax <= etaMin){
    cerr << "Invalid Eta Limits " << etaMin << " -> " << etaMax << "Passed to TimingJetFinder" << endl;
    exit(20);
  }
  _etaMax= etaMax;
  PixelSize(pixelSize);
  _R=R;

  const double grid_spacing(0.6);

  //suppress fastjet banner
  fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);

  jetDef.reset(new JetDefinition(fastjet::antikt_algorithm, _R, fastjet::E_scheme, fastjet::Best));
  active_area.reset(new AreaDefinition(fastjet::active_area));
  bge.reset(new GridMedianBackgroundEstimator(_etaMax, grid_spacing));
  select_fwd.reset(new Selector(SelectorAbsRapRange(_etaMin,_etaMax)));
}

void TimingJetFinder::PixelSize(double pixelSize){
  _pixelSize = (pixelSize > 1e-6) ? pixelSize : 1e-6;
}

void TimingJetFinder::selectJets(JetVector &particlesForJets, JetVector &selectedJets){

  try{
    fastjet::ClusterSequenceArea clustSeq(particlesForJets, (*jetDef), (*active_area));
    bge->set_particles(particlesForJets);

    try{
      fastjet::Subtractor subtractor(bge.get());
      
      JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
      JetVector subtractedJets = subtractor(inclusiveJets);
      
      selectedJets.clear();
      selectedJets = (*select_fwd)(subtractedJets);
    }
    catch(...){
      cerr << "Error caught" << endl;
      exit(20);
    }
  }
  catch(...){
    cerr << "Error caught here" << endl;
    exit(20);
  }
}

void TimingJetFinder::selectSegmentedJets(JetVector &particlesForJets, JetVector &selectedJets){
  selectJets(particlesForJets,selectedJets);
}
